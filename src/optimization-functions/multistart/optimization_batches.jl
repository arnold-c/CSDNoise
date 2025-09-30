using StructArrays: StructVector
using FLoops: FLoops
using BangBang: BangBang
using DataFrames: DataFrames
using ProgressMeter: ProgressMeter

export optimize_scenarios_in_batches_structvector

"""
    optimize_scenarios_in_batches_structvector(
        missing_scenarios, data_arrs, bounds, optim_config;
        batch_size=nothing, executor=FLoops.SequentialEx(), save_results=true, checkpoint_dir="", verbose=true
    )

Optimize scenarios in batches with thread safety and incremental saving.

# Behavior:
- If `batch_size` is specified: Process scenarios in batches and save a checkpoint after each batch
- If `batch_size` is `nothing`: Run all simulations multithreaded without batching or checkpoints

# Thread Safety Design
This function follows a fork-join pattern with clear phase separation:

1. **Parallel Phase**: Each batch is processed with FLoops.@floop where:
   - Each task writes to a unique, pre-allocated index in batch_results[idx]
   - No shared mutable state between tasks
   - No use of threadid() or thread-local storage

2. **Serial Phase**: Between batches, all I/O and state management happens single-threaded:
   - DataFrame creation and merging
   - Progress updates
   - Checkpoint saving
   - Result accumulation

This avoids the threadid() anti-pattern described in:
https://julialang.org/blog/2023/07/PSA-dont-use-threadid/

The design is robust to task migration, thread adoption, and other Julia threading
implementation details because it reasons about tasks, not threads.
"""
function optimize_scenarios_in_batches_structvector(
        missing_scenarios::StructVector{OptimizationScenario},
        data_arrs::T1,
        bounds::T2,
        optim_config::T3;
        executor = FLoops.SequentialEx(),
        batch_size::Union{Int64, Nothing} = nothing,
        save_results = true,
        checkpoint_dir::String = "",
        verbose::Bool = true
    ) where {T1 <: NamedTuple, T2 <: NamedTuple, T3 <: NamedTuple}

    n_missing = length(missing_scenarios)

    if n_missing == 0
        return StructVector(OptimizationResult[])
    end

    # If no batch_size specified, run all scenarios multithreaded without checkpointing
    if isnothing(batch_size)
        verbose && @info "Running all $n_missing scenarios without batching"

        # Setup progress tracking
        if verbose
            prog = ProgressMeter.Progress(n_missing; desc = "Optimizing scenarios: ", showspeed = true)
        end

        # Process all scenarios in parallel
        FLoops.@floop executor for scenario in missing_scenarios
            FLoops.@reduce all_results = vcat(
                OptimizationResult[],
                optimize_single_scenario(
                    scenario,
                    data_arrs,
                    bounds,
                    optim_config
                )
            )
        end

        verbose && ProgressMeter.finish!(prog)
        return StructVector(all_results)
    end

    # Batched execution with checkpointing
    verbose && @info "Running $n_missing scenarios in batches of $batch_size with checkpointing"

    # Ensure batch size is reasonable for threading
    if executor isa FLoops.ThreadedEx # Default value
        # Choose a bigger batch size than the number of threads so Julia can
        # allocate tasks to threads appropriately and can use at least as many
        # tasks as there are threads
        batch_size < Threads.nthreads() && @warn "Consider setting the batch size to at least the number of threads ($(Threads.nthreads()))"
        batch_size % Threads.nthreads() == 0 && @warn "Consider setting the batch size to a multiple of the number of threads ($(Threads.nthreads()))"
    end

    # Initialize results storage
    all_results = OptimizationResult[]

    # Setup progress tracking
    if verbose
        prog = ProgressMeter.Progress(n_missing; desc = "Optimizing scenarios: ", showspeed = true)
    end

    # Process scenarios in batches
    scenario_batches = collect(Iterators.partition(1:n_missing, batch_size))

    for (batch_idx, batch_indices) in pairs(scenario_batches)
        batch_scenarios = @view missing_scenarios[batch_indices]
        batch_size_actual = length(batch_indices)

        verbose && println(styled"{green:Processing batch $batch_idx/$(length(scenario_batches)) ($(batch_size_actual) scenarios)}")

        FLoops.@floop executor for scenario in batch_scenarios
            # Direct struct access - no conversion needed!
            FLoops.@reduce batch_results = BangBang.append!!(
                OptimizationResult[],
                [
                    optimize_single_scenario(
                        scenario,
                        data_arrs,
                        bounds,
                        optim_config
                    ),
                ]
            )
        end
        BangBang.append!!(all_results, batch_results)

        # Update progress
        verbose && ProgressMeter.update!(prog, sum(length.(scenario_batches[1:batch_idx])))

        # Save checkpoint after each batch when batch_size is specified
        if save_results && !isempty(checkpoint_dir)
            save_checkpoint_structvector(
                StructVector(all_results),
                checkpoint_dir,
                batch_idx
            )
            verbose && @info "Saved checkpoint after batch $batch_idx"
        end
    end

    return StructVector(all_results)
end

# Keep old function for backward compatibility
function optimize_scenarios_in_batches(
        missing_scenarios_df::DataFrame,
        data_arrs::T1,
        bounds::T2,
        optim_config::T3;
        batch_size::Int = 10,
        executor = FLoops.SequentialEx(),
        save_every_n::Int = 10,
        save_results = true,
        checkpoint_dir::String = "",
        verbose::Bool = true
    ) where {T1 <: NamedTuple, T2 <: NamedTuple, T3 <: NamedTuple}
    n_missing = nrow(missing_scenarios_df)

    if n_missing == 0
        return DataFrame()
    end

    # Auto-configure batch size for threaded execution
    # Note: We use nthreads() only as a heuristic, not for correctness
    if executor isa FLoops.ThreadedEx
        if batch_size == 10  # Default value
            # Use nthreads() as a heuristic for load balancing, but cap at reasonable bounds
            suggested_batch_size = max(1, n_missing ÷ (Threads.nthreads() * 4))
            batch_size = clamp(suggested_batch_size, 1, 100)  # Cap between 1 and 100
            verbose && @info "Auto-configured batch size for threading: $batch_size (suggested: $suggested_batch_size)"
        end
    end

    # Initialize results storage
    all_results = DataFrame[]

    # Setup progress tracking
    if verbose
        prog = Progress(n_missing; desc = "Optimizing scenarios: ", showspeed = true)
    end

    # Process scenarios in batches
    scenario_batches = collect(Iterators.partition(1:n_missing, batch_size))

    for (batch_idx, batch_indices) in pairs(scenario_batches)
        batch_scenarios = missing_scenarios_df[batch_indices, :]
        batch_size_actual = length(batch_indices)

        verbose && println(styled"{green:Processing batch $batch_idx/$(length(scenario_batches)) ($(batch_size_actual) scenarios)}")

        # Process batch in parallel (thread-safe within batch)
        # THREAD SAFETY: Pre-allocate fixed-size array where each task writes to unique index
        # This avoids the threadid() anti-pattern and follows task-based parallelism principles
        batch_results = Vector{OptimizedValues}(undef, batch_size_actual)

        FLoops.@floop executor for (idx, scenario_row) in pairs(eachrow(batch_scenarios))
            # Convert DataFrame row to NamedTuple
            scenario = dataframe_row_to_scenario(scenario_row)

            # Run optimization (each task works independently)
            result = optimize_single_scenario(
                scenario,
                data_arrs,
                bounds,
                optim_config
            )

            # Store result at unique index (thread-safe: no overlapping writes)
            batch_results[idx] = result
        end

        # SERIAL PHASE: Convert batch results to DataFrame (single-threaded)
        # All parallel tasks have completed at this point (implicit barrier)
        batch_scenarios_vec = StructVector([dataframe_row_to_scenario(row) for row in eachrow(batch_scenarios)])
        batch_df = create_results_dataframe(
            batch_scenarios_vec,
            StructVector(batch_results),
        )

        # Add to results (single-threaded - no race conditions)
        push!(all_results, batch_df)

        # Update progress (single-threaded)
        verbose && update!(prog, sum(length.(scenario_batches[1:batch_idx])))

        # Save checkpoint periodically (single-threaded I/O)
        if batch_idx % save_every_n == 0 && !isempty(checkpoint_dir) && save_results
            combined_df = vcat(all_results...; cols = :union)
            save_checkpoint_atomic(combined_df, checkpoint_dir, batch_idx)
        end
    end

    # Combine all results
    final_results_df = if !isempty(all_results)
        vcat(all_results...; cols = :union)
    else
        create_empty_results_dataframe()
    end

    return final_results_df
end
