export ews_multistart_optimization

"""
    ews_multistart_optimization(specification_vecs, data_arrs; kwargs...)

Optimize EWS hyperparameters using multistart optimization with Sobol sequences.
Alternative to grid search that scales better with parameter dimensionality.

# Arguments
- `specification_vecs`: Named tuple containing specification vectors
- `data_arrs`: Named tuple containing ensemble data arrays

# Keyword Arguments
- `filedir`: Directory for saving results
- `optimization_filename_base`: Base filename for results
- `n_sobol_points`: Number of Sobol sequence points per scenario (default: 100)
- `local_algorithm`: NLopt local optimization algorithm (default: NLopt.LN_BOBYQA)
- `maxeval`: Maximum evaluations per local optimization (default: 1000)
- `xtol_rel`: Relative tolerance for parameter convergence (default: 1e-3)
- `xtol_abs`: Absolute tolerance for parameter convergence (default: 1e-3)
- `ftol_rel`: Relative tolerance for function convergence (default: 1e-4)
- `executor`: Executor for loops (default: `FLoops.SequentialEx()`)
- `quantile_bounds`: Bounds for threshold quantile (default: (0.5, 0.99))
- `consecutive_bounds`: Bounds for consecutive thresholds (default: (1.0, 10.0))
- `force`: Force recomputation even if results exist (default: false)
- `return_df`: Return DataFrame of results (default: true)
- `save_results`: Save results to file (default: true)
- `verbose`: Print progress information (default: true)
- `batch_size`: Batch size for processing scenarios. If specified, saves checkpoint after each batch. If `nothing`, runs all multithreaded (default: 10)
- `disable_time_check`: Skip time estimation confirmation (default: false)

# Returns
- DataFrame with optimal parameters for each scenario
"""
function ews_multistart_optimization(
        specification_vecs,
        data_arrs;
        # File management
        filedir = outdir("ensemble", "ews-multistart-optimization"),
        optimization_filename_base = "ews-multistart-optimization.jld2",
        optimization_output_filepath = joinpath(
            filedir,
            string(Dates.now()) * "_" * optimization_filename_base,
        ),
        # Optimization configuration
        n_sobol_points = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
        executor = FLoops.SequentialEx(),
        # Parameter bounds
        quantile_bounds = (0.5, 0.99),
        consecutive_bounds = (1.0, 10.0),
        # Control options
        force = false,
        return_df = true,
        save_results = true,
        verbose = true,
        # Caching options
        batch_size = 10,
        disable_time_check = false,
    )
    if !isdir(filedir)
        mkpath(filedir)
    end

    # Setup checkpoint directory
    checkpoint_dir = joinpath(filedir, "checkpoints")

    # Create all scenario combinations as StructVector
    all_scenarios = create_scenarios_structvector(specification_vecs)

    if verbose
        n_total_scenarios = length(all_scenarios)
        println(StyledStrings.styled"{green:Starting Multistart Optimization with Caching}")
        println(StyledStrings.styled"Total scenarios defined: {yellow:$n_total_scenarios}")
        println(StyledStrings.styled"Sobol points per scenario: {blue:$n_sobol_points}")
    end

    # Load existing results (including checkpoints)
    existing_results = if force
        StructVector(OptimizationResult[])
    else
        load_previous_multistart_results_structvector(filedir, optimization_filename_base)
    end

    n_existing = length(existing_results)
    if verbose && n_existing > 0
        println(StyledStrings.styled"Found {cyan:$n_existing} existing results")
    end

    # Find missing scenarios
    missing_scenarios = find_missing_scenarios(all_scenarios, existing_results)
    n_missing = length(missing_scenarios)

    if verbose
        println(StyledStrings.styled"Missing scenarios to optimize: {yellow:$n_missing}")
    end

    # Check with user if needed
    if n_missing > 0
        if !confirm_optimization_run_structvector(
                missing_scenarios,
                n_sobol_points;
                disable_time_check = disable_time_check
            )
            @info "Optimization cancelled by user"
            return return_df ? DF.DataFrame(existing_results) : existing_results
        end
    else
        @info "All scenarios already optimized"
        return return_df ? DF.DataFrame(existing_results) : existing_results
    end

    # Define parameter bounds (only quantile and consecutive thresholds)
    bounds = (
        lowers = [quantile_bounds[1], consecutive_bounds[1]],
        uppers = [quantile_bounds[2], consecutive_bounds[2]],
    )

    # Optimization configuration
    optim_config = (
        n_sobol_points = n_sobol_points,
        local_algorithm = local_algorithm,
        maxeval = maxeval,
        xtol_rel = xtol_rel,
        xtol_abs = xtol_abs,
        ftol_rel = ftol_rel,
    )

    # Run optimization for missing scenarios in batches
    new_results = optimize_scenarios_in_batches_structvector(
        missing_scenarios,
        data_arrs,
        bounds,
        optim_config;
        batch_size = batch_size,
        executor = executor,
        save_results = save_results,
        checkpoint_dir = checkpoint_dir,
        verbose = verbose
    )

    # Combine existing and new results
    final_results = if !isempty(existing_results)
        vcat(existing_results, new_results)
    else
        new_results
    end

    # Save final results
    if save_results
        save_results_structvector(final_results, optimization_output_filepath)

        # Clean up checkpoints after successful save
        cleanup_checkpoints(checkpoint_dir)
    end

    if verbose
        n_final = length(final_results)
        println(StyledStrings.styled"{green:Optimization complete! Total results: {yellow:$n_final}}")
    end

    return return_df ? DF.DataFrame(final_results) : final_results
end
