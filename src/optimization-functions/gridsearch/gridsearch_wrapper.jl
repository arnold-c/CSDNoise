export ews_hyperparam_gridsearch_structvector

"""
    ews_hyperparam_gridsearch_structvector(specification_vecs, data_arrs; kwargs...)

Grid search optimization using StructVector for efficient parallel processing.
Reuses functions from multistart optimization where possible.

# Arguments
- `specification_vecs`: Named tuple containing specification vectors including grid parameters
- `data_arrs`: Named tuple containing ensemble data arrays

# Keyword Arguments
- `filedir`: Directory for saving results
- `gridsearch_filename_base`: Base filename for results
- `executor`: Executor for loops (default: `FLoops.ThreadedEx()`)
- `batch_size`: Batch size for processing scenarios (default: 50)
- `force`: Force recomputation even if results exist (default: false)
- `return_results` Return the results (default: true)
- `save_results`: Save results to file (default: true)
- `verbose`: Print progress information (default: true)
- `disable_time_check`: Skip time estimation confirmation (default: false)
"""
function ews_hyperparam_gridsearch_structvector(
        specification_vecs::GridSearchSpecificationVecs;
        # File management
        filedir = outdir("ensemble", "ews-hyperparam-gridsearch"),
        gridsearch_filename_base = "ews-hyperparam-gridsearch-structvector.jld2",
        gridsearch_output_filepath = joinpath(
            filedir,
            string(Dates.now()) * "_" * gridsearch_filename_base,
        ),
        # Execution configuration
        executor = FLoops.ThreadedEx(),
        # Control options
        force = false,
        return_results = true,
        save_results = true,
        save_checkpoints = false,
        verbose = true,
        disable_time_check = false,
    )
    if !isdir(filedir)
        mkpath(filedir)
    end

    # Setup checkpoint directory
    checkpoint_dir = joinpath(filedir, "checkpoints")

    # Create all grid search scenarios as StructVector
    # This includes all combinations with grid parameters
    all_scenarios = create_gridsearch_scenarios_structvector(specification_vecs)

    if verbose
        n_total_scenarios = length(all_scenarios)
        println(StyledStrings.styled"{green:Starting Grid Search with StructVector}")
        println(StyledStrings.styled"Total grid points: {yellow:$n_total_scenarios}")
    end

    # Load existing results (including checkpoints)
    existing_results = if force
        StructVector(OptimizationResult[])
    else
        load_previous_gridsearch_results_structvector(filedir, gridsearch_filename_base)
    end

    n_existing = length(existing_results)
    if verbose && n_existing > 0
        println(StyledStrings.styled"Found {cyan:$n_existing} existing results")
    end

    # Find missing scenarios - reuse function from multistart
    missing_scenarios = find_missing_scenarios(all_scenarios, existing_results)
    n_missing = length(missing_scenarios)

    if verbose
        println(StyledStrings.styled"Missing grid points to evaluate: {yellow:$n_missing}")
    end

    # Check with user if needed
    if n_missing > 0
        if !confirm_gridsearch_run_structvector(
                missing_scenarios;
                disable_time_check = disable_time_check
            )
            @info "Grid search cancelled by user"
            return return_results ? existing_results : nothing
        end
    else
        @info "All grid points already evaluated"
        return return_results ? existing_results : nothing
    end

    # Run grid search for missing scenarios in batches
    new_results = evaluate_gridsearch_scenarios(
        missing_scenarios,
        data_arrs;
        executor = executor,
        save_results = save_results,
        save_checkpoints = save_checkpoints,
        checkpoint_dir = checkpoint_dir,
        verbose = verbose
    )

    # Combine existing and new results
    if !isempty(existing_results)
        BangBang.append!!(existing_results, new_results)
    end

    # Save final results - reuse function from multistart
    if save_results
        save_results_structvector(existing_results, gridsearch_output_filepath)

        # Clean up checkpoints after successful save
        cleanup_checkpoints(checkpoint_dir)
    end

    if verbose
        n_final = length(existing_results)
        println(StyledStrings.styled"{green:Grid search complete! Total results: {yellow:$n_final}}")
    end

    return return_results ? existing_results : nothing
end
