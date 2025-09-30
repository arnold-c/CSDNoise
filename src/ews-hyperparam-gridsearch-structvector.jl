using Dates: Dates, Day
using ProgressMeter
using UnPack: @unpack
using FLoops: FLoops
using BangBang: BangBang
using MicroCollections: MicroCollections
using StyledStrings
using Try: Try
using DrWatson: @tagsave
using JLD2
using StructArrays: StructVector
using Bumper
using MultistartOptimization
using NLopt

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
        specification_vecs,
        data_arrs;
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
    n_total_scenarios = length(all_scenarios)

    if verbose
        println(styled"{green:Starting Grid Search with StructVector}")
        println(styled"Total grid points: {yellow:$n_total_scenarios}")
    end

    # Load existing results (including checkpoints)
    existing_results = if force
        StructVector(OptimizationResult[])
    else
        load_previous_gridsearch_results_structvector(filedir, gridsearch_filename_base)
    end

    n_existing = length(existing_results)
    if verbose && n_existing > 0
        println(styled"Found {cyan:$n_existing} existing results")
    end

    # Find missing scenarios - reuse function from multistart
    missing_scenarios = find_missing_scenarios(all_scenarios, existing_results)
    n_missing = length(missing_scenarios)

    if verbose
        println(styled"Missing grid points to evaluate: {yellow:$n_missing}")
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
        println(styled"{green:Grid search complete! Total results: {yellow:$n_final}}")
    end

    return return_results ? existing_results : nothing
end

"""
    create_gridsearch_scenarios_structvector(
    	specification_vecs;
		executor=ThreadedEx()
    )

Create all grid search scenarios including parameter combinations.
Returns StructVector{GridSearchScenario} with all combinations.
"""
function create_gridsearch_scenarios_structvector(specification_vecs)
    @unpack ensemble_specification_vec,
        null_specification_vec,
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec,
        ews_threshold_quantile_vec,
        ews_consecutive_thresholds_vec = specification_vecs

    combinations = Iterators.product(
        ensemble_specification_vec,
        null_specification_vec,
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec,
        ews_threshold_quantile_vec,
        ews_consecutive_thresholds_vec
    )
    n_combinations = length(combinations)

    # Create all combinations including grid parameters
    # Use explicit type annotation and avoid collect(...) splat pattern
    CombinationType = Tuple{
        EnsembleSpecification,
        EnsembleSpecification,
        NoiseSpecification,
        IndividualTestSpecification,
        Float64,
        EWSMetricSpecification,
        EWSEndDateType,
        EWSThresholdWindowType,
        Year,
        String,
        Float64,
        Int64,
    }

    scenarios_vec = Vector{GridSearchScenario}(undef, n_combinations)

    for (
            i, (
                ensemble_spec,
                null_spec,
                noise_spec,
                test_spec,
                percent_tested,
                ews_metric_spec,
                ews_enddate_type,
                ews_window,
                ews_burnin,
                ews_metric,
                threshold_quantile,
                consecutive_thresholds,
            ),
        ) in enumerate(combinations)

        scenarios_vec[i] = GridSearchScenario(
            ensemble_spec,
            null_spec,
            noise_spec,
            test_spec,
            percent_tested,
            ews_metric_spec,
            ews_enddate_type,
            ews_window,
            ews_burnin,
            ews_metric,
            threshold_quantile,
            consecutive_thresholds,
        )
    end

    return StructVector(scenarios_vec)
end


"""
    load_previous_gridsearch_results_structvector(filedir, filename_base)

Load previous grid search results as StructVector, including checkpoints.
Reuses checkpoint loading logic from multistart.
"""
function load_previous_gridsearch_results_structvector(filedir::String, filename_base::String)
    if !isdir(filedir)
        return StructVector(OptimizationResult[])
    end

    # Get most recent file - reuse function from multistart
    load_filepath = get_most_recent_hyperparam_filepath(filename_base, filedir)

    if Try.iserr(load_filepath)
        # Try to load checkpoint files
        return load_checkpoint_results_structvector(filedir)
    end

    try
        data = JLD2.load(Try.unwrap(load_filepath))

        # Try to load as StructVector first (new format)
        if haskey(data, "gridsearch_results")
            return data["gridsearch_results"]::StructVector{OptimizationResult}
        end

        # Fall back to DataFrame format (old format) and convert
        if haskey(data, "ews_df")
            df = data["ews_df"]
            return df_to_structvector(df, OptimizationResult)
        end

        return StructVector(OptimizationResult[])
    catch
        # Try to load checkpoint files as fallback
        return load_checkpoint_results_structvector(filedir)
    end
end

"""
    confirm_gridsearch_run_structvector(missing_scenarios; disable_time_check)

Confirm with user before running grid search on StructVector of scenarios.
"""
function confirm_gridsearch_run_structvector(
        missing_scenarios::StructVector{GridSearchScenario};
        disable_time_check::Bool = false
    )
    n_missing = length(missing_scenarios)

    if n_missing == 0
        @info "No missing grid points to evaluate"
        return false
    end

    if disable_time_check
        return true
    end

    # Estimate time (simpler than multistart since no optimization)
    estimated_time = n_missing * 0.5  # ~0.5 seconds per grid point

    if estimated_time > 300  # 5 minutes
        println(styled"{yellow:Warning:} Estimated grid search time: {red:$(round(estimated_time/60, digits=1))} minutes")
        print("Continue? (y/N): ")
        response = readline()
        return lowercase(strip(response)) in ["y", "yes"]
    end

    return true
end

"""
    evaluate_gridsearch_scenarios(
    	missing_scenarios,
    	data_arrs;
    	executor=FLoops.ThreadedEx(),
    	save_results=true,
		save_checkpoints=false,
    	checkpoint_dir="",
    	verbose=true
    )

Evaluate grid search scenarios with thread safety and incremental saving after all quantiles are run for a scenario.
Reuses the batching pattern from multistart optimization.
"""
function evaluate_gridsearch_scenarios(
        missing_scenarios::StructVector{GridSearchScenario},
        data_arrs::NamedTuple{T};
        executor = FLoops.ThreadedEx(),
        save_results = true,
        save_checkpoints = false,
        checkpoint_dir::String = "",
        verbose::Bool = true
    ) where {T}

    n_missing = length(missing_scenarios)

    if n_missing == 0
        return StructVector(OptimizationResult[])
    end

    all_results = OptimizationResult[]

    # Setup progress tracking
    if verbose
        prog = Progress(n_missing; desc = "Evaluating grid points: ", showspeed = true)
    end

    checkpoint_num = 0

    @unpack emergent_incidence_arr,
        null_incidence_arr = data_arrs
    # Creating the incidence, noise, testing arrays and EWS metrics is expensive
    # Loop through all unique combinations in a nested loop to allow re-use within
    # the optimization and grid search steps
    unique_ensemble_specs = unique(missing_scenarios.ensemble_specification)
    for ensemble_specification in unique_ensemble_specs
        ensemble_specific_scenarios = filter(
            s -> s.ensemble_specification == ensemble_specification,
            missing_scenarios
        )

        unique_null_specs = unique(ensemble_specific_scenarios.null_specification)
        for null_specification in unique_null_specs
            null_specific_scenarios = filter(
                s -> s.null_specification == null_specification,
                ensemble_specific_scenarios
            )
            unique_noise_specs = unique(null_specific_scenarios.noise_specification)
            for noise_specification in unique_noise_specs
                noise_specific_scenarios = filter(
                    s -> s.noise_specification == noise_specification,
                    null_specific_scenarios
                )

                noisearr = create_noise_arr(
                    noise_specification,
                    emergent_incidence_arr,
                    ensemble_specification;
                    seed = 1234,
                )[1]

                unique_perc_tested = unique(noise_specific_scenarios.percent_tested)
                for percent_tested in unique_perc_tested
                    perc_tested_specific_scenarios = filter(
                        s -> s.percent_tested == percent_tested,
                        noise_specific_scenarios
                    )

                    unique_test_specs = unique(perc_tested_specific_scenarios.test_specification)
                    for test_specification in unique_test_specs
                        test_specific_scenarios = filter(
                            s -> s.test_specification == test_specification,
                            perc_tested_specific_scenarios
                        )

                        testarr = create_testing_arrs(
                            emergent_incidence_arr,
                            noisearr,
                            percent_tested,
                            test_specification,
                        )

                        null_testarr = create_testing_arrs(
                            null_incidence_arr,
                            noisearr,
                            percent_tested,
                            test_specification,
                        )

                        unique_enddate_types = unique(test_specific_scenarios.ews_enddate_type)
                        for ews_enddate_type in unique_enddate_types
                            ews_enddate_specific_scenarios = filter(
                                s -> s.ews_enddate_type == ews_enddate_type,
                                test_specific_scenarios
                            )

                            unique_ews_metric_specs = unique(ews_enddate_specific_scenarios.ews_metric_specification)

                            for ews_metric_specification in unique_ews_metric_specs
                                metric_spec_specific_scenarios = filter(
                                    s -> s.ews_metric_specification == ews_metric_specification,
                                    ews_enddate_specific_scenarios
                                )

                                ews_metrics, null_ews_metrics = generate_ensemble_ews_metrics(
                                    data_arrs,
                                    testarr,
                                    null_testarr,
                                    ews_metric_specification,
                                    ews_enddate_type,
                                )

                                isempty(ews_metrics) && error("No valid EWS metrics could be computed for any emergent simulation")
                                isempty(null_ews_metrics) && error("No valid EWS metrics could be computed for any null simulation")

                                unique_ews_metrics = unique(metric_spec_specific_scenarios.ews_metric)
                                for ews_metric in unique_ews_metrics
                                    ews_metric_specific_scenarios = filter(
                                        s -> s.ews_metric == ews_metric,
                                        ews_enddate_specific_scenarios
                                    )

                                    unique_threshold_windows = unique(ews_metric_specific_scenarios.ews_threshold_window)
                                    for ews_threshold_window in unique_threshold_windows
                                        window_specific_scenarios = filter(
                                            s -> s.ews_threshold_window == ews_threshold_window,
                                            ews_metric_specific_scenarios
                                        )

                                        unique_threshold_burnin = unique(window_specific_scenarios.ews_threshold_burnin)
                                        for ews_threshold_burnin in unique_threshold_burnin
                                            burnin_specific_scenarios = filter(
                                                s -> s.ews_threshold_burnin == ews_threshold_burnin,
                                                window_specific_scenarios
                                            )

                                            unique_consecutive_thresholds = unique(burnin_specific_scenarios.consecutive_thresholds)
                                            for consecutive_thresholds in unique_consecutive_thresholds
                                                consecutive_thresh_specific_scenarios = filter(
                                                    s -> s.consecutive_thresholds == consecutive_thresholds,
                                                    burnin_specific_scenarios
                                                )

                                                unique_threshold_quantiles = unique(consecutive_thresh_specific_scenarios.threshold_quantile)

                                                optimization_results = Vector{OptimizationResult}(undef, length(unique_threshold_quantiles))
                                                FLoops.@floop executor for (i, threshold_quantile) in pairs(unique_threshold_quantiles)
                                                    local grid_scenario = GridSearchScenario(
                                                        ensemble_specification,
                                                        null_specification,
                                                        noise_specification,
                                                        test_specification,
                                                        percent_tested,
                                                        ews_metric_specification,
                                                        ews_enddate_type,
                                                        ews_threshold_window,
                                                        ews_threshold_burnin,
                                                        ews_metric,
                                                        threshold_quantile,
                                                        consecutive_thresholds
                                                    )

                                                    optimization_results[i] = gridsearch_optimization_function(
                                                        grid_scenario,
                                                        ews_metrics,
                                                        null_ews_metrics
                                                    )
                                                end

                                                BangBang.append!!(all_results, optimization_results)

                                                # Save checkpoint after each batch
                                                if save_results && save_checkpoints && !isempty(checkpoint_dir)
                                                    save_checkpoint_structvector(
                                                        StructVector(all_results),
                                                        checkpoint_dir,
                                                        checkpoint_num
                                                    )
                                                    verbose && @info "Saved checkpoint after batch $checkpoint_num"
                                                    checkpoint_num += 1
                                                end

                                                if verbose
                                                    next!(prog)
                                                end


                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return StructVector(all_results)
end
