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

# Computation key structs for hierarchical grouping to minimize allocations
struct NoiseComputationKey
    ensemble_specification
    null_specification
    noise_specification
end

struct TestingComputationKey
    noise_key::NoiseComputationKey
    percent_tested
    test_specification
end

struct EWSComputationKey
    testing_key::TestingComputationKey
    ews_metric_specification
    ews_enddate_type
end

struct ScenarioComputationKey
    ews_key::EWSComputationKey
    ews_metric
    ews_threshold_window
    ews_threshold_burnin
    consecutive_thresholds
end

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

function evaluate_gridsearch_scenarios2(
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
    unique_null_specs = unique(missing_scenarios.null_specification)
    unique_noise_specs = unique(missing_scenarios.noise_specification)
    unique_perc_tested = unique(missing_scenarios.percent_tested)
    unique_test_specs = unique(missing_scenarios.test_specification)
    unique_enddate_types = unique(missing_scenarios.ews_enddate_type)
    unique_ews_metric_specs = unique(missing_scenarios.ews_metric_specification)
    unique_ews_metrics = unique(missing_scenarios.ews_metric)
    unique_threshold_windows = unique(missing_scenarios.ews_threshold_window)
    unique_threshold_burnin = unique(missing_scenarios.ews_threshold_burnin)
    unique_consecutive_thresholds = unique(missing_scenarios.consecutive_thresholds)
    unique_threshold_quantiles = unique(missing_scenarios.threshold_quantile)

    for ensemble_specification in unique_ensemble_specs
        for null_specification in unique_null_specs
            for noise_specification in unique_noise_specs

                noisearr = create_noise_arr(
                    noise_specification,
                    emergent_incidence_arr,
                    ensemble_specification;
                    seed = 1234,
                )[1]

                for percent_tested in unique_perc_tested
                    for test_specification in unique_test_specs

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

                        for ews_enddate_type in unique_enddate_types
                            for ews_metric_specification in unique_ews_metric_specs

                                ews_metrics, null_ews_metrics = generate_ensemble_ews_metrics(
                                    data_arrs,
                                    testarr,
                                    null_testarr,
                                    ews_metric_specification,
                                    ews_enddate_type,
                                )

                                isempty(ews_metrics) && error("No valid EWS metrics could be computed for any emergent simulation")
                                isempty(null_ews_metrics) && error("No valid EWS metrics could be computed for any null simulation")

                                n_emergent_sims = length(ews_metrics)
                                n_null_sims = length(null_ews_metrics)

                                for ews_metric in unique_ews_metrics
                                    for ews_threshold_window in unique_threshold_windows
                                        for ews_threshold_burnin in unique_threshold_burnin
                                            for consecutive_thresholds in unique_consecutive_thresholds
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

function evaluate_gridsearch_scenarios_optimized(
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

    all_results = Vector{OptimizationResult}()
    sizehint!(all_results, n_missing)

    if verbose
        prog = Progress(n_missing; desc = "Evaluating grid points (optimized): ", showspeed = true)
    end

    checkpoint_num = 0
    processed_scenarios = 0

    @unpack emergent_incidence_arr, null_incidence_arr = data_arrs

    noise_groups = group_by_noise_key(missing_scenarios)

    for (noise_key, noise_scenarios) in noise_groups
        noisearr = create_noise_arr(
            noise_key.noise_specification,
            emergent_incidence_arr,
            noise_key.ensemble_specification;
            seed = 1234,
        )[1]

        testing_groups = group_by_testing_key(noise_scenarios)

        for (testing_key, testing_scenarios) in testing_groups
            testarr = create_testing_arrs(
                emergent_incidence_arr,
                noisearr,
                testing_key.percent_tested,
                testing_key.test_specification,
            )

            null_testarr = create_testing_arrs(
                null_incidence_arr,
                noisearr,
                testing_key.percent_tested,
                testing_key.test_specification,
            )

            ews_groups = group_by_ews_key(testing_scenarios)

            for (ews_key, ews_scenarios) in ews_groups
                ews_metrics, null_ews_metrics = generate_ensemble_ews_metrics(
                    data_arrs,
                    testarr,
                    null_testarr,
                    ews_key.ews_metric_specification,
                    ews_key.ews_enddate_type,
                )

                isempty(ews_metrics) && error("No valid EWS metrics could be computed for any emergent simulation")
                isempty(null_ews_metrics) && error("No valid EWS metrics could be computed for any null simulation")

                scenario_groups = group_by_scenario_key(ews_scenarios)

                for (scenario_key, final_scenarios) in scenario_groups
                    unique_quantiles = unique(s.threshold_quantile for s in final_scenarios)
                    group_results = Vector{OptimizationResult}(undef, length(unique_quantiles))

                    FLoops.@floop executor for (i, threshold_quantile) in pairs(unique_quantiles)
                        grid_scenario = GridSearchScenario(
                            scenario_key.ews_key.testing_key.noise_key.ensemble_specification,
                            scenario_key.ews_key.testing_key.noise_key.null_specification,
                            scenario_key.ews_key.testing_key.noise_key.noise_specification,
                            scenario_key.ews_key.testing_key.test_specification,
                            scenario_key.ews_key.testing_key.percent_tested,
                            scenario_key.ews_key.ews_metric_specification,
                            scenario_key.ews_key.ews_enddate_type,
                            scenario_key.ews_threshold_window,
                            scenario_key.ews_threshold_burnin,
                            scenario_key.ews_metric,
                            threshold_quantile,
                            scenario_key.consecutive_thresholds
                        )

                        group_results[i] = gridsearch_optimization_function(
                            grid_scenario,
                            ews_metrics,
                            null_ews_metrics
                        )
                    end

                    append!(all_results, group_results)
                    processed_scenarios += length(group_results)

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
                        for _ in 1:length(group_results)
                            next!(prog)
                        end
                    end
                end
            end
        end
    end

    return StructVector(all_results)
end

"""
    evaluate_gridsearch_scenarios_multistart(
        missing_scenarios,
        data_arrs;
        executor=FLoops.ThreadedEx(),
        save_results=true,
        save_checkpoints=false,
        checkpoint_dir="",
        verbose=true,
        n_sobol_points=100,
        local_algorithm=NLopt.LN_BOBYQA,
        maxeval=1000,
        xtol_rel=1.0e-3,
        xtol_abs=1.0e-3,
        ftol_rel=1.0e-4,
        quantile_bounds=(0.5, 0.99)
    )

Evaluate grid search scenarios with multistart optimization for quantile thresholds.
Keeps the nested loop structure but replaces the inner loop over quantiles with optimization.
Reuses multistart optimization code for the optimization logic.
"""
function evaluate_gridsearch_scenarios_multistart(
        missing_scenarios::StructVector{GridSearchScenario},
        data_arrs::NamedTuple{T};
        executor = FLoops.ThreadedEx(),
        save_results = true,
        save_checkpoints = false,
        checkpoint_dir::String = "",
        verbose::Bool = true,
        # Optimization configuration
        n_sobol_points = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
        quantile_bounds = (0.5, 0.99)
    ) where {T}

    n_missing = length(missing_scenarios)

    if n_missing == 0
        return StructVector(OptimizationResult[])
    end

    all_results = OptimizationResult[]

    # Setup progress tracking
    if verbose
        prog = Progress(n_missing; desc = "Optimizing scenarios with multistart: ", showspeed = true)
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

                                n_emergent_sims = length(ews_metrics)
                                n_null_sims = length(null_ews_metrics)

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

                                                # MULTISTART OPTIMIZATION REPLACES THE INNER LOOP
                                                # Instead of looping over unique_threshold_quantiles, optimize them
                                                optimization_result = optimize_quantile_threshold_multistart(
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
                                                    consecutive_thresholds,
                                                    ews_metrics,
                                                    null_ews_metrics;
                                                    n_sobol_points = n_sobol_points,
                                                    local_algorithm = local_algorithm,
                                                    maxeval = maxeval,
                                                    xtol_rel = xtol_rel,
                                                    xtol_abs = xtol_abs,
                                                    ftol_rel = ftol_rel,
                                                    quantile_bounds = quantile_bounds
                                                )

                                                push!(all_results, optimization_result)

                                                # Save checkpoint after each optimization
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

"""
    optimize_quantile_threshold_multistart(
        ensemble_specification, null_specification, noise_specification,
        test_specification, percent_tested, ews_metric_specification,
        ews_enddate_type, ews_threshold_window, ews_threshold_burnin,
        ews_metric, consecutive_thresholds, ews_metrics, null_ews_metrics;
        n_sobol_points=100, local_algorithm=NLopt.LN_BOBYQA, maxeval=1000,
        xtol_rel=1.0e-3, xtol_abs=1.0e-3, ftol_rel=1.0e-4,
        quantile_bounds=(0.5, 0.99)
    )

Optimize the quantile threshold for a single scenario using multistart optimization.
Reuses the multistart optimization logic from ews-multistart-optimization.jl.
"""
function optimize_quantile_threshold_multistart(
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
        consecutive_thresholds,
        ews_metrics,
        null_ews_metrics;
        n_sobol_points = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
        quantile_bounds = (0.5, 0.99)
    )

    # Create optimization scenario (similar to OptimizationScenario but for grid search)
    scenario = GridSearchScenario(
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
        0.5,  # placeholder quantile - will be optimized
        consecutive_thresholds
    )

    # Create tracker instance for this scenario
    tracker = OptimizationTracker()

    # Create objective function closure that updates tracker
    objective = quantile_vec -> quantile_only_objective_function_with_tracking(
        quantile_vec,
        scenario,
        consecutive_thresholds,
        ews_metrics,
        null_ews_metrics,
        tracker
    )

    # Setup multistart problem (only optimizing quantile, so 1D)
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [quantile_bounds[1]],  # lower bound
        [quantile_bounds[2]]   # upper bound
    )

    # Configure local optimization method
    local_method = MultistartOptimization.NLopt_local_method(
        local_algorithm;
        xtol_rel = xtol_rel,
        xtol_abs = xtol_abs,
        maxeval = maxeval,
    )

    # Configure multistart method (TikTak uses Sobol sequences)
    multistart_method = MultistartOptimization.TikTak(n_sobol_points)

    # Run optimization
    MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    # Extract optimal parameters from tracker (which has the best metrics)
    optimal_quantile = tracker.best_params[1]

    return OptimizationResult(
        # From scenario
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
        # From optimized values
        optimal_quantile,
        consecutive_thresholds,
        tracker.best_accuracy,
        tracker.best_sensitivity,
        tracker.best_specificity
    )
end


"""
    quantile_only_objective_function_with_tracking(
        quantile_vec, scenario, consecutive_thresholds, ews_metrics, null_ews_metrics, tracker
    )

Objective function for quantile threshold optimization that tracks the best solution.
Returns 1 - accuracy for minimization and updates tracker with best metrics.
Only optimizes the quantile threshold, keeping consecutive thresholds fixed.
Adapted from ews_objective_function_with_tracking in ews-multistart-optimization.jl.
"""
function quantile_only_objective_function_with_tracking(
        quantile_vec::Vector{Float64},
        scenario::GridSearchScenario,
        consecutive_thresholds::Int64,
        ews_metrics,
        null_ews_metrics,
        tracker::OptimizationTracker
    )

    # Extract quantile from vector (1D optimization)
    threshold_quantile = quantile_vec[1]

    # Pre-compute values outside the loop for type stability
    ews_metric_symbol = Symbol(scenario.ews_metric)

    # Calculate accuracy
    n_emergent_sims = length(ews_metrics)
    n_null_sims = length(null_ews_metrics)
    true_positives = 0
    true_negatives = 0

    for sim in eachindex(ews_metrics)
        # Use pre-computed EWS metrics
        ews_vals = ews_metrics[sim]
        null_ews_vals = null_ews_metrics[sim]

        # Check threshold exceedances
        exceeds_threshold = exceeds_ews_threshold(
            ews_vals,
            ews_metric_symbol,
            scenario.ews_threshold_window,
            threshold_quantile,
            scenario.ews_threshold_burnin,
        )

        detection_index = calculate_ews_trigger_index(
            exceeds_threshold,
            consecutive_thresholds,
        )

        null_exceeds_threshold = exceeds_ews_threshold(
            null_ews_vals,
            ews_metric_symbol,
            scenario.ews_threshold_window,
            threshold_quantile,
            scenario.ews_threshold_burnin,
        )

        null_detection_index = calculate_ews_trigger_index(
            null_exceeds_threshold,
            consecutive_thresholds,
        )

        # Update counts
        if Try.isok(detection_index)
            true_positives += 1
        end
        if Try.iserr(null_detection_index)
            true_negatives += 1
        end
    end

    # Calculate metrics
    sensitivity = calculate_sensitivity(true_positives, n_emergent_sims)
    specificity = calculate_specificity(true_negatives, n_null_sims)
    accuracy = calculate_balanced_accuracy(sensitivity, specificity)
    loss = 1.0 - accuracy

    # Update tracker if this is the best solution so far
    if loss < tracker.best_loss
        tracker.best_loss = loss
        tracker.best_accuracy = accuracy
        tracker.best_sensitivity = sensitivity
        tracker.best_specificity = specificity
        tracker.best_params = copy(quantile_vec)
    end

    return loss  # Return scalar loss for optimizer
end

function evaluate_gridsearch_scenarios_bumper(
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

    # Pre-allocate results
    all_results = Vector{OptimizationResult}(undef, n_missing)
    result_idx = 0

    if verbose
        prog = Progress(n_missing; desc = "Evaluating grid points (Bumper): ", showspeed = true)
    end

    checkpoint_num = 0

    # Sort scenarios once to group related computations
    # Create sortable keys by converting structs to strings
    function make_sort_key(s)
        return (
            string(s.ensemble_specification),
            string(s.null_specification),
            string(s.noise_specification),
            s.percent_tested,
            string(s.test_specification),
            string(s.ews_metric_specification),
            string(s.ews_enddate_type),
            s.ews_metric,
            string(s.ews_threshold_window),
            string(s.ews_threshold_burnin),
            s.consecutive_thresholds,
            s.threshold_quantile,
        )
    end

    indexed_scenarios = [(i, make_sort_key(s)) for (i, s) in enumerate(missing_scenarios)]
    sort!(indexed_scenarios, by = x -> x[2])
    sorted_indices = [i for (i, _) in indexed_scenarios]

    @unpack emergent_incidence_arr, null_incidence_arr = data_arrs

    # Track current computation state to avoid redundant work
    current_noise_key = nothing
    current_testing_key = nothing
    current_ews_key = nothing
    current_scenario_key = nothing

    current_noisearr = nothing
    current_testarr = nothing
    current_null_testarr = nothing
    current_ews_metrics = nothing
    current_null_ews_metrics = nothing

    # Process sorted scenarios
    i = 1
    while i <= n_missing
        idx = sorted_indices[i]
        scenario = missing_scenarios[idx]

        # Check if we need new noise array
        noise_key = (
            scenario.ensemble_specification, scenario.null_specification, scenario.noise_specification,
        )
        if noise_key != current_noise_key
            current_noisearr = create_noise_arr(
                scenario.noise_specification,
                emergent_incidence_arr,
                scenario.ensemble_specification;
                seed = 1234,
            )[1]
            current_noise_key = noise_key
            current_testing_key = nothing
        end

        # Check if we need new testing arrays
        testing_key = (noise_key..., scenario.percent_tested, scenario.test_specification)
        if testing_key != current_testing_key
            current_testarr = create_testing_arrs(
                emergent_incidence_arr,
                current_noisearr,
                scenario.percent_tested,
                scenario.test_specification,
            )
            current_null_testarr = create_testing_arrs(
                null_incidence_arr,
                current_noisearr,
                scenario.percent_tested,
                scenario.test_specification,
            )
            current_testing_key = testing_key
            current_ews_key = nothing
        end

        # Check if we need new EWS metrics
        ews_key = (testing_key..., scenario.ews_metric_specification, scenario.ews_enddate_type)
        if ews_key != current_ews_key
            current_ews_metrics, current_null_ews_metrics = generate_ensemble_ews_metrics(
                data_arrs,
                current_testarr,
                current_null_testarr,
                scenario.ews_metric_specification,
                scenario.ews_enddate_type,
            )

            isempty(current_ews_metrics) &&
                error("No valid EWS metrics could be computed for any emergent simulation")
            isempty(current_null_ews_metrics) &&
                error("No valid EWS metrics could be computed for any null simulation")

            current_ews_key = ews_key
            current_scenario_key = nothing
        end

        # Find all scenarios with same expensive computations
        scenario_key = (
            ews_key..., scenario.ews_metric, scenario.ews_threshold_window,
            scenario.ews_threshold_burnin, scenario.consecutive_thresholds,
        )

        if scenario_key != current_scenario_key
            # Count matching scenarios first
            n_matching = 0
            j = i
            while j <= n_missing
                s = missing_scenarios[sorted_indices[j]]
                if (
                        s.ensemble_specification, s.null_specification, s.noise_specification,
                        s.percent_tested, s.test_specification, s.ews_metric_specification,
                        s.ews_enddate_type, s.ews_metric, s.ews_threshold_window,
                        s.ews_threshold_burnin, s.consecutive_thresholds,
                    ) == scenario_key
                    n_matching += 1
                    j += 1
                else
                    break
                end
            end

            # Use Bumper only for isbits types (Int indices)
            # OptimizationResult is not isbits, so use regular allocation
            batch_results = @no_escape begin
                batch_indices = @alloc(Int, n_matching)

                for k in 1:n_matching
                    batch_indices[k] = sorted_indices[i + k - 1]
                end

                # Use regular allocation for OptimizationResult since it's not isbits
                batch_results_temp = Vector{OptimizationResult}(undef, n_matching)

                for k in 1:n_matching
                    local s = missing_scenarios[batch_indices[k]]
                    local grid_scenario = GridSearchScenario(
                        s.ensemble_specification,
                        s.null_specification,
                        s.noise_specification,
                        s.test_specification,
                        s.percent_tested,
                        s.ews_metric_specification,
                        s.ews_enddate_type,
                        s.ews_threshold_window,
                        s.ews_threshold_burnin,
                        s.ews_metric,
                        s.threshold_quantile,
                        s.consecutive_thresholds,
                    )

                    batch_results_temp[k] = gridsearch_optimization_function(
                        grid_scenario, current_ews_metrics, current_null_ews_metrics
                    )
                end

                batch_results_temp
            end

            # Copy results to main array
            for k in 1:n_matching
                all_results[result_idx + k] = batch_results[k]
            end
            result_idx += n_matching

            # Save checkpoint after each batch
            if save_results && save_checkpoints && !isempty(checkpoint_dir)
                save_checkpoint_structvector(
                    StructVector(@view all_results[1:result_idx]), checkpoint_dir, checkpoint_num
                )
                verbose && @info "Saved checkpoint after batch $checkpoint_num"
                checkpoint_num += 1
            end

            if verbose
                for _ in 1:n_matching
                    next!(prog)
                end
            end

            i += n_matching
            current_scenario_key = scenario_key
        else
            i += 1
        end
    end

    return StructVector(@view all_results[1:result_idx])
end

"""
    evaluate_single_gridsearch_scenario(scenario, data_arrs)

Evaluate a single grid search scenario without optimization.
Returns an OptimizationResult struct with the evaluated metrics.
"""
function evaluate_single_gridsearch_scenario(
        scenario::GridSearchScenario,
        data_arrs::NamedTuple
    )
    # Create OptimizationScenario from GridSearchScenario (without grid params)
    opt_scenario = OptimizationScenario(
        scenario.ensemble_specification,
        scenario.null_specification,
        scenario.noise_specification,
        scenario.test_specification,
        scenario.percent_tested,
        scenario.ews_metric_specification,
        scenario.ews_enddate_type,
        scenario.ews_threshold_window,
        scenario.ews_threshold_burnin,
        scenario.ews_metric,
    )

    # Pre-compute expensive simulation data once per scenario
    # Reuse function from multistart optimization
    cached_data = create_cached_simulation_data(opt_scenario, data_arrs)

    return gridsearch_optimization_function(scenario, cached_data)
end

function gridsearch_optimization_function(
        scenario::GridSearchScenario,
        ews_metrics,
        null_ews_metrics,
    )
    ews_classification_results = calculate_ews_classification_results(
        scenario,
        ews_metrics,
        null_ews_metrics,
    )

    sensitivity = calculate_sensitivity(ews_classification_results)
    specificity = calculate_specificity(ews_classification_results)
    balanced_accuracy = calculate_balanced_accuracy(sensitivity, specificity)

    return OptimizationResult(
        # From scenario
        scenario.ensemble_specification,
        scenario.null_specification,
        scenario.noise_specification,
        scenario.test_specification,
        scenario.percent_tested,
        scenario.ews_metric_specification,
        scenario.ews_enddate_type,
        scenario.ews_threshold_window,
        scenario.ews_threshold_burnin,
        scenario.ews_metric,
        scenario.threshold_quantile,
        scenario.consecutive_thresholds,
        # Calculated values
        balanced_accuracy,
        sensitivity,
        specificity,
    )
end

# Grouping helper functions for hierarchical computation optimization

function group_by_noise_key(scenarios)
    groups = Dict{NoiseComputationKey, Vector{GridSearchScenario}}()
    for scenario in scenarios
        key = NoiseComputationKey(
            scenario.ensemble_specification,
            scenario.null_specification,
            scenario.noise_specification
        )
        push!(get!(groups, key, GridSearchScenario[]), scenario)
    end
    return groups
end

function group_by_testing_key(scenarios)
    groups = Dict{TestingComputationKey, Vector{GridSearchScenario}}()
    for scenario in scenarios
        noise_key = NoiseComputationKey(
            scenario.ensemble_specification,
            scenario.null_specification,
            scenario.noise_specification
        )
        key = TestingComputationKey(
            noise_key,
            scenario.percent_tested,
            scenario.test_specification
        )
        push!(get!(groups, key, GridSearchScenario[]), scenario)
    end
    return groups
end

function group_by_ews_key(scenarios)
    groups = Dict{EWSComputationKey, Vector{GridSearchScenario}}()
    for scenario in scenarios
        noise_key = NoiseComputationKey(
            scenario.ensemble_specification,
            scenario.null_specification,
            scenario.noise_specification
        )
        testing_key = TestingComputationKey(
            noise_key,
            scenario.percent_tested,
            scenario.test_specification
        )
        key = EWSComputationKey(
            testing_key,
            scenario.ews_metric_specification,
            scenario.ews_enddate_type
        )
        push!(get!(groups, key, GridSearchScenario[]), scenario)
    end
    return groups
end

function group_by_scenario_key(scenarios)
    groups = Dict{ScenarioComputationKey, Vector{GridSearchScenario}}()
    for scenario in scenarios
        noise_key = NoiseComputationKey(
            scenario.ensemble_specification,
            scenario.null_specification,
            scenario.noise_specification
        )
        testing_key = TestingComputationKey(
            noise_key,
            scenario.percent_tested,
            scenario.test_specification
        )
        ews_key = EWSComputationKey(
            testing_key,
            scenario.ews_metric_specification,
            scenario.ews_enddate_type
        )
        key = ScenarioComputationKey(
            ews_key,
            scenario.ews_metric,
            scenario.ews_threshold_window,
            scenario.ews_threshold_burnin,
            scenario.consecutive_thresholds
        )
        push!(get!(groups, key, GridSearchScenario[]), scenario)
    end
    return groups
end
