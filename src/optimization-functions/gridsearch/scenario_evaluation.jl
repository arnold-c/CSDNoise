export evaluate_gridsearch_scenarios

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
        missing_scenarios::StructVector{GridSearchScenario};
        use_threads = false,
        save_results = true,
        save_checkpoints = false,
        checkpoint_dir::String = "",
        verbose::Bool = true
    )

    n_missing = length(missing_scenarios)

    if n_missing == 0
        return StructVector(OptimizationResult[])
    end

    all_results = OptimizationResult[]

    # Setup progress tracking
    if verbose
        prog = ProgressMeter.Progress(n_missing; desc = "Evaluating grid points: ", showspeed = true)
    end

    checkpoint_num = 0

    # Creating the incidence, noise, testing arrays and EWS metrics is expensive
    # Loop through all unique combinations in a nested loop to allow re-use within
    # the optimization and grid search steps
    unique_ensemble_specs = unique(missing_scenarios.ensemble_specification)
    for ensemble_specification in unique_ensemble_specs
        ensemble_specific_scenarios = filter(
            s -> s.ensemble_specification == ensemble_specification,
            missing_scenarios
        )

        # create the target disease simulations

        # create a loop that zips the noise level and enddate_type to first trim the
        # target disease simulation and then creates the noise simulations

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
                                            for (i, threshold_quantile) in pairs(unique_threshold_quantiles)
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
                                                ProgressMeter.next!(prog)
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
