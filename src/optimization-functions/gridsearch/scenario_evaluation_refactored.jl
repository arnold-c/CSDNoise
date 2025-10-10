export evaluate_gridsearch_scenarios_refactored

function evaluate_gridsearch_scenarios_refactored(
        missing_scenarios::StructVector{GridSearchScenario};
        use_threads = false,
        save_results = true,
        save_checkpoints = false,
        checkpoint_dir::String = "",
        verbose::Bool = true,
        seed = 1234
    )

    n_missing = length(missing_scenarios)
    n_missing == 0 && return StructVector(OptimizationResult[])

    all_results = OptimizationResult[]

    if verbose
        prog = ProgressMeter.Progress(n_missing; desc = "Evaluating grid points: ", showspeed = true)
    end

    checkpoint_num = 0

    ensemble_groups = group_structvector(missing_scenarios, :ensemble_specification)

    for (ensemble_key, ensemble_scenarios) in ensemble_groups
        ensemble_simulation = generate_single_ensemble(
            ensemble_key.ensemble_specification;
            seed = seed
        )

        Reff_thresholds_vec = calculate_all_Reff_thresholds(ensemble_simulation)

        noise_trim_groups = group_structvector(
            ensemble_scenarios,
            :noise_level, :noise_type_description, :ews_enddate_type
        )

        for (noise_trim_key, noise_trim_scenarios) in noise_trim_groups

            enddates_vec = calculate_all_ews_enddates(
                Reff_thresholds_vec,
                noise_trim_key.ews_enddate_type
            )

            trimmed_ensemble = trim_ensemble_simulations(
                ensemble_simulation,
                enddates_vec
            )


            noise_specification = create_noise_specification(
                noise_trim_key.noise_level,
                noise_trim_key.noise_type_description,
                trimmed_ensemble
            )

            noisearr = create_noise_arr(
                noise_specification,
                trimmed_ensemble.emergent_incidence_arr,
                ensemble_key.ensemble_specification;
                seed = seed
            )[1]

            test_groups = group_structvector(
                noise_trim_scenarios,
                :percent_tested, :test_specification
            )

            for (test_key, test_scenarios) in test_groups
                testarr = create_testing_arrs(
                    trimmed_ensemble.emergent_incidence_arr,
                    noisearr,
                    test_key.percent_tested,
                    test_key.test_specification
                )

                null_testarr = create_testing_arrs(
                    trimmed_ensemble.null_incidence_arr,
                    noisearr,
                    test_key.percent_tested,
                    test_key.test_specification
                )

                ews_spec_groups = group_structvector(
                    test_scenarios,
                    :ews_metric_specification
                )

                for (ews_spec_key, ews_spec_scenarios) in ews_spec_groups
                    ews_metrics, null_ews_metrics = generate_ensemble_ews_metrics(
                        trimmed_ensemble.data_arrs,
                        testarr,
                        null_testarr,
                        ews_spec_key.ews_metric_specification,
                        noise_trim_key.ews_enddate_type
                    )

                    isempty(ews_metrics) && error("No valid EWS metrics for emergent")
                    isempty(null_ews_metrics) && error("No valid EWS metrics for null")

                    opt_groups = group_structvector(
                        ews_spec_scenarios,
                        :ews_metric, :ews_threshold_window, :ews_threshold_burnin,
                        :consecutive_thresholds
                    )

                    for (opt_key, opt_scenarios) in opt_groups
                        unique_quantiles = unique(opt_scenarios.threshold_quantile)

                        optimization_results = map(unique_quantiles) do threshold_quantile
                            grid_scenario = first(
                                filter(
                                    s -> s.threshold_quantile == threshold_quantile,
                                    opt_scenarios
                                )
                            )

                            gridsearch_optimization_function(
                                grid_scenario,
                                ews_metrics,
                                null_ews_metrics
                            )
                        end

                        BangBang.append!!(all_results, optimization_results)

                        if save_results && save_checkpoints && !isempty(checkpoint_dir)
                            save_checkpoint_structvector(
                                StructVector(all_results),
                                checkpoint_dir,
                                checkpoint_num
                            )
                            verbose && @info "Saved checkpoint $checkpoint_num"
                            checkpoint_num += 1
                        end

                        verbose && ProgressMeter.update!(prog, length(all_results))
                    end
                end
            end
        end
    end

    return StructVector(all_results)
end
