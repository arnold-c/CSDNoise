export evaluate_gridsearch_scenarios_refactored

function evaluate_gridsearch_scenarios_refactored(
        missing_scenarios::StructVector{GridSearchScenario};
        scheduler = :dynamic, #:serial, :greedy, :static, :dynamic
        save_results = true,
        save_checkpoints = false,
        save_checkpoint_num = 5,
        checkpoint_dir::String = "",
        checkpoint_output_filename_base = joinpath(
            filedir,
            string(Dates.now()) * "_" * "checkpoint_batch_",
        ),
        verbose::Bool = true,
        verbose_noise_optimization = false,
        seed = 1234,
        dynamic_noise_optimization_parameters::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters()

    )
    @assert scheduler in [:dynamic, :static, :greedy, :serial]

    n_missing = length(missing_scenarios)
    n_missing == 0 && return StructVector(OptimizationResult[])

    # Intialize with empty 0-length vector as could be unspecified number of quantiles
    # for a given scenario
    all_results = OptimizationResult[]

    if verbose
        prog = ProgressMeter.Progress(n_missing; desc = "Evaluating grid points: ", showspeed = true)
    end

    checkpoint_num = 1

    ensemble_groups = group_structvector(missing_scenarios, :ensemble_specification)

    for (ensemble_key, ensemble_scenarios) in ensemble_groups
        verbose && println("Ensemble specification: $(ensemble_key.ensemble_specification.label)")

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
            verbose && println("\tEWS enddate type: $(noise_trim_key.ews_enddate_type)\n\tNoise level: $(noise_trim_key.noise_level)\n\tNoise type: $(noise_trim_key.noise_type_description)")

            enddates_vec_result = calculate_all_ews_enddates(
                Reff_thresholds_vec,
                noise_trim_key.ews_enddate_type
            )
            if Try.iserr(enddates_vec_result)
                error_message =
                """
                Failed to calculate EWS enddates for ensemble with parameters:
                  Ensemble specification: $(ensemble_key.ensemble_specification.label)
                  EWS enddate type: $(noise_trim_key.ews_enddate_type)
                  Noise level: $(noise_trim_key.noise_level)
                  Noise type: $(noise_trim_key.noise_type_description)
                  Reff thresholds type: $(typeof(Reff_thresholds_vec))

                Original error: $(Try.unwrap_err(enddates_vec_result))

                This typically occurs when Reff doesn't cross 1.0 within the simulation time.
                Consider:
                  - Lowering the vaccination rate post-burnin
                  - Increasing the simulation length
                  - Adjusting the dynamics parameters
                """
                error(error_message)
            end
            enddates_vec = Try.unwrap(enddates_vec_result)

            trimmed_ensemble = trim_ensemble_simulations(
                ensemble_simulation,
                enddates_vec
            )

            noise_vecs = if noise_trim_key.noise_type_description == :static
                create_noise_vecs(
                    PoissonNoise(noise_trim_key.noise_level),
                    ensemble_key.ensemble_specification,
                    enddates_vec,
                    trimmed_ensemble.emergent_seir_run,
                    seed = seed
                )
            else
                # Call the original optimization function with the computed mean
                optim_res = optimize_dynamic_noise_params_wrapper(
                    ensemble_key.ensemble_specification,
                    trimmed_ensemble.emergent_seir_run,
                    enddates_vec,
                    noise_trim_key.noise_level,
                    dynamic_noise_optimization_parameters;
                    verbose = verbose_noise_optimization,
                    seed = seed
                )

                recreate_noise_vecs(
                    ensemble_key.ensemble_specification,
                    enddates_vec,
                    optim_res.location[1],
                )
            end

            test_groups = group_structvector(
                noise_trim_scenarios,
                :percent_tested, :test_specification
            )

            for (test_key, test_scenarios) in test_groups
                verbose && println("\t\tTest Specification: $(test_key.test_specification)\n\t\tTest Percentage: $(test_key.percent_tested)")

                ensemble_test_positives = create_test_positive_vecs(
                    trimmed_ensemble,
                    noise_vecs,
                    test_key.percent_tested,
                    test_key.test_specification,
                )

                ews_spec_groups = group_structvector(
                    test_scenarios,
                    :ews_metric_specification
                )

                for (ews_spec_key, ews_spec_scenarios) in ews_spec_groups
                    verbose && println("\t\t\tEWS specification: $(ews_spec_key.ews_metric_specification)")
                    ensemble_ews_metrics = generate_ensemble_ews_metrics(
                        ews_spec_key.ews_metric_specification,
                        ensemble_test_positives,
                    )

                    opt_groups = group_structvector(
                        ews_spec_scenarios,
                        :ews_metric,
                        :ews_threshold_window,
                        :consecutive_thresholds,
                        :threshold_quantile
                    )

                    opt_scenarios_vec = collect(values(opt_groups))

                    results_batch = OhMyThreads.tmap(
                        opt_scenarios_vec; scheduler = scheduler
                    ) do opt_scenario_sv
                        @assert length(opt_scenario_sv) == 1
                        optimization_scenario = opt_scenario_sv[1]

                        optimization_result = gridsearch_optimization(
                            optimization_scenario,
                            ensemble_ews_metrics,
                        )

                        return optimization_result
                    end

                    BangBang.append!!(all_results, results_batch)

                    if save_results && save_checkpoints && checkpoint_num % save_checkpoint_num == 0 && !isempty(checkpoint_dir)
                        save_checkpoint_structvector(
                            StructVector(all_results),
                            checkpoint_dir,
                            checkpoint_output_filename_base,
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

    return StructVector(all_results)
end
