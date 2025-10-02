export create_gridsearch_scenarios_structvector

"""
    create_gridsearch_scenarios_structvector(
    	specification_vecs;
		executor=ThreadedEx()
    )

Create all grid search scenarios including parameter combinations.
Returns StructVector{GridSearchScenario} with all combinations.
"""
function create_gridsearch_scenarios_structvector(specification_vecs)
    @unpack ensemble_noise_spec_pairs_vec,
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
        ensemble_noise_spec_pairs_vec,
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
