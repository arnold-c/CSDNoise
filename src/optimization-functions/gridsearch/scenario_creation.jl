export create_gridsearch_scenarios_structvector

"""
    create_gridsearch_scenarios_structvector(
    	specification_vecs;
		executor=ThreadedEx()
    )

Create all grid search scenarios including parameter combinations.
Returns StructVector{GridSearchScenario} with all combinations.
"""
function create_gridsearch_scenarios_structvector(specification_vecs::GridSearchSpecificationVecs)
    @unpack ensemble_specification_vec,
        noise_level_vec,
        noise_type_description_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_metric_vec,
        ews_threshold_quantile_vec,
        ews_consecutive_thresholds_vec = specification_vecs

    combinations = Iterators.product(
        ensemble_specification_vec,
        noise_level_vec,
        noise_type_description_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_metric_vec,
        ews_threshold_quantile_vec,
        ews_consecutive_thresholds_vec
    )
    n_combinations = length(combinations)

    scenarios_vec = Vector{GridSearchScenario}(undef, n_combinations)

    for (
            i, (
                ensemble_spec,
                noise_level,
                noise_type_description,
                test_spec,
                percent_tested,
                ews_metric_spec,
                ews_enddate_type,
                ews_window,
                ews_metric,
                threshold_quantile,
                consecutive_thresholds,
            ),
        ) in enumerate(combinations)

        scenarios_vec[i] = GridSearchScenario(
            ensemble_specification = ensemble_spec,
            noise_level = noise_level,
            noise_type_description = noise_type_description,
            test_specification = test_spec,
            percent_tested = percent_tested,
            ews_metric_specification = ews_metric_spec,
            ews_enddate_type = ews_enddate_type,
            ews_threshold_window = ews_window,
            ews_metric = ews_metric,
            threshold_quantile = threshold_quantile,
            consecutive_thresholds = consecutive_thresholds,
        )
    end

    return StructVector(scenarios_vec)
end
