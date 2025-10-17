export gridsearch_optimization

function gridsearch_optimization(
        optimization_scenario::GridSearchScenario,
        ensemble_ews_metrics::EnsembleEWSMetrics,
    )
    ews_classification_result = calculate_ews_classification_results(
        optimization_scenario,
        ensemble_ews_metrics
    )

    sensitivity = calculate_sensitivity(ews_classification_result)
    specificity = calculate_specificity(ews_classification_result)
    balanced_accuracy = calculate_balanced_accuracy(sensitivity, specificity)

    optimization_result = OptimizationResult(
        ensemble_specification = optimization_scenario.ensemble_specification,
        noise_level = optimization_scenario.noise_level,
        noise_type_description = optimization_scenario.noise_type_description,
        test_specification = optimization_scenario.test_specification,
        percent_tested = optimization_scenario.percent_tested,
        ews_metric_specification = optimization_scenario.ews_metric_specification,
        ews_enddate_type = optimization_scenario.ews_enddate_type,
        ews_threshold_window = optimization_scenario.ews_threshold_window,
        ews_metric = optimization_scenario.ews_metric,
        threshold_quantile = optimization_scenario.threshold_quantile,
        consecutive_thresholds = optimization_scenario.consecutive_thresholds,
        accuracy = balanced_accuracy,
        sensitivity = sensitivity,
        specificity = specificity
    )

    return optimization_result
end
