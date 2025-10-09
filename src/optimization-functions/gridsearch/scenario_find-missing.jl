export find_missing_scenarios

"""
    find_missing_scenarios(all_scenarios, completed_results)

Find optimization (grid search or optimization) scenarios that haven't been computed yet.
"""
function find_missing_scenarios(
        all_scenarios::StructVector{T},
        completed_results::StructVector{OptimizationResult}
    )::StructVector{T} where {T <: Union{GridSearchScenario, OptimizationScenario}}
    if isempty(completed_results)
        return all_scenarios
    end

    # Filter to find missing scenarios using the helper function
    missing_mask = map(all_scenarios) do scenario
        !scenario_in_results(scenario, completed_results)
    end

    return all_scenarios[missing_mask]
end

"""
    scenario_in_results(scenario, results)

Check if a GridSearchScenario or OptimizationScenario is present in a StructVector of OptimizationResult.
Returns true if any OptimizationResult matches the scenario.
"""
function scenario_in_results(
        scenario::T,
        results::StructVector{OptimizationResult}
    )::Bool where {T <: Union{GridSearchScenario, OptimizationScenario}}
    return any(result -> scenario_equals_optimization_result(scenario, result), results)
end

"""
    scenario_equals_optimization_result(scenario, result)

Check if a GridSearchScenario of OptimizationScenario instance equals an OptimizationResult instance
by comparing all shared fields.
"""
function scenario_equals_optimization_result(
        scenario::GridSearchScenario,
        result::OptimizationResult
    )::Bool
    return scenario.ensemble_specification == result.ensemble_specification &&
        scenario.noise_level == result.noise_level &&
        scenario.test_specification == result.test_specification &&
        scenario.percent_tested == result.percent_tested &&
        scenario.ews_metric_specification == result.ews_metric_specification &&
        scenario.ews_enddate_type == result.ews_enddate_type &&
        scenario.ews_threshold_window == result.ews_threshold_window &&
        scenario.ews_metric == result.ews_metric &&
        scenario.threshold_quantile == result.threshold_quantile &&
        scenario.consecutive_thresholds == result.consecutive_thresholds
end

function scenario_equals_optimization_result(
        scenario::OptimizationScenario,
        result::OptimizationResult
    )::Bool
    return scenario.ensemble_specification == result.ensemble_specification &&
        scenario.noise_level == result.noise_level &&
        scenario.test_specification == result.test_specification &&
        scenario.percent_tested == result.percent_tested &&
        scenario.ews_metric_specification == result.ews_metric_specification &&
        scenario.ews_enddate_type == result.ews_enddate_type &&
        scenario.ews_threshold_window == result.ews_threshold_window &&
        scenario.ews_metric == result.ews_metric
end
