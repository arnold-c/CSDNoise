export find_missing_scenarios

"""
    find_missing_scenarios(all_scenarios, completed_results)

Find optimization (grid search or optimization) scenarios that haven't been computed yet.
"""
function find_missing_scenarios(
        all_scenarios::StructVector{T},
        completed_results::StructVector{OptimizationResult}
    ) where {T <: Union{GridSearchScenario, OptimizationScenario}}
    # )::StructVector{T} where {T <: Union{GridSearchScenario, OptimizationScenario}}
    if isempty(completed_results)
        return all_scenarios
    end

    # Build a dictionary of completed scenarios for O(1) lookup
    completed_dict = build_scenario_dict(completed_results, eltype(all_scenarios))

    # Filter to find missing scenarios using dictionary lookup
    missing_mask = map(all_scenarios) do scenario
        !haskey(completed_dict, scenario_key(scenario))
    end

    return all_scenarios[missing_mask]
end

"""
    build_scenario_dict(results, scenario_type)

Build a dictionary mapping scenario keys to true for all completed results.
"""
function build_scenario_dict(
        results::StructVector{OptimizationResult},
        ::Type{T}
    )::Dict{Tuple, Bool} where {T <: Union{GridSearchScenario, OptimizationScenario}}
    dict = Dict{Tuple, Bool}()

    for result in results
        dict[result_key(result, T)] = true
    end
    return dict
end

"""
    scenario_key(scenario)

Generate a tuple key from a scenario for dictionary lookup.
"""
function scenario_key(scenario::GridSearchScenario)
    return (
        scenario.ensemble_specification,
        scenario.noise_level,
        scenario.noise_type_description,
        scenario.test_specification,
        scenario.percent_tested,
        scenario.ews_metric_specification,
        scenario.ews_enddate_type,
        scenario.ews_threshold_window,
        scenario.ews_metric,
        scenario.threshold_quantile,
        scenario.consecutive_thresholds,
    )
end

function scenario_key(scenario::OptimizationScenario)
    return (
        scenario.ensemble_specification,
        scenario.noise_level,
        scenario.noise_type_description,
        scenario.test_specification,
        scenario.percent_tested,
        scenario.ews_metric_specification,
        scenario.ews_enddate_type,
        scenario.ews_threshold_window,
        scenario.ews_metric,
    )
end

"""
    result_key(result, scenario_type)

Generate a tuple key from an OptimizationResult matching the scenario type.
"""
function result_key(
        result::OptimizationResult,
        ::Type{GridSearchScenario}
    )
    return (
        result.ensemble_specification,
        result.noise_level,
        result.noise_type_description,
        result.test_specification,
        result.percent_tested,
        result.ews_metric_specification,
        result.ews_enddate_type,
        result.ews_threshold_window,
        result.ews_metric,
        result.threshold_quantile,
        result.consecutive_thresholds,
    )
end

function result_key(
        result::OptimizationResult,
        ::Type{OptimizationScenario}
    )
    return (
        result.ensemble_specification,
        result.noise_level,
        result.noise_type_description,
        result.test_specification,
        result.percent_tested,
        result.ews_metric_specification,
        result.ews_enddate_type,
        result.ews_threshold_window,
        result.ews_metric,
    )
end
