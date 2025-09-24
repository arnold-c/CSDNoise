using Try: Try

export find_missing_scenarios,
    scenario_in_results,
    scenario_equals_optimization_result,
    create_cached_simulation_data

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
        scenario.null_specification == result.null_specification &&
        scenario.noise_specification == result.noise_specification &&
        scenario.test_specification == result.test_specification &&
        scenario.percent_tested == result.percent_tested &&
        scenario.ews_metric_specification == result.ews_metric_specification &&
        scenario.ews_enddate_type == result.ews_enddate_type &&
        scenario.ews_threshold_window == result.ews_threshold_window &&
        scenario.ews_threshold_burnin == result.ews_threshold_burnin &&
        scenario.ews_metric == result.ews_metric &&
        scenario.threshold_quantile == result.threshold_quantile &&
        scenario.consecutive_thresholds == result.consecutive_thresholds
end

function scenario_equals_optimization_result(
        scenario::OptimizationScenario,
        result::OptimizationResult
    )::Bool
    return scenario.ensemble_specification == result.ensemble_specification &&
        scenario.null_specification == result.null_specification &&
        scenario.noise_specification == result.noise_specification &&
        scenario.test_specification == result.test_specification &&
        scenario.percent_tested == result.percent_tested &&
        scenario.ews_metric_specification == result.ews_metric_specification &&
        scenario.ews_enddate_type == result.ews_enddate_type &&
        scenario.ews_threshold_window == result.ews_threshold_window &&
        scenario.ews_threshold_burnin == result.ews_threshold_burnin &&
        scenario.ews_metric == result.ews_metric
end

"""
    create_cached_simulation_data(scenario, data_arrs)

Pre-compute expensive simulation data once per scenario to avoid repeated computation
during parameter optimization. This includes noise arrays, test arrays, and EWS metrics.
"""
function create_cached_simulation_data(
        scenario::OptimizationScenario,
        data_arrs::NamedTuple
    )
    @unpack noise_specification, test_specification, percent_tested,
        ews_enddate_type, ews_metric_specification = scenario

    @unpack ensemble_single_incarr, null_single_incarr,
        ensemble_specification, ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs = data_arrs

    # Create noise array once per scenario (expensive operation)
    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr,
        ensemble_specification;
        seed = 1234,
    )[1]

    # Create test arrays once per scenario (expensive operation)
    testarr = create_testing_arrs(
        ensemble_single_incarr,
        noisearr,
        percent_tested,
        test_specification,
    )

    null_testarr = create_testing_arrs(
        null_single_incarr,
        noisearr,
        percent_tested,
        test_specification,
    )

    # Must be the same size as need to pair each null series to a
    # emergent series when calculating the end date for the EWS
    # metric calculation
    @assert size(testarr, 3) == size(null_testarr, 3)

    # Select appropriate thresholds
    thresholds = get_enddate_thresholds(ews_enddate_type, data_arrs)

    # Pre-compute EWS metrics for all simulations
    ensemble_nsims = size(ensemble_single_incarr, 3)
    ews_metrics = Vector{EWSMetrics}()
    null_ews_metrics = Vector{EWSMetrics}()

    for sim in 1:ensemble_nsims
        ews_enddate = calculate_ews_enddate(thresholds[sim], ews_enddate_type)
        if Try.isok(ews_enddate)
            ews_enddate_val = Try.unwrap(ews_enddate)

            ews_metric = calculate_ews_metrics_for_simulation(
                ews_metric_specification,
                @view(testarr[:, :, sim]),
                thresholds[sim],
                ews_enddate_val
            )
            push!(ews_metrics, ews_metric)

            null_ews_metric = calculate_ews_metrics_for_simulation(
                ews_metric_specification,
                @view(null_testarr[:, :, sim]),
                thresholds[sim],
                ews_enddate_val
            )
            push!(null_ews_metrics, null_ews_metric)
        end
    end

    # Ensure we have at least some valid simulations
    if isempty(ews_metrics)
        error("No valid EWS metrics could be computed for any simulation")
    end

    return CachedSimulationData(
        testarr,
        null_testarr,
        thresholds,
        ews_metrics,
        null_ews_metrics
    )
end

function generate_ensemble_ews_metrics(
        data_arrs::NamedTuple{T},
        testarr,
        null_testarr,
        ews_metric_specification,
        ews_enddate_type,
    ) where {T}
    thresholds = get_enddate_thresholds(ews_enddate_type, data_arrs)

    ensemble_nsims = size(testarr, 3)
    ews_metrics = Vector{EWSMetrics}()
    null_ews_metrics = Vector{EWSMetrics}()

    for sim in 1:ensemble_nsims
        ews_enddate = calculate_ews_enddate(thresholds[sim], ews_enddate_type)
        if Try.isok(ews_enddate)
            ews_enddate_val = Try.unwrap(ews_enddate)

            ews_metric = calculate_ews_metrics_for_simulation(
                ews_metric_specification,
                @view(testarr[:, :, sim]),
                thresholds[sim],
                ews_enddate_val
            )
            push!(ews_metrics, ews_metric)

            null_ews_metric = calculate_ews_metrics_for_simulation(
                ews_metric_specification,
                @view(null_testarr[:, :, sim]),
                thresholds[sim],
                ews_enddate_val
            )
            push!(null_ews_metrics, null_ews_metric)
        end
    end

    return ews_metrics, null_ews_metrics
end
