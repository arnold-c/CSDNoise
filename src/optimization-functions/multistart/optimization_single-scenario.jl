using MultistartOptimization: MultistartOptimization

export optimize_single_scenario

"""
    optimize_single_scenario(scenario, data_arrs, bounds, config)

Optimize EWS parameters for a single scenario using multistart optimization.
Returns an OptimizationResult struct containing both scenario parameters and optimized values.
"""
function optimize_single_scenario(
        scenario::OptimizationScenario,
        data_arrs::T1,
        bounds::T2,
        config::T3
    ) where {T1 <: NamedTuple, T2 <: NamedTuple, T3 <: NamedTuple}
    @unpack ews_metric = scenario

    # Pre-compute expensive simulation data once per scenario
    cached_data = create_cached_simulation_data(scenario, data_arrs)

    # Create tracker instance for this scenario
    tracker = OptimizationTracker()

    # Create objective function closure that updates tracker
    objective = params -> ews_objective_function_with_tracking(
        params,
        scenario,
        cached_data,
        tracker
    )

    # Setup multistart problem
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        bounds.lowers,
        bounds.uppers
    )

    # Configure local optimization method
    local_method = MultistartOptimization.NLopt_local_method(
        config.local_algorithm;
        xtol_rel = config.xtol_rel,
        xtol_abs = config.xtol_abs,
        maxeval = config.maxeval,
    )

    # Configure multistart method (TikTak uses Sobol sequences)
    multistart_method = MultistartOptimization.TikTak(config.n_sobol_points)

    # Run optimization
    MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    # Extract optimal parameters from tracker (which has the best metrics)
    optimal_params = map_continuous_to_ews_parameters(tracker.best_params)

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
        # From optimized values
        optimal_params.threshold_quantile,
        optimal_params.consecutive_thresholds,
        tracker.best_accuracy,
        tracker.best_sensitivity,
        tracker.best_specificity
    )
end
