export ews_objective_function_with_tracking,
    map_continuous_to_ews_parameters

"""
    ews_objective_function_with_tracking(params_vec, scenario, cached_data, tracker)

Objective function for EWS parameter optimization that tracks the best solution.
Returns 1 - accuracy for minimization and updates tracker with best metrics.
Uses pre-computed cached simulation data for efficiency.
"""
function ews_objective_function_with_tracking(
        params_vec::Vector{Float64},
        scenario::OptimizationScenario,
        cached_data::CachedSimulationData,
        tracker::OptimizationTracker
    )
    @unpack ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window, ensemble_specification, ews_metric = scenario

    @unpack burnin = ensemble_specification.time_parameters

    @unpack testarr,
        null_testarr,
        thresholds,
        ews_metrics,
        null_ews_metrics = cached_data

    # Map continuous parameters to EWS parameters
    # Note: We use continuous optimization with rounding for integer parameters.
    # This is appropriate since consecutive_thresholds has a small integer range (1-10)
    # and mixed-integer optimization would be more complex without significant benefit.
    ews_params = map_continuous_to_ews_parameters(params_vec)

    # Pre-compute values outside the loop for type stability
    ews_metric_symbol = Symbol(ews_metric)
    threshold_quantile = ews_params.threshold_quantile

    # Calculate accuracy
    n_emergent_sims = length(ews_metrics)
    n_null_sims = length(null_ews_metrics)
    true_positives = 0
    true_negatives = 0

    for sim in eachindex(ews_metrics)
        # Use pre-computed EWS metrics
        ews_vals = ews_metrics[sim]
        null_ews_vals = null_ews_metrics[sim]

        # Check threshold exceedances
        exceeds_threshold = exceeds_ews_threshold(
            ews_vals,
            ews_metric_symbol,
            ews_threshold_window,
            threshold_quantile,
            burnin,
        )

        detection_index = calculate_ews_trigger_index(
            exceeds_threshold,
            ews_params.consecutive_thresholds,
        )

        null_exceeds_threshold = exceeds_ews_threshold(
            null_ews_vals,
            ews_metric_symbol,
            ews_threshold_window,
            threshold_quantile,
            burnin,
        )

        null_detection_index = calculate_ews_trigger_index(
            null_exceeds_threshold,
            ews_params.consecutive_thresholds,
        )

        # Update counts
        if Try.isok(detection_index)
            true_positives += 1
        end
        if Try.iserr(null_detection_index)
            true_negatives += 1
        end
    end

    # Calculate metrics
    sensitivity = calculate_sensitivity(true_positives, n_emergent_sims)
    specificity = calculate_specificity(true_negatives, n_null_sims)
    accuracy = calculate_balanced_accuracy(sensitivity, specificity)
    loss = 1.0 - accuracy

    # Update tracker if this is the best solution so far
    if loss < tracker.best_loss
        tracker.best_loss = loss
        tracker.best_accuracy = accuracy
        tracker.best_sensitivity = sensitivity
        tracker.best_specificity = specificity
        tracker.best_params = copy(params_vec)
    end

    return loss  # Return scalar loss for optimizer
end


"""
    map_continuous_to_ews_parameters(params_vec)

Map continuous optimization parameters to discrete EWS parameters.
"""
function map_continuous_to_ews_parameters(params_vec::Vector{Float64})
    return (
        threshold_quantile = params_vec[1],  # Already in correct range
        consecutive_thresholds = round(Int, params_vec[2]),
    )
end
