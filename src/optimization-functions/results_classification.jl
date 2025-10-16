using UnPack: @unpack
using Try: Try

export calculate_ews_classification_results,
    calculate_sensitivity,
    calculate_specificity,
    calculate_balanced_accuracy

"""
    calculate_ews_classification_results(gridsearch_scenario, ensemble_ews_metrics)

Calculate EWS classification performance metrics from ensemble EWS metrics.

This function evaluates the classification performance of early warning signals
by computing true positives, true negatives, false positives, and false negatives
based on threshold exceedances and trigger detection. It serves as a convenience
wrapper that extracts emergent and null metrics from an ensemble object.

# Arguments
- `gridsearch_scenario::GridSearchScenario`: Configuration specifying EWS parameters and thresholds
- `ensemble_ews_metrics::EnsembleEWSMetrics`: Container with both emergent and null EWS metrics

# Returns
- `EWSClassificationResults`: Object containing classification counts and totals

# Implementation
This method extracts the emergent and null EWS metrics from the ensemble object
and delegates to the main implementation method.

# Example
```julia
# Evaluate classification performance for a scenario
results = calculate_ews_classification_results(scenario, ensemble_metrics)
sensitivity = calculate_sensitivity(results)
specificity = calculate_specificity(results)
```

# See Also
- [`calculate_ews_classification_results`](@ref): Main implementation method
- [`EnsembleEWSMetrics`](@ref): Container type for ensemble metrics
- [`GridSearchScenario`](@ref): Configuration type for scenarios
"""
function calculate_ews_classification_results(
        gridsearch_scenario::GridSearchScenario,
        ensemble_ews_metrics::EnsembleEWSMetrics
    )
    return calculate_ews_classification_results(
        gridsearch_scenario,
        ensemble_ews_metrics.emergent_ews_metrics,
        ensemble_ews_metrics.null_ews_metrics
    )
end


"""
    calculate_ews_classification_results(gridsearch_scenario, emergent_ews_metrics, null_ews_metrics)

Calculate EWS classification performance metrics from separate emergent and null EWS metrics.

This function implements the core algorithm for evaluating early warning signal
classification performance. It processes paired emergent and null simulations,
applies threshold detection logic, and computes classification metrics including
true positives, true negatives, false positives, and false negatives.

# Arguments
- `gridsearch_scenario::GridSearchScenario`: Configuration containing:
  - `ews_metric`: The EWS metric to evaluate (e.g., "variance", "autocovariance")
  - `ews_threshold_window`: Window type for threshold calculation
  - `threshold_quantile`: Quantile level for threshold determination
  - `consecutive_thresholds`: Number of consecutive exceedances required for trigger
  - `ensemble_specification`: Contains burn-in period and other parameters
- `emergent_ews_metrics::StructVector{EWSMetrics}`: EWS metrics from emergent simulations
- `null_ews_metrics::StructVector{EWSMetrics}`: EWS metrics from null simulations

# Returns
- `EWSClassificationResults`: Object containing:
  - `true_positives`: Emergent simulations correctly detected
  - `true_negatives`: Null simulations correctly not detected
  - `false_positives`: Null simulations incorrectly detected
  - `false_negatives`: Emergent simulations incorrectly not detected
  - `n_emergent_sims`: Total number of emergent simulations
  - `n_null_sims`: Total number of null simulations

# Algorithm
For each paired simulation:
1. Extract the specified EWS metric from both emergent and null scenarios
2. Apply threshold detection using the configured window and quantile
3. Check for trigger using consecutive threshold exceedances
4. Update classification counts:
   - True positive: Emergent simulation triggers detection
   - True negative: Null simulation does not trigger detection
   - False positive: Null simulation triggers detection
   - False negative: Emergent simulation does not trigger detection

# Performance Considerations
- Processes simulations in pairs to ensure matched evaluation
- Uses pre-computed EWS metrics for efficiency
- Leverages Try.jl for robust error handling in trigger detection

# Example
```julia
# Configure scenario parameters
scenario = GridSearchScenario(
    ews_metric="variance",
    ews_threshold_window=EWSThresholdWindowType(ExpandingThresholdWindow()),
    threshold_quantile=0.95,
    consecutive_thresholds=2,
    ensemble_specification=ensemble_spec
)

# Calculate classification results
results = calculate_ews_classification_results(
    scenario,
    emergent_metrics,
    null_metrics
)

# Extract performance metrics
sensitivity = calculate_sensitivity(results)
specificity = calculate_specificity(results)
balanced_accuracy = calculate_balanced_accuracy(sensitivity, specificity)
```

# See Also
- [`exceeds_ews_threshold`](@ref): Determines threshold exceedances
- [`calculate_ews_trigger_index`](@ref): Detects consecutive exceedances
- [`EWSClassificationResults`](@ref): Result container type
- [`calculate_sensitivity`](@ref), [`calculate_specificity`](@ref): Performance metrics
"""
function calculate_ews_classification_results(
        gridsearch_scenario::GridSearchScenario,
        emergent_ews_metrics::StructVector{EWSMetrics},
        null_ews_metrics::StructVector{EWSMetrics},
    )
    @unpack ews_metric,
        ews_threshold_window,
        ensemble_specification,
        threshold_quantile,
        consecutive_thresholds = gridsearch_scenario
    @unpack burnin = ensemble_specification.time_parameters

    ews_metric_symbol = Symbol(ews_metric)

    n_emergent_sims = length(emergent_ews_metrics)
    n_null_sims = length(null_ews_metrics)

    true_positives = 0
    true_negatives = 0

    for sim in eachindex(emergent_ews_metrics)
        # Use pre-computed EWS metrics
        ews_vals = emergent_ews_metrics[sim]
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
            consecutive_thresholds,
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
            consecutive_thresholds,
        )

        # Update counts
        if Try.isok(detection_index)
            true_positives += 1
        end
        if Try.iserr(null_detection_index)
            true_negatives += 1
        end
    end

    return EWSClassificationResults(
        true_positives,
        true_negatives,
        n_null_sims - true_negatives,
        n_emergent_sims - true_positives,
        n_emergent_sims,
        n_null_sims
    )
end

"""
    calculate_balanced_accuracy(sensitivity, specificity)

Calculate the balanced accuracy from sensitivity and specificity values.

Balanced accuracy is the arithmetic mean of sensitivity (true positive rate) and
specificity (true negative rate), providing a single metric that accounts for
performance on both positive and negative cases. It is balanced because emergent
and null simulations are paired to ensure matching evaluation lengths, so there
are equal numbers of emergent and null simulations, and therefore equal weighting
in the accuracy calculation.

# Arguments
- `sensitivity`: The sensitivity (true positive rate) value, typically between 0 and 1
- `specificity`: The specificity (true negative rate) value, typically between 0 and 1

# Returns
- `Float64`: The balanced accuracy as (sensitivity + specificity) / 2

# Examples
```julia
balanced_accuracy = calculate_balanced_accuracy(0.8, 0.9)  # Returns 0.85
```
"""
function calculate_balanced_accuracy(sensitivity, specificity)
    return (sensitivity + specificity) / 2
end

function calculate_sensitivity(classification_results::EWSClassificationResults)
    return calculate_sensitivity(classification_results.true_positives, classification_results.n_emergent_sims)
end

"""
    calculate_sensitivity(true_positives, n_emergent_sims)

Calculate the sensitivity (true positive rate) from the number of true positives
and total number of emergent simulations.

Sensitivity measures the proportion of actual positive cases that are correctly
identified as positive. In the context of early warning systems, this represents
the fraction of emergent simulations where a critical transition was correctly detected.

# Arguments
- `true_positives`: Number of true positive detections
- `n_emergent_sims`: Total number of emergent simulations

# Returns
- `Float64`: The sensitivity value as true_positives / n_emergent_sims

# Examples
```julia
sensitivity = calculate_sensitivity(80, 100)  # Returns 0.8
```
"""
function calculate_sensitivity(true_positives, n_emergent_sims)
    return true_positives / n_emergent_sims

end

function calculate_specificity(classification_results::EWSClassificationResults)
    return calculate_specificity(classification_results.true_negatives, classification_results.n_null_sims)
end

"""
    calculate_specificity(true_negatives, n_null_sims)

Calculate the specificity (true negative rate) from the number of true negatives
and total number of null simulations.

Specificity measures the proportion of actual negative cases that are correctly
identified as negative. In the context of early warning systems, this represents
the fraction of null simulations where no critical transition occurred and no false
alarm was raised.

# Arguments
- `true_negatives`: Number of true negative detections
- `n_null_sims`: Total number of null simulations

# Returns
- `Float64`: The specificity value as true_negatives / n_null_sims

# Examples
```julia
specificity = calculate_specificity(90, 100)  # Returns 0.9
```
"""
function calculate_specificity(true_negatives, n_null_sims)
    return true_negatives / n_null_sims
end
