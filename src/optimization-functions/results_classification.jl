using UnPack: @unpack
using Try: Try

export calculate_ews_classification_results,
    calculate_sensitivity,
    calculate_specificity,
    calculate_balanced_accuracy


function calculate_ews_classification_results(
        scenario,
        ews_metrics,
        null_ews_metrics,
    )
    @unpack ews_metric,
        ews_threshold_window,
        ensemble_specification,
        threshold_quantile,
        consecutive_thresholds = scenario
    @unpack burnin = ensemble_specification.time_parameters

    ews_metric_symbol = Symbol(ews_metric)

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
