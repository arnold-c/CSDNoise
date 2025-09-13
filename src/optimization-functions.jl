export calculate_sensitivity, calculate_specificity, calculate_balanced_accuracy

"""
    calculate_balanced_accuracy(sensitivity, specificity)

Calculate the balanced accuracy from sensitivity and specificity values.

Balanced accuracy is the arithmetic mean of sensitivity (true positive rate) and
specificity (true negative rate), providing a single metric that accounts for
performance on both positive and negative cases.

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
