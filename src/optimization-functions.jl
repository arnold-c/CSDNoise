export calculate_sensitivity, calculate_specificity, calculate_accuracy

"""
    calculate_accuracy(sensitivity, specificity)

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
accuracy = calculate_accuracy(0.8, 0.9)  # Returns 0.85
```
"""
function calculate_accuracy(sensitivity, specificity)
    return (sensitivity + specificity) / 2
end

"""
    calculate_sensitivity(true_positives, ensemble_nsims)

Calculate the sensitivity (true positive rate) from the number of true positives 
and total number of simulations.

Sensitivity measures the proportion of actual positive cases that are correctly 
identified as positive. In the context of early warning systems, this represents 
the fraction of simulations where a critical transition was correctly detected.

# Arguments
- `true_positives`: Number of true positive detections
- `ensemble_nsims`: Total number of simulations in the ensemble

# Returns
- `Float64`: The sensitivity value as true_positives / ensemble_nsims

# Examples
```julia
sensitivity = calculate_sensitivity(80, 100)  # Returns 0.8
```
"""
function calculate_sensitivity(true_positives, ensemble_nsims)
    return true_positives / ensemble_nsims

end

"""
    calculate_specificity(true_negatives, ensemble_nsims)

Calculate the specificity (true negative rate) from the number of true negatives 
and total number of simulations.

Specificity measures the proportion of actual negative cases that are correctly 
identified as negative. In the context of early warning systems, this represents 
the fraction of simulations where no critical transition occurred and no false 
alarm was raised.

# Arguments
- `true_negatives`: Number of true negative detections
- `ensemble_nsims`: Total number of simulations in the ensemble

# Returns
- `Float64`: The specificity value as true_negatives / ensemble_nsims

# Examples
```julia
specificity = calculate_specificity(90, 100)  # Returns 0.9
```
"""
function calculate_specificity(true_negatives, ensemble_nsims)
    return true_negatives / ensemble_nsims
end
