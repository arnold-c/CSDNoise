export calculate_auc,
    get_tau

"""
    calculate_auc(emergent_tau, null_tau)

Calculate the Area Under the Curve (AUC) for distinguishing between emergent and null tau values.

This function computes the AUC statistic using the Mann-Whitney U test approach,
which measures how well tau values from emergent scenarios can be distinguished
from tau values from null scenarios. An AUC of 0.5 indicates no discrimination,
while values closer to 1.0 indicate better discrimination ability.

# Arguments
- `emergent_tau`: Vector of tau values from emergent (signal-present) scenarios
- `null_tau`: Vector of tau values from null (signal-absent) scenarios

# Returns
- `Float64`: AUC value between 0 and 1, where:
  - 0.5 = no discrimination ability (random performance)
  - > 0.5 = emergent scenarios tend to have higher tau values
  - < 0.5 = null scenarios tend to have higher tau values
  - 1.0 = perfect discrimination

# Algorithm
The AUC is calculated using the Mann-Whitney U statistic:
1. Combine emergent and null tau vectors
2. Rank all values in descending order (using tied ranks for equal values)
3. Sum ranks corresponding to null scenarios
4. Apply Mann-Whitney U formula to convert to AUC

# Mathematical Formula
```
AUC = (sum_null_ranks - n_null * (n_null + 1) / 2) / (n_emergent * n_null)
```

# Example
```julia
# Tau values from different scenario types
emergent_taus = [0.8, 0.9, 0.7, 0.85]
null_taus = [0.3, 0.4, 0.2, 0.35]

auc = calculate_auc(emergent_taus, null_taus)
# Returns value close to 1.0, indicating good discrimination
```

# Implementation Notes
- Uses `StatsBase.tiedrank` with negative values to rank in descending order
- Handles tied values appropriately through tied ranking
- Equivalent to the probability that a randomly chosen emergent tau exceeds a randomly chosen null tau

# See Also
- [`get_tau`](@ref): Extracts tau values from EWS metrics for AUC calculation
- [`StatsBase.tiedrank`](@ref): Ranking function used in the calculation
"""
function calculate_auc(
        emergent_tau,
        null_tau,
    )
    combined_taus = vcat(null_tau, emergent_tau)

    ranks = StatsBase.tiedrank(-combined_taus)
    n_emergent = length(emergent_tau)
    n_null = length(null_tau)
    sum_null_ranks = sum(ranks[1:n_null])
    return (sum_null_ranks - n_null * (n_null + 1) / 2) /
        (n_emergent * n_null)
end

"""
    get_tau(ews_metrics; tau_metric=:variance_tau, statistic_function=StatsBase.mean)

Extract and summarize tau values from EWS metrics using a specified statistic.

This function retrieves a specific tau metric from an EWS metrics object and
applies a statistical summary function to obtain a single representative value.
Tau metrics represent the strength or magnitude of early warning signals.

# Arguments
- `ews_metrics`: EWS metrics object containing various tau measurements
- `tau_metric::Symbol=:variance_tau`: The specific tau metric to extract
- `statistic_function::Function=StatsBase.mean`: Function to summarize the tau vector

# Common Tau Metrics
The following tau metrics are typically available:
- `:variance_tau`: Tau values based on variance trends
- `:autocovariance_tau`: Tau values based on autocovariance trends
- `:autocorrelation_tau`: Tau values based on autocorrelation trends
- `:coefficient_of_variation_tau`: Tau values based on coefficient of variation trends
- `:skewness_tau`: Tau values based on skewness trends
- `:kurtosis_tau`: Tau values based on kurtosis trends

# Common Statistic Functions
- `StatsBase.mean`: Average tau value (default)
- `StatsBase.median`: Median tau value
- `maximum`: Maximum tau value
- `minimum`: Minimum tau value
- `StatsBase.std`: Standard deviation of tau values

# Returns
- Scalar value: Result of applying the statistic function to the tau vector

# Example
```julia
# Get mean variance tau
mean_var_tau = get_tau(ews_metrics, tau_metric=:variance_tau)

# Get median autocovariance tau
median_autocov_tau = get_tau(
    ews_metrics,
    tau_metric=:autocovariance_tau,
    statistic_function=StatsBase.median
)

# Get maximum coefficient of variation tau
max_cv_tau = get_tau(
    ews_metrics,
    tau_metric=:coefficient_of_variation_tau,
    statistic_function=maximum
)
```

# Implementation Notes
- Uses `getproperty` to dynamically access the specified tau metric
- Flexible design allows any summary statistic function
- Commonly used with `calculate_auc` for discrimination analysis

# See Also
- [`calculate_auc`](@ref): Uses tau values from this function for AUC calculation
- [`StatsBase.mean`](@ref), [`StatsBase.median`](@ref): Common statistic functions
"""
function get_tau(
        ews_metrics;
        tau_metric = :variance_tau,
        statistic_function = StatsBase.mean,
    )
    tau_vector = getproperty(ews_metrics, tau_metric)

    return statistic_function(tau_vector)
end
