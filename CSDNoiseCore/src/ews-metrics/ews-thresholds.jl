export exceeds_ews_threshold,
    get_ews_metric_vec

"""
    exceeds_ews_threshold(ewsmetrics, metric, window_type, quantile=0.95, burn_in=Dates.Day(10))

Determine whether EWS metric values exceed dynamically calculated thresholds.

This function evaluates whether early warning signal (EWS) metric values exceed
thresholds calculated from historical data using either expanding or rolling windows.
The threshold calculation method is determined by the `window_type` parameter.

# Arguments
- `ewsmetrics::EWSMetrics`: Container with calculated EWS metrics
- `metric::Symbol`: The specific EWS metric to evaluate (e.g., `:variance`, `:autocovariance`)
- `window_type::EWSThresholdWindowType`: Type of threshold window calculation
- `quantile::Float64=0.95`: Quantile level for threshold calculation (default: 0.95)
- `burn_in::Dates.Period=Dates.Day(10)`: Initial period to exclude from threshold calculation

# Returns
- `Vector{Bool}`: Boolean vector indicating whether each time point exceeds the threshold

# Example
```julia
# Using expanding window thresholds
exceeds = exceeds_ews_threshold(
    ews_metrics,
    :variance,
    EWSThresholdWindowType(ExpandingThresholdWindow()),
    0.95,
    Dates.Day(14)
)
```

# See Also
- [`get_ews_metric_vec`](@ref): Extracts specific metric vectors from EWSMetrics
- [`ExpandingThresholdWindow`](@ref): Expanding window threshold calculation
- [`RollingThresholdWindow`](@ref): Rolling window threshold calculation
"""
function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::EWSThresholdWindowType,
        quantile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    return exceeds_ews_threshold(
        ewsmetrics,
        metric,
        LightSumTypes.variant(window_type),
        quantile,
        burn_in,
    )
end

"""
    exceeds_ews_threshold(ewsmetrics, metric, window_type::RollingThresholdWindow, quantile=0.95, burn_in=Dates.Day(10))

Rolling window threshold calculation for EWS metrics (placeholder implementation).

This method is intended to calculate thresholds using a rolling window approach,
where thresholds are computed from a fixed-size moving window of historical data.
Currently returns a placeholder result for type stability.

# Arguments
- `ewsmetrics::EWSMetrics`: Container with calculated EWS metrics
- `metric::Symbol`: The specific EWS metric to evaluate
- `window_type::RollingThresholdWindow`: Rolling window threshold specification
- `quantile::Float64=0.95`: Quantile level for threshold calculation
- `burn_in::Dates.Period=Dates.Day(10)`: Initial period to exclude

# Returns
- `Vector{Bool}`: Placeholder boolean vector (currently `fill(false, 2)`)

# Implementation Status
This method is not yet implemented and serves as a placeholder for type stability.
The actual rolling window threshold calculation logic needs to be added.

# See Also
- [`exceeds_ews_threshold`](@ref): Main threshold evaluation function
- [`ExpandingThresholdWindow`](@ref): Alternative expanding window approach
"""
function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::RollingThresholdWindow,
        quantile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    return fill(false, 2)
end

"""
    exceeds_ews_threshold(ewsmetrics, metric, window_type::ExpandingThresholdWindow, quantile=0.95, burn_in=Dates.Day(10))

Calculate threshold exceedances using expanding window approach.

This method implements expanding window threshold calculation, where thresholds
are computed from all historical data up to each time point. The threshold at
each time point is the specified quantile of all previous observations after
the burn-in period.

# Arguments
- `ewsmetrics::EWSMetrics`: Container with calculated EWS metrics
- `metric::Symbol`: The specific EWS metric to evaluate
- `window_type::ExpandingThresholdWindow`: Expanding window threshold specification
- `quantile::Float64=0.95`: Quantile level for threshold calculation
- `burn_in::Dates.Period=Dates.Day(10)`: Initial period to exclude from threshold calculation

# Returns
- `Vector{Bool}`: Boolean vector indicating threshold exceedances for each time point

# Implementation Details
- Converts burn-in period to index based on aggregation interval
- Uses expanding quantile calculation for threshold determination
- Handles NaN values in the EWS metric vector appropriately
- Compares current values against previous time point's threshold

# Example
```julia
exceeds = exceeds_ews_threshold(
    ews_metrics,
    :autocovariance,
    ExpandingThresholdWindow(),
    0.90,
    Dates.Day(7)
)
```

# See Also
- [`_expanding_ews_thresholds`](@ref): Internal implementation function
- [`get_ews_metric_vec`](@ref): Extracts metric vectors from EWSMetrics
"""
function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::ExpandingThresholdWindow,
        quantile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    ews_vec = get_ews_metric_vec(ewsmetrics, metric)

    @unpack aggregation = ewsmetrics.ews_specification
    burn_in_index = Int64(Dates.days(burn_in) ÷ Dates.days(aggregation))

    @assert burn_in_index >= 1 && burn_in_index <= length(ews_vec)

    return _expanding_ews_thresholds(
        ews_vec,
        quantile,
        burn_in_index,
    )
end

"""
    get_ews_metric_vec(ewsmetrics, metric)

Extract a specific EWS metric vector from an EWSMetrics container.

This function retrieves the time series vector for a specified early warning
signal metric from an `EWSMetrics` object. It validates that the requested
metric is one of the supported EWS metrics and returns the corresponding
vector of Float64 values.

# Arguments
- `ewsmetrics::EWSMetrics`: Container with calculated EWS metrics
- `metric::Symbol`: The specific metric to extract

# Supported Metrics
The following metrics can be extracted:
- `:mean`: Mean of the time series
- `:variance`: Variance of the time series
- `:coefficient_of_variation`: Coefficient of variation
- `:index_of_dispersion`: Index of dispersion (variance-to-mean ratio)
- `:skewness`: Third moment (skewness)
- `:kurtosis`: Fourth moment (kurtosis)
- `:autocovariance`: Lag-1 autocovariance
- `:autocorrelation`: Lag-1 autocorrelation

# Returns
- `Vector{Float64}`: Time series vector of the requested metric

# Throws
- `AssertionError`: If the requested metric is not in the supported list

# Example
```julia
variance_vec = get_ews_metric_vec(ews_metrics, :variance)
autocorr_vec = get_ews_metric_vec(ews_metrics, :autocorrelation)
```

# See Also
- [`EWSMetrics`](@ref): Container type for EWS metric calculations
- [`exceeds_ews_threshold`](@ref): Uses this function to extract metrics for threshold evaluation
"""
function get_ews_metric_vec(
        ewsmetrics::T1,
        metric::T2,
    ) where {T1 <: EWSMetrics, T2 <: Symbol}
    @assert metric in [
        :mean,
        :variance,
        :coefficient_of_variation,
        :index_of_dispersion,
        :skewness,
        :kurtosis,
        :autocovariance,
        :autocorrelation,
    ]

    ews_vec = getproperty(ewsmetrics, metric)::Vector{Float64}
    return ews_vec
end


"""
    _expanding_ews_thresholds(ews_vec, quantile, burn_in_index)

Internal function to calculate expanding window threshold exceedances.

This function implements the core logic for expanding window threshold calculation.
It iteratively builds up the distribution of historical values and determines
threshold exceedances based on the specified quantile level.

# Arguments
- `ews_vec::Vector{Float64}`: Vector of EWS metric values
- `quantile::Float64`: Quantile level for threshold calculation (e.g., 0.95)
- `burn_in_index::Int64`: Index after which to start threshold calculations

# Returns
- `Vector{Bool}`: Boolean vector indicating threshold exceedances

# Implementation Details
- Handles NaN values by skipping them in calculations but preserving array structure
- Uses a worker vector to avoid computing quantiles on arrays containing NaNs
- Calculates expanding quantiles using all available historical data
- Compares current values against the threshold from the previous time point

# Algorithm
1. Initialize result vectors and working arrays
2. For each time point after burn-in:
   - Add non-NaN values to the working distribution
   - Calculate the expanding quantile threshold
   - Compare current value against previous threshold
3. Return boolean vector of exceedances

# See Also
- [`exceeds_ews_threshold`](@ref): Public interface that calls this function
- [`StatsBase.quantile`](@ref): Used for threshold calculation
"""
function _expanding_ews_thresholds(
        ews_vec::Vector{Float64},
        quantile::Float64,
        burn_in_index::Int64,
    )
    ews_vec_len = length(ews_vec)
    ews_distributions = fill(NaN, ews_vec_len)
    ews_worker_vec = fill(
        NaN, sum((!isnan).(ews_vec))
    )
    exceeds_thresholds = zeros(Bool, ews_vec_len)

    worker_ind = 0
    for i in eachindex(ews_vec)
        if isnan(ews_vec[i])
            if i > burn_in_index
                ews_distributions[i] = ews_distributions[(i - 1)]
            end
            continue
        end
        worker_ind += 1
        # use online stats to build up new distribution to avoid computing quantiles for vectors containing NaNs
        ews_worker_vec[worker_ind] = ews_vec[i]
        if i > burn_in_index
            ews_distributions[i] = StatsBase.quantile(
                @view(ews_worker_vec[1:worker_ind]),
                quantile
            )

            exceeds_thresholds[i] = ews_vec[i] >= ews_distributions[(i - 1)]
        end
    end

    @assert worker_ind == length(ews_worker_vec)

    return exceeds_thresholds
end
