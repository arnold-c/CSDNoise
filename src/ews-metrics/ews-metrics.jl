using StatsBase: StatsBase
using LightSumTypes: variant
using DataFrames: DataFrames
using Bumper: @no_escape, @alloc
using Printf: @sprintf
using Dates: Dates

export EWSMetrics

"""
    EWSMetrics(ews_spec::EWSMetricSpecification, timeseries)

Compute Early Warning Signal (EWS) metrics from a time series using the specified configuration.

This function calculates various statistical metrics that serve as early warning signals
for critical transitions in dynamical systems. The metrics include mean, variance,
coefficient of variation, index of dispersion, skewness, kurtosis, autocovariance,
and autocorrelation, along with their Kendall tau-B correlation coefficients with time.

# Arguments
- `ews_spec::EWSMetricSpecification`: Configuration specifying the EWS method, bandwidth, aggregation, and lag parameters
- `timeseries`: Input time series data to analyze

# Returns
- `EWSMetrics`: Struct containing all computed EWS metrics and their temporal correlations

# Throws
- `ErrorException`: If the aggregated time series length is insufficient for the specified bandwidth

# Examples
```julia
# Create EWS specification
spec = EWSMetricSpecification(Backward(), 50, Day(1), 1)

# Compute EWS metrics from time series
metrics = EWSMetrics(spec, my_timeseries)
```
"""
function EWSMetrics(
        ews_spec::EWSMetricSpecification, timeseries
    )
    aggregated_timeseries = aggregate_timeseries(
        timeseries, ews_spec.aggregation
    )

    aggregated_bandwidth = aggregate_bandwidth(ews_spec)

    if length(aggregated_timeseries) < aggregated_bandwidth
        error(
            "Not enough data for bandwidth: bandwidth = $(ews_spec.bandwidth), aggregated bandwidth = $(aggregated_bandwidth), aggregation = $(ews_spec.aggregation), aggregated time series length = $(length(aggregated_timeseries))\n",
            "ews_specification = $(ews_spec)",
        )
    end

    mean_vec = ews_mean(
        ews_spec.method, aggregated_timeseries, aggregated_bandwidth
    )
    var_vec = ews_var(
        ews_spec.method, mean_vec, aggregated_timeseries, aggregated_bandwidth
    )
    var2_vec = var_vec .^ 2
    sd_vec = sqrt.(var_vec)
    sd3_vec = sd_vec .^ 3
    m3_vec = _ews_moment(
        ews_spec.method, mean_vec, aggregated_timeseries, 3,
        aggregated_bandwidth,
    )
    m4_vec = _ews_moment(
        ews_spec.method, mean_vec, aggregated_timeseries, 4,
        aggregated_bandwidth,
    )
    autocov_vec = ews_autocov(
        ews_spec.method,
        mean_vec,
        aggregated_timeseries,
        aggregated_bandwidth;
        lag = ews_spec.lag,
    )

    cov_vec = ews_cov(sd_vec, mean_vec)
    iod_vec = ews_iod(var_vec, mean_vec)
    skew_vec = ews_skew(m3_vec, sd3_vec)
    kurtosis_vec = ews_kurtosis(m4_vec, var2_vec)
    autocor_vec = ews_autocor(autocov_vec, sd_vec)

    return EWSMetrics(
        ews_spec,
        mean_vec,
        var_vec,
        cov_vec,
        iod_vec,
        skew_vec,
        kurtosis_vec,
        autocov_vec,
        autocor_vec,
        kendall_tau(mean_vec),
        kendall_tau(var_vec),
        kendall_tau(cov_vec),
        kendall_tau(iod_vec),
        kendall_tau(skew_vec),
        kendall_tau(kurtosis_vec),
        kendall_tau(autocov_vec),
        kendall_tau(autocor_vec),
    )
end

"""
    ews_mean(method::EWSMethod, timeseries, bandwidth::Integer)

Calculate rolling mean of a time series using the specified EWS method.

# Arguments
- `method::EWSMethod`: EWS calculation method (Centered or Backward)
- `timeseries`: Input time series data
- `bandwidth::Integer`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling mean values for each time point

# Notes
- For Centered method: uses symmetric window around each point
- For Backward method: uses only past values
"""
function ews_mean(
        method::EWSMethod,
        timeseries,
        bandwidth::T1,
    ) where {T1 <: Integer}
    mean_vec = zeros(Float64, length(timeseries))
    ews_mean!(mean_vec, method, timeseries, bandwidth)
    return mean_vec
end

"""
    ews_mean!(mean_vec, method::EWSMethod, timeseries, bandwidth::Integer)

In-place calculation of rolling mean, storing results in pre-allocated vector.

This mutating version avoids memory allocation by writing results directly to the
provided `mean_vec` array, making it more efficient for repeated calculations.

# Arguments
- `mean_vec`: Pre-allocated output vector to store mean values
- `method::EWSMethod`: EWS calculation method (Centered or Backward)
- `timeseries`: Input time series data
- `bandwidth::Integer`: Window size for rolling calculation

# Returns
- `nothing` (results stored in `mean_vec`)

# Notes
- `mean_vec` must be pre-allocated with correct length
"""
function ews_mean!(
        mean_vec,
        method::EWSMethod,
        timeseries,
        bandwidth::T1,
    ) where {T1 <: Integer}
    mean_func!(mean_vec, method, timeseries, bandwidth)
    return nothing
end

"""
    mean_func!(mean_vec, method::EWSMethod, timeseries, bandwidth)

Dispatch function for method-specific mean calculations.

This function extracts the variant from the EWSMethod sum type and dispatches
to the appropriate implementation (Centered or Backward).

# Arguments
- `mean_vec`: Pre-allocated output vector
- `method::EWSMethod`: EWS method containing the calculation variant
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `nothing` (results stored in `mean_vec`)
"""
function mean_func!(
        mean_vec,
        method::EWSMethod,
        timeseries,
        bandwidth,
    )
    mean_func!(mean_vec, variant(method), timeseries, bandwidth)
    return nothing
end

"""
    mean_func!(mean_vec, method::Centered, timeseries, bandwidth)

Calculate rolling mean using centered windows around each time point.

This implementation uses symmetric windows that extend both forward and backward
from each point. At boundaries, it adapts the window size to use available data:
- Near start: extends window from the start to the index + bandwidth - 1
- Near end: extends window from the end to the index - bandwidth + 1
- At extremes: uses entire available time series

# Arguments
- `mean_vec`: Pre-allocated output vector
- `method::Centered`: Centered windowing method
- `timeseries`: Input time series data
- `bandwidth`: Half-width of the symmetric window

# Returns
- `nothing` (results stored in `mean_vec`)

# Notes
- Window boundaries are handled adaptively to maximize data usage
- Provides symmetric treatment around each point when possible
"""
function mean_func!(
        mean_vec,
        method::Centered,
        timeseries,
        bandwidth,
    )
    tlength = length(timeseries)
    @inbounds for i in eachindex(timeseries)
        if i < bandwidth && i + bandwidth <= tlength
            mean_vec[i] = mean(@view(timeseries[begin:(i + bandwidth - 1)]))
        elseif i < bandwidth
            mean_vec[i] = mean(@view(timeseries[begin:end]))
        elseif i + bandwidth > tlength
            mean_vec[i] = mean(@view(timeseries[(i - bandwidth + 1):end]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bandwidth + 1):(i + bandwidth - 1)]))
        end
    end
    return nothing
end

"""
    mean_func!(mean_vec, method::Backward, timeseries, bandwidth)

Calculate rolling mean using backward-looking windows.

This implementation only uses past and current values, making it suitable for
real-time applications where future data is not available. The window extends
backward from each point:
- Early points: uses all available past data (expanding window)
- Later points: uses fixed bandwidth of past values

# Arguments
- `mean_vec`: Pre-allocated output vector
- `method::Backward`: Backward-looking windowing method
- `timeseries`: Input time series data
- `bandwidth`: Number of past time points to include in window

# Returns
- `nothing` (results stored in `mean_vec`)

# Notes
- Early time points use expanding windows with all available history
- For longer time series, produces functionally equivalent values as backwards
windowing methods that have offset indices
"""
function mean_func!(
        mean_vec,
        method::Backward,
        timeseries,
        bandwidth,
    )
    @inbounds for i in eachindex(timeseries)
        if i < bandwidth
            mean_vec[i] = mean(@view(timeseries[begin:i]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bandwidth + 1):i]))
        end
    end
    return nothing
end

"""
    _ews_moment(method::EWSMethod, timeseries, moment, bandwidth)

Calculate rolling central moments by first computing the rolling mean.

This variant computes the rolling mean internally and then calculates the
specified central moment. Less efficient when mean is already available.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `moment`: Order of the central moment (e.g., 2 for variance, 3 for skewness)
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling central moment values

# Notes
- Computes mean internally, so use the other variant if mean is pre-computed
- Central moments are calculated as E[(X - μ)^k] where k is the moment order
"""
function _ews_moment(
        method::EWSMethod, timeseries, moment, bandwidth
    )
    return _ews_moment(
        method,
        ews_mean(method, timeseries, bandwidth),
        timeseries,
        moment,
        bandwidth,
    )
end

"""
    _ews_moment(method::EWSMethod, mean_timeseries, timeseries, moment, bandwidth)

Calculate rolling central moments using pre-computed rolling mean values.

This more efficient variant uses pre-computed mean values to calculate central
moments, avoiding redundant mean calculations. Uses memory-efficient allocation
with the Bumper.jl `@no_escape` macro.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `mean_timeseries`: Pre-computed rolling mean values
- `timeseries`: Input time series data
- `moment`: Order of the central moment (e.g., 2 for variance, 3 for skewness)
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling central moment values

# Notes
- More efficient when mean is already computed
- Uses stack allocation for temporary arrays via `@no_escape`
- Central moments: E[(X - μ)^k] where μ is the rolling mean
"""
function _ews_moment(
        method::EWSMethod, mean_timeseries, timeseries, moment, bandwidth
    )
    return @no_escape begin
        diff = @alloc(Float64, length(timeseries))
        diff .= (timeseries .- mean_timeseries) .^ moment
        ews_mean(method, diff, bandwidth)
    end
end

"""
    ews_var(method::EWSMethod, timeseries, bandwidth)

Calculate rolling variance of a time series.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling variance values

# Notes
- Variance is calculated as the second central moment
"""
function ews_var(method::EWSMethod, timeseries, bandwidth)
    mean_vec = ews_mean(method, timeseries, bandwidth)
    return _ews_moment(method, mean_vec, timeseries, 2, bandwidth)
end

"""
    ews_var(method::EWSMethod, mean_vec, timeseries, bandwidth)

Calculate rolling variance using pre-computed mean values.

More efficient variant when mean values are already available.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `mean_vec`: Pre-computed rolling mean values
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling variance values
"""
function ews_var(method::EWSMethod, mean_vec, timeseries, bandwidth)
    return _ews_moment(method, mean_vec, timeseries, 2, bandwidth)
end

"""
    ews_cov(method::EWSMethod, timeseries, bandwidth)

Calculate rolling coefficient of variation by computing mean and variance internally.

This variant performs the complete calculation from raw time series data,
computing rolling mean, variance, and standard deviation before calculating
the coefficient of variation. Uses memory-efficient stack allocation.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling coefficient of variation values (σ/μ)

# Notes
- Computes all intermediate statistics internally
- Less efficient than the variant using pre-computed values
- Uses `@no_escape` for memory-efficient temporary allocations
"""
function ews_cov(method::EWSMethod, timeseries, bandwidth)
    return @no_escape begin
        mean_vec = @alloc(Float64, length(timeseries))
        var_vec = @alloc(Float64, length(timeseries))
        mean_vec = ews_mean(method, timeseries, bandwidth)
        var_vec = ews_var(method, mean_vec, timeseries, bandwidth)
        sd_vec = @alloc(Float64, length(var_vec))
        sd_vec .= sqrt.(var_vec)
        ews_cov(sd_vec, mean_vec)
    end
end

"""
    ews_cov(sd_vec, mean_vec)

Calculate coefficient of variation from standard deviation and mean vectors.

The coefficient of variation is a normalized measure of dispersion, calculated as
the ratio of standard deviation to mean.

# Arguments
- `sd_vec`: Vector of standard deviation values
- `mean_vec`: Vector of mean values

# Returns
- `Vector{Float64}`: Coefficient of variation values (σ/μ)

# Notes
- Higher values indicate increased relative variability
- Important EWS metric for detecting critical transitions
"""
function ews_cov(sd_vec, mean_vec)
    return sd_vec ./ mean_vec
end

"""
    ews_iod(method::EWSMethod, timeseries, bandwidth)

Calculate rolling index of dispersion by computing variance and mean internally.

This variant performs the complete calculation from raw time series data,
computing both rolling variance and mean before calculating their ratio.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling index of dispersion values (σ²/μ)

# Notes
- Computes variance and mean internally
- Less efficient than the variant using pre-computed values
- Useful when only raw time series data is available
"""
function ews_iod(method::EWSMethod, timeseries, bandwidth)
    var_vec = ews_var(method, timeseries, bandwidth)
    mean_vec = ews_mean(method, timeseries, bandwidth)
    return ews_iod(var_vec, mean_vec)
end

"""
    ews_iod(var_vec, mean_vec)

Calculate index of dispersion from variance and mean vectors.

The index of dispersion is the ratio of variance to mean, useful for detecting
overdispersion in count data.

# Arguments
- `var_vec`: Vector of variance values
- `mean_vec`: Vector of mean values

# Returns
- `Vector{Float64}`: Index of dispersion values (σ²/μ)

# Notes
- Values > 1 indicate overdispersion
- Values < 1 indicate underdispersion
- Values = 1 indicate Poisson-like behavior
"""
function ews_iod(var_vec, mean_vec)
    return var_vec ./ mean_vec
end

"""
    ews_skew(method::EWSMethod, timeseries, bandwidth)

Calculate rolling skewness by computing all required statistics internally.

This variant performs the complete skewness calculation from raw time series,
computing rolling mean, third central moment, variance, and standard deviation.
Uses a single worker vector for memory efficiency.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling skewness values (μ₃/σ³)

# Notes
- Computes mean, variance, and third moment internally
- Uses memory-efficient allocation with shared worker vector
- Skewness = E[(X-μ)³]/σ³ where μ is mean and σ is standard deviation
"""
function ews_skew(method::EWSMethod, timeseries, bandwidth)
    return @no_escape begin
        # use for both mean and sd vec as only one needed at a time
        worker_vec = @alloc(Float64, length(timeseries))
        m3_vec = @alloc(Float64, length(timeseries))
        sd3_vec = @alloc(Float64, length(timeseries))
        worker_vec .= ews_mean(method, timeseries, bandwidth)
        m3_vec .= _ews_moment(method, worker_vec, timeseries, 3, bandwidth)
        sd3_vec .=
            sqrt.(ews_var(method, worker_vec, timeseries, bandwidth)) .^ 3
        ews_skew(m3_vec, sd3_vec)
    end
end

"""
    ews_skew(m3_vec, sd3_vec)

Calculate skewness from third moment and cubed standard deviation vectors.

Skewness measures the asymmetry of the probability distribution.

# Arguments
- `m3_vec`: Vector of third central moments
- `sd3_vec`: Vector of cubed standard deviations

# Returns
- `Vector{Float64}`: Skewness values (μ₃/σ³)

# Notes
- Positive values indicate right-skewed distributions
- Negative values indicate left-skewed distributions
- Zero indicates symmetric distributions
"""
function ews_skew(m3_vec, sd3_vec)
    return m3_vec ./ sd3_vec
end

"""
    ews_kurtosis(method::EWSMethod, timeseries, bandwidth)

Calculate rolling kurtosis by computing all required statistics internally.

This variant performs the complete kurtosis calculation from raw time series,
computing rolling mean, fourth central moment, and variance before calculating
the kurtosis ratio.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation

# Returns
- `Vector{Float64}`: Rolling kurtosis values (μ₄/σ⁴)

# Notes
- Computes mean, variance, and fourth moment internally
- Uses memory-efficient stack allocation for temporary arrays
- Kurtosis = E[(X-μ)⁴]/σ⁴ where μ is mean and σ² is variance
"""
function ews_kurtosis(method::EWSMethod, timeseries, bandwidth)
    return @no_escape begin
        mean_vec = @alloc(Float64, length(timeseries))
        m4_vec = @alloc(Float64, length(timeseries))
        var2_vec = @alloc(Float64, length(timeseries))
        mean_vec .= ews_mean(method, timeseries, bandwidth)
        m4_vec .= _ews_moment(method, mean_vec, timeseries, 4, bandwidth)
        var2_vec .= ews_var(method, mean_vec, timeseries, bandwidth) .^ 2
        ews_kurtosis(m4_vec, var2_vec)
    end
end

"""
    ews_kurtosis(m4_vec, var2_vec)

Calculate kurtosis from fourth moment and squared variance vectors.

Kurtosis measures the "tailedness" of the probability distribution.

# Arguments
- `m4_vec`: Vector of fourth central moments
- `var2_vec`: Vector of squared variances

# Returns
- `Vector{Float64}`: Kurtosis values (μ₄/σ⁴)

# Notes
- Values > 3 indicate heavy-tailed distributions (leptokurtic)
- Values < 3 indicate light-tailed distributions (platykurtic)
- Values = 3 indicate normal distribution-like tails (mesokurtic)
"""
function ews_kurtosis(m4_vec, var2_vec)
    return m4_vec ./ var2_vec
end

"""
    ews_autocov(method::EWSMethod, timeseries, bandwidth; lag=1)

Calculate rolling autocovariance by computing the rolling mean internally.

This variant computes the rolling mean from the raw time series before
calculating autocovariance. Less efficient when mean is already available.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation
- `lag`: Time lag for autocovariance calculation (default: 1)

# Returns
- `Vector{Float64}`: Rolling autocovariance values

# Notes
- Computes rolling mean internally using stack allocation
- Less efficient than variant with pre-computed mean
- First `lag` values are set to NaN due to insufficient data
"""
function ews_autocov(method::EWSMethod, timeseries, bandwidth; lag = 1)
    return @no_escape begin
        mean_vec = @alloc(Float64, length(timeseries))
        mean_vec .= ews_mean(method, timeseries, bandwidth)
        ews_autocov(
            method, mean_vec, timeseries, bandwidth; lag = lag
        )
    end
end

"""
    ews_autocov(method::EWSMethod, mean_vec, timeseries, bandwidth; lag=1)

Calculate rolling autocovariance of a time series.

Autocovariance measures the covariance between a time series and a lagged version of itself.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `mean_vec`: Pre-computed rolling mean values
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation
- `lag`: Time lag for autocovariance calculation (default: 1)

# Returns
- `Vector{Float64}`: Rolling autocovariance values

# Notes
- First `lag` values are set to NaN due to insufficient data
- Higher values indicate stronger temporal correlation
"""
function ews_autocov(
        method::EWSMethod,
        mean_vec,
        timeseries,
        bandwidth;
        lag = 1,
    )
    meandiff = timeseries .- mean_vec
    autocov_vec = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i <= lag
            continue
        end
        autocov_vec[i] = meandiff[i] * meandiff[i - lag]
    end
    autocov_vec = ews_mean(method, autocov_vec, bandwidth)
    autocov_vec[begin:lag] .= NaN
    return autocov_vec
end

"""
    ews_autocor(method::EWSMethod, timeseries, bandwidth; lag=1)

Calculate rolling autocorrelation by computing all required statistics internally.

This variant performs the complete autocorrelation calculation from raw time series,
computing variance, standard deviation, autocovariance, and the lagged standard
deviation product needed for normalization.

# Arguments
- `method::EWSMethod`: EWS calculation method
- `timeseries`: Input time series data
- `bandwidth`: Window size for rolling calculation
- `lag`: Time lag for autocorrelation calculation (default: 1)

# Returns
- `Vector{Float64}`: Rolling autocorrelation values (range: [-1, 1])

# Notes
- Computes variance, autocovariance, and lagged products internally
- Uses memory-efficient stack allocation for temporary arrays
- Autocorrelation = autocovariance / (σₜ × σₜ₋ₗₐ𝓰)
"""
function ews_autocor(method::EWSMethod, timeseries, bandwidth; lag = 1)
    sd2 = zeros(Float64, length(timeseries))
    @no_escape begin
        var_vec = @alloc(Float64, length(timeseries))
        sd_vec = @alloc(Float64, length(timeseries))
        autocov_vec = @alloc(Float64, length(timeseries))
        var_vec .= ews_var(method, timeseries, bandwidth)
        sd_vec .= sqrt.(var_vec)
        autocov_vec .= ews_autocov(method, timeseries, bandwidth; lag = lag)
        lagged_sd = @alloc(Float64, length(sd_vec))
        _lagged_vector(lagged_sd, sd_vec, lag)
        sd2 .= sd_vec .* lagged_sd
    end
    return autocov_vec ./ sd2
end

"""
    ews_autocor(autocov_vec, sd_vec; lag=1)

Calculate autocorrelation from autocovariance and standard deviation vectors.

Autocorrelation is the normalized autocovariance, providing a scale-free measure
of temporal correlation.

# Arguments
- `autocov_vec`: Vector of autocovariance values
- `sd_vec`: Vector of standard deviation values
- `lag`: Time lag used in autocovariance calculation (default: 1)

# Returns
- `Vector{Float64}`: Autocorrelation values (range: [-1, 1])

# Notes
- Values close to 1 indicate strong positive temporal correlation
- Values close to -1 indicate strong negative temporal correlation
- Values close to 0 indicate weak temporal correlation
"""
function ews_autocor(autocov_vec, sd_vec; lag = 1)
    sd2 = zeros(Float64, length(sd_vec))
    @no_escape begin
        lagged_sd = @alloc(Float64, length(sd_vec))
        _lagged_vector(lagged_sd, sd_vec, lag)
        sd2 .= sd_vec .* lagged_sd
    end
    return autocov_vec ./ sd2
end

"""
    _lagged_vector(lagged_vec, vec, lag)

Create a lagged version of a vector by shifting values backward in time.

This helper function creates a time-lagged version of the input vector,
where each element at position i contains the value from position i-lag.
Elements without sufficient history are set to NaN.

# Arguments
- `lagged_vec`: Pre-allocated output vector to store lagged values
- `vec`: Input vector to lag
- `lag`: Number of time steps to lag (must be positive)

# Returns
- `nothing` (results stored in `lagged_vec`)

# Notes
- First `lag` elements are set to NaN due to insufficient history
- Used for autocorrelation calculations requiring lagged standard deviations
- In-place operation for memory efficiency
"""
function _lagged_vector(lagged_vec, vec, lag)
    @inbounds for i in eachindex(lagged_vec)
        if i <= lag
            lagged_vec[i] = NaN
            continue
        end
        lagged_vec[i] = vec[i - lag]
    end
    return nothing
end

"""
    kendall_tau(ews_vec::Vector{F}) where {F <: AbstractFloat}

Calculate Kendall's tau-B correlation coefficient between an EWS metric and time.

This function computes the rank correlation between the EWS metric values and their
temporal ordering, providing a measure of monotonic trend over time.

# Arguments
- `ews_vec::Vector{F}`: Vector of EWS metric values

# Returns
- `Float64`: Kendall's tau-B correlation coefficient (range: [-1, 1])

# Notes
- Positive values indicate increasing trend over time
- Negative values indicate decreasing trend over time
- NaN values are filtered out before calculation
- Used to detect temporal trends in EWS metrics approaching critical transitions
"""
function kendall_tau(ews_vec::Vector{F}) where {F <: AbstractFloat}
    filtered_ews_vec = filter(x -> !isnan(x), ews_vec)
    return StatsBase.corkendall(
        collect(1:length(filtered_ews_vec)), filtered_ews_vec
    )
end
