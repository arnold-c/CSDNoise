export EWSMethod,
    Backward,
    Centered,
    EWSMetricSpecification,
    EWSMetrics,
    EnsembleEWSMetrics,
    EWSThresholdWindowType,
    ExpandingThresholdWindow,
    RollingThresholdWindow,
    EWSEndDateType,
    ReffStart,
    ReffEnd,
    OutbreakStart,
    OutbreakMiddle,
    OutbreakEnd

"""
    AbstractEWSMethod

Abstract base type for Early Warning Signal (EWS) calculation methods.

This type serves as the parent for all EWS calculation approaches, defining the interface
for different windowing strategies used in time series analysis for outbreak detection.
"""
abstract type AbstractEWSMethod end

"""
    Backward

EWS calculation method using backward-looking windows.

This method calculates Early Warning Signals by looking backward from each time point,
using only historical data available up to that point. This approach is more realistic
for real-time outbreak detection scenarios.
"""
struct Backward end

"""
    Centered

EWS calculation method using centered windows.

This method calculates Early Warning Signals using centered windows around each time point,
incorporating both past and future data. While not suitable for real-time detection,
this approach can provide better signal quality for retrospective analysis.
"""
struct Centered end

"""
    EWSMethod

Sum type representing the available Early Warning Signal calculation methods.

This type can be either `Backward` or `Centered`, determining how time windows
are positioned relative to each calculation point in the time series.

# Variants
- `Backward`: Uses backward-looking windows (realistic for real-time detection)
- `Centered`: Uses centered windows (better for retrospective analysis)

# Example
```julia
method = EWSMethod(Backward())
# or
method = EWSMethod(Centered())
```
"""
LightSumTypes.@sumtype EWSMethod(Backward, Centered) <: AbstractEWSMethod

"""
    EWSMetricSpecification

Configuration for Early Warning Signal metric calculations.

This struct defines all parameters needed to compute EWS metrics from time series data,
including the calculation method, temporal aggregation settings, and file system organization.

# Fields
- `method::EWSMethod`: The windowing method (Backward or Centered) for EWS calculations
- `aggregation::Dates.Day`: Time period for aggregating raw data before EWS calculation
- `bandwidth::Dates.Day`: Window size for computing rolling statistics in EWS metrics
- `lag::Int64`: Number of time steps to lag the EWS calculation (for lead time analysis)
- `dirpath::String`: File system path for storing/retrieving cached EWS results

# Constructors
Multiple constructors are available:
- Full specification with explicit `dirpath`
- Automatic path generation from parameters (recommended)
- Support for general `DatePeriod` types that get converted to `Day`

# Example
```julia
# Automatic path generation
spec = EWSMetricSpecification(
    method = EWSMethod(Backward()),
    aggregation = Dates.Day(7),
    bandwidth = Dates.Day(30),
    lag = 0
)

# Using DatePeriod types
spec = EWSMetricSpecification(
    EWSMethod(Centered()),
    Dates.Week(1),      # Converted to Day(7)
    Dates.Month(1),     # Converted to Day(30)
    5
)
```
"""
Base.@kwdef struct EWSMetricSpecification
    method::EWSMethod
    aggregation::Dates.Day
    bandwidth::Dates.Day
    lag::Int64
    dirpath::String
end

"""
    EWSMetricSpecification(method, aggregation, bandwidth, lag)

Constructor for `EWSMetricSpecification` with automatic path generation.

This constructor automatically generates the `dirpath` field based on the provided
parameters, creating a standardized file system organization for EWS results.

# Arguments
- `method::EWSMethod`: The windowing method for EWS calculations
- `aggregation::Dates.Day`: Time period for data aggregation
- `bandwidth::Dates.Day`: Window size for rolling statistics
- `lag::Int64`: Number of time steps to lag the calculation

# Returns
- `EWSMetricSpecification`: Complete specification with auto-generated path

# Example
```julia
spec = EWSMetricSpecification(
    EWSMethod(Backward()),
    Dates.Day(7),
    Dates.Day(30),
    0
)
```
"""
function EWSMetricSpecification(
        method::EWSMethod,
        aggregation::Dates.Day,
        bandwidth::Dates.Day,
        lag::Int64,
    )
    aggregation_days_val = Dates.value(aggregation)
    bandwidth_days_val = Dates.value(bandwidth)

    return EWSMetricSpecification(
        method,
        aggregation,
        bandwidth,
        lag,
        _EWSMetricSpecification_path(
            method,
            aggregation_days_val,
            bandwidth_days_val,
            lag,
        ),
    )
end

"""
    EWSMetricSpecification(method, aggregation, bandwidth, lag)

Constructor for `EWSMetricSpecification` accepting general `DatePeriod` types.

This constructor provides flexibility by accepting any `DatePeriod` type for
aggregation and bandwidth parameters, automatically converting them to `Day`
periods for internal consistency.

# Arguments
- `method::EWSMethod`: The windowing method for EWS calculations
- `aggregation::Dates.DatePeriod`: Time period for data aggregation (converted to days)
- `bandwidth::Dates.DatePeriod`: Window size for rolling statistics (converted to days)
- `lag::Int64`: Number of time steps to lag the calculation

# Returns
- `EWSMetricSpecification`: Complete specification with converted periods and auto-generated path

# Example
```julia
spec = EWSMetricSpecification(
    EWSMethod(Backward()),
    Dates.Week(1),      # Becomes Day(7)
    Dates.Month(1),     # Becomes Day(30)
    0
)
```
"""
function EWSMetricSpecification(
        method::EWSMethod, aggregation::Dates.DatePeriod, bandwidth::Dates.DatePeriod, lag::Int64
    )
    aggregation_days_val = Dates.days(aggregation)
    bandwidth_days_val = Dates.days(bandwidth)
    aggregation_days = Dates.Day(aggregation_days_val)
    bandwidth_days = Dates.Day(bandwidth_days_val)

    return EWSMetricSpecification(
        method,
        aggregation_days,
        bandwidth_days,
        lag,
        _EWSMetricSpecification_path(
            method,
            aggregation_days_val,
            bandwidth_days_val,
            lag,
        ),
    )
end

"""
    _EWSMetricSpecification_path(method, aggregation, bandwidth, lag)

Internal function to generate standardized directory paths for EWS metric storage.

Creates a hierarchical directory structure that organizes EWS results by their
calculation parameters, enabling efficient caching and retrieval of computed metrics.

# Arguments
- `method::EWSMethod`: The windowing method used
- `aggregation::Int64`: Aggregation period in days
- `bandwidth::Int64`: Bandwidth window size in days
- `lag::Int64`: Lag value for the calculation

# Returns
- `String`: Hierarchical path string for file system organization

# Path Structure
The generated path follows the pattern:
```
ews-method_<method>/ews-aggregation-days_<days>/ews-bandwidth-days_<days>/ews-lag_<lag>
```

# Example
```julia
path = _EWSMetricSpecification_path(EWSMethod(Backward()), 7, 30, 0)
# Returns: "ews-method_backward/ews-aggregation-days_7/ews-bandwidth-days_30/ews-lag_0"
```
"""
function _EWSMetricSpecification_path(
        method::EWSMethod,
        aggregation::Int64,
        bandwidth::Int64,
        lag::Int64,
    )
    return joinpath(
        "ews-method_$(method_string(method))",
        "ews-aggregation-days_$(aggregation)",
        "ews-bandwidth-days_$(bandwidth)",
        "ews-lag_$(lag)",
    )
end


"""
    EWSMetrics

Container for computed Early Warning Signal metrics and their Kendall's tau trend statistics.

This struct stores the complete set of EWS metrics calculated from time series data,
including both the raw metric time series and their corresponding trend statistics
(Kendall's tau) that quantify the strength of temporal trends leading up to critical transitions.

# Fields

## Configuration
- `ews_specification::EWSMetricSpecification`: Parameters used to compute these metrics

## Time Series Metrics
- `mean::Vector{Float64}`: Rolling mean values over time
- `variance::Vector{Float64}`: Rolling variance values over time
- `coefficient_of_variation::Vector{Float64}`: Rolling CV (std/mean) values over time
- `index_of_dispersion::Vector{Float64}`: Rolling variance-to-mean ratio over time
- `skewness::Vector{Float64}`: Rolling skewness values over time
- `kurtosis::Vector{Float64}`: Rolling kurtosis values over time
- `autocovariance::Vector{Float64}`: Rolling lag-1 autocovariance over time
- `autocorrelation::Vector{Float64}`: Rolling lag-1 autocorrelation over time

## Trend Statistics (Kendall's Tau)
- `mean_tau::Float64`: Kendall's tau for the mean time series trend
- `variance_tau::Float64`: Kendall's tau for the variance time series trend
- `coefficient_of_variation_tau::Float64`: Kendall's tau for the CV time series trend
- `index_of_dispersion_tau::Float64`: Kendall's tau for the index of dispersion trend
- `skewness_tau::Float64`: Kendall's tau for the skewness time series trend
- `kurtosis_tau::Float64`: Kendall's tau for the kurtosis time series trend
- `autocovariance_tau::Float64`: Kendall's tau for the autocovariance time series trend
- `autocorrelation_tau::Float64`: Kendall's tau for the autocorrelation time series trend

# Usage
The tau values represent the strength and direction of trends in each metric, with values
closer to ±1 indicating stronger monotonic trends. Positive tau values suggest increasing
trends (potential early warning signals), while negative values suggest decreasing trends.

# Example
```julia
metrics = EWSMetrics(
    ews_specification = spec,
    mean = [1.0, 1.1, 1.2, 1.3],
    variance = [0.1, 0.15, 0.2, 0.3],
    # ... other time series
    mean_tau = 0.67,
    variance_tau = 0.89,
    # ... other tau values
)
```
"""
Base.@kwdef struct EWSMetrics
    ews_specification::EWSMetricSpecification
    mean::Vector{Float64}
    variance::Vector{Float64}
    coefficient_of_variation::Vector{Float64}
    index_of_dispersion::Vector{Float64}
    skewness::Vector{Float64}
    kurtosis::Vector{Float64}
    autocovariance::Vector{Float64}
    autocorrelation::Vector{Float64}
    mean_tau::Float64
    variance_tau::Float64
    coefficient_of_variation_tau::Float64
    index_of_dispersion_tau::Float64
    skewness_tau::Float64
    kurtosis_tau::Float64
    autocovariance_tau::Float64
    autocorrelation_tau::Float64
end

"""
    AbstractEWSThresholdWindowType

Abstract base type for Early Warning Signal threshold window strategies.

This type serves as the parent for different approaches to defining the baseline
period used for calculating EWS alert thresholds.
"""
Base.@kwdef struct EnsembleEWSMetrics
    emergent_ews_metrics::StructVector{EWSMetrics}
    null_ews_metrics::StructVector{EWSMetrics}
end

abstract type AbstractEWSThresholdWindowType end

"""
    ExpandingThresholdWindow

Threshold window strategy that expands over time.

This approach uses an expanding window for threshold calculation, where the baseline
period grows to include all available historical data up to each time point. This
provides more stable thresholds as more data becomes available but may be less
sensitive to recent changes in system dynamics.
"""
struct ExpandingThresholdWindow end

"""
    RollingThresholdWindow

Threshold window strategy that uses a fixed-size rolling window.

This approach uses a rolling window of fixed size for threshold calculation,
maintaining a constant-length baseline period that moves forward in time. This
provides more adaptive thresholds that can respond to recent changes in system
dynamics but may be less stable with limited data.
"""
struct RollingThresholdWindow end

"""
    EWSThresholdWindowType

Sum type representing the available threshold window strategies for EWS calculations.

This type determines how the baseline period is defined when calculating alert
thresholds for Early Warning Signals.

# Variants
- `ExpandingThresholdWindow`: Uses expanding windows (more stable, less adaptive)
- `RollingThresholdWindow`: Uses rolling windows (more adaptive, potentially less stable)

# Example
```julia
window_type = EWSThresholdWindowType(ExpandingThresholdWindow())
# or
window_type = EWSThresholdWindowType(RollingThresholdWindow())
```
"""
LightSumTypes.@sumtype EWSThresholdWindowType(
    ExpandingThresholdWindow,
    RollingThresholdWindow
) <: AbstractEWSThresholdWindowType

"""
    AbstractEWSEndDateType

Abstract base type for Early Warning Signal calculation end date strategies.

This type serves as the parent for different approaches to determining when to stop
EWS calculations relative to outbreak or epidemic phases.
"""
abstract type AbstractEWSEndDateType end

"""
    ReffStart

End date strategy that stops EWS calculations at the start of the R_eff > 1 period.

This strategy terminates EWS calculations when the effective reproduction number
first exceeds 1, marking the beginning of exponential growth phase.
"""
struct ReffStart end

"""
    ReffEnd

End date strategy that stops EWS calculations at the end of the R_eff > 1 period.

This strategy terminates EWS calculations when the effective reproduction number
returns to or below 1, marking the end of the exponential growth phase.
"""
struct ReffEnd end

"""
    OutbreakStart

End date strategy that stops EWS calculations at the start of the outbreak period.

This strategy terminates EWS calculations at the beginning of the defined outbreak
period, which may be based on case count thresholds or other outbreak criteria.
"""
struct OutbreakStart end

"""
    OutbreakMiddle

End date strategy that stops EWS calculations at the middle of the outbreak period.

This strategy terminates EWS calculations at the midpoint of the defined outbreak
period, providing a balance between early detection and outbreak progression.
"""
struct OutbreakMiddle end

"""
    OutbreakEnd

End date strategy that stops EWS calculations at the end of the outbreak period.

This strategy terminates EWS calculations at the conclusion of the defined outbreak
period, allowing for the full outbreak trajectory to be analyzed.
"""
struct OutbreakEnd end

"""
    EWSEndDateType

Sum type representing the available end date strategies for EWS calculations.

This type determines when to stop calculating Early Warning Signals relative to
outbreak dynamics, affecting the temporal scope of the analysis and the lead time
available for early detection.

# Variants
- `ReffStart`: Stop at the beginning of R_eff > 1 period
- `ReffEnd`: Stop at the end of R_eff > 1 period
- `OutbreakStart`: Stop at the beginning of outbreak period
- `OutbreakMiddle`: Stop at the middle of outbreak period
- `OutbreakEnd`: Stop at the end of outbreak period

# Example
```julia
end_date_type = EWSEndDateType(OutbreakStart())
# or
end_date_type = EWSEndDateType(ReffStart())
```
"""
LightSumTypes.@sumtype EWSEndDateType(
    ReffStart,
    ReffEnd,
    OutbreakStart,
    OutbreakMiddle,
    OutbreakEnd
) <: AbstractEWSEndDateType
