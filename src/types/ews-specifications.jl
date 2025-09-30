using Dates: Dates
using LightSumTypes: @sumtype

export EWSMethod,
    Backward,
    Centered,
    EWSMetricSpecification,
    EWSMetrics,
    EWSThresholdWindowType,
    ExpandingThresholdWindow,
    RollingThresholdWindow,
    EWSEndDateType,
    Reff_start,
    Reff_end,
    Outbreak_start,
    Outbreak_middle,
    Outbreak_end

abstract type AbstractEWSMethod end
struct Backward end
struct Centered end
@sumtype EWSMethod(Backward, Centered) <: AbstractEWSMethod

struct EWSMetricSpecification
    method::EWSMethod
    aggregation::Dates.Day
    bandwidth::Dates.Day
    lag::Int64
    dirpath::String
end

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

function get_ews_metric_specification_description(ews_metric_specification)
    return "Method: $(method_string(ews_metric_specification.method)), Aggregation: $(ews_metric_specification.aggregation), Bandwidth: $(ews_metric_specification.bandwidth), Lag: $(ews_metric_specification.lag)"
end

method_string(method::EWSMethod) = lowercase(split(string(method), "::")[1])

struct EWSMetrics
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

abstract type AbstractEWSThresholdWindowType end
struct ExpandingThresholdWindow end
struct RollingThresholdWindow end
@sumtype EWSThresholdWindowType(ExpandingThresholdWindow, RollingThresholdWindow) <: AbstractEWSThresholdWindowType

abstract type AbstractEWSEndDateType end

struct Reff_start end
struct Reff_end end
struct Outbreak_start end
struct Outbreak_middle end
struct Outbreak_end end

@sumtype EWSEndDateType(Reff_start, Reff_end, Outbreak_start, Outbreak_middle, Outbreak_end) <: AbstractEWSEndDateType
