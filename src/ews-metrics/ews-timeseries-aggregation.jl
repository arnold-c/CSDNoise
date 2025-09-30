using Dates: Dates
using StatsBase: StatsBase

export aggregate_bandwidth,
    aggregate_thresholds_vec,
    aggregate_Reff_vec,
    aggregate_timeseries

"""
    aggregate_bandwidth(ews_spec::EWSMetricSpecification)

Calculate the effective bandwidth after temporal aggregation.

# Arguments
- `ews_spec::EWSMetricSpecification`: EWS specification containing bandwidth and aggregation parameters

# Returns
- `Int`: Aggregated bandwidth (original bandwidth divided by aggregation period)

# Throws
- `AssertionError`: If the aggregated bandwidth is not a positive integer
"""
function aggregate_bandwidth(ews_spec::EWSMetricSpecification)
    aggregate_bandwidth = ews_spec.bandwidth ÷ ews_spec.aggregation
    @assert aggregate_bandwidth > 0
    @assert isinteger(aggregate_bandwidth)

    return aggregate_bandwidth
end

"""
    aggregate_thresholds_vec(thresholdsvec, aggregation)

Aggregate a vector of threshold indicators over time periods.

Uses a binary aggregation function where any threshold crossing within the aggregation
period results in a positive indicator for that period.

# Arguments
- `thresholdsvec`: Vector of threshold indicators (typically binary)
- `aggregation`: Temporal aggregation period

# Returns
- Aggregated threshold vector where each element indicates if any threshold was crossed in that period
"""
function aggregate_thresholds_vec(thresholdsvec, aggregation)
    return aggregate_timeseries(thresholdsvec, aggregation, x -> sum(x) >= 1)
end

"""
    aggregate_Reff_vec(Reff_vec, aggregation)

Aggregate effective reproduction number (Reff) values over time periods using mean aggregation.

# Arguments
- `Reff_vec`: Vector of effective reproduction numbers
- `aggregation`: Temporal aggregation period

# Returns
- Aggregated Reff vector with mean values for each aggregation period
"""
function aggregate_Reff_vec(Reff_vec, aggregation)
    return aggregate_timeseries(Reff_vec, aggregation, StatsBase.mean)
end


"""
    aggregate_timeseries(timeseries, aggregation::DatePeriod, stat_function=sum)

Aggregate time series data over specified date periods using a statistical function.

# Arguments
- `timeseries`: Input time series data
- `aggregation::DatePeriod`: Temporal aggregation period (e.g., Day(7) for weekly)
- `stat_function`: Statistical function to apply during aggregation (default: sum)

# Returns
- Aggregated time series with reduced temporal resolution

# Notes
- Returns original time series unchanged if aggregation is Day(1)
- Uses sum as default aggregation function, suitable for count data
"""
function aggregate_timeseries(
        timeseries,
        aggregation::T1,
        stat_function = sum,
    ) where {T1 <: Dates.DatePeriod}
    if aggregation == Dates.Day(1)
        return timeseries
    end
    return _aggregate_timeseries(timeseries, aggregation, stat_function)
end

"""
    _aggregate_timeseries(timeseries, aggregation::DatePeriod, stat_function=mean)

Internal implementation for aggregating time series data over date periods.

This function performs the actual aggregation work by dividing the time series into
chunks based on the aggregation period and applying the statistical function to each chunk.

# Arguments
- `timeseries`: Input time series data
- `aggregation::DatePeriod`: Temporal aggregation period
- `stat_function`: Statistical function to apply (default: mean, unlike public version which defaults to sum)

# Returns
- Aggregated time series with length `length(timeseries) ÷ aggregation_days`

# Notes
- Uses mean as default (different from public `aggregate_timeseries` which uses sum)
- Assumes aggregation period can be evenly divided into the time series length. If not, the end of the time series
will be truncated
"""
function _aggregate_timeseries(
        timeseries,
        aggregation::T1,
        stat_function = StatsBase.mean,
    ) where {T1 <: Dates.DatePeriod}
    aggregation_days = Dates.days(aggregation)
    aggregate_timeseries = zeros(
        eltype(stat_function(@view(timeseries[1:2]))),
        length(timeseries) ÷ aggregation_days,
    )
    for i in eachindex(aggregate_timeseries)
        aggregate_timeseries[i] = stat_function(
            @view(
                timeseries[((i - 1) * aggregation_days + 1):(i * aggregation_days)]
            )
        )
    end
    return aggregate_timeseries
end
