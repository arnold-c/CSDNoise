export get_enddate_thresholds,
    exceeds_ews_threshold,
    get_ews_metric_vec

get_enddate_thresholds(data_arrs, enddate_type::EWSEndDateType) = get_enddate_thresholds(data_arrs, LightSumTypes.variant(enddate_type))

function get_enddate_thresholds(data_arrs, enddate_type::Union{Reff_start, Reff_end})
    return data_arrs.ensemble_single_Reff_thresholds_vec
end

function get_enddate_thresholds(data_arrs, enddate_type::Union{Outbreak_start, Outbreak_middle, Outbreak_end})
    return data_arrs.emergent_outbreak_threshold_vecs
end

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

# TODO: Implement this method. Place holder as needed for type stability
function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::RollingThresholdWindow,
        quantile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    return fill(false, 2)
end

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


function _expanding_ews_thresholds(
        ews_vec::Vector{Float64},
        quantile::Float64,
        burn_in_index::Int64,
    )
    # ) where {T1 <: Integer, F <: AbstractFloat}
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
