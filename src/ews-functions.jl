using SumTypes
using StatsBase: StatsBase

@sum_type EWSThresholdWindow begin
    Expanding
    Rolling
end

function expanding_ews_thresholds(
    ewsmetrics::T1,
    metric::T2,
    window_type::T3;
    percentiles = (0.8, 0.95),
    burn_in = 10,
) where {T1<:EWSMetrics,T2<:Symbol,T3<:EWSThresholdWindow}
    ews_vec = getproperty(ewsmetrics, metric)

    ews_distributions = fill(NaN, length(ews_vec), length(percentiles))
    exceeds_thresholds = Array{Bool}(
        undef, length(ews_vec), length(percentiles)
    )

    for i in eachindex(ews_vec)
        if i <= burn_in
            continue
        end
        ews_distributions[i, :] .= map(
            p -> StatsBase.quantile(@view(ews_vec[1:i]), p),
            percentiles,
        )

        exceeds_thresholds[i, :] .= ews_vec[i] .>= ews_distributions[i, :]
    end

    return ews_distributions, exceeds_thresholds
end
