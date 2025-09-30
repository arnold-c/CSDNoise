using StatsBase: StatsBase

export calculate_auc,
    get_tau

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

function get_tau(
        ews_metrics;
        tau_metric = :variance_tau,
        statistic_function = StatsBase.mean,
    )
    tau_vector = getproperty(ews_metrics, tau_metric)

    return statistic_function(tau_vector)
end
