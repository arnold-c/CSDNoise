using StatsBase: StatsBase

export simulation_tau_heatmap_df!

function simulation_tau_heatmap_df!(
        ews_df,
        test_ewsmetrics,
        ews_metric;
        individual_test_specification = IndividualTestSpecification(1.0, 1.0, 0),
        ews_metric_specification = EWSMetricSpecification(Centered, 7, 52, 1),
        ews_enddate_type = Reff_start::EWSEndDateType,
        statistic_function = StatsBase.mean,
    )
    @assert names(ews_df) == [
        "ews_metric",
        "test_specification",
        "ews_metric_specification",
        "ews_enddate_type",
        "ews_metric_value",
        "ews_metric_vector",
    ]

    ews_metric_tau_sym = Symbol(ews_metric, "_tau")
    ews_tau = get_tau(
        test_ewsmetrics;
        tau_metric = ews_metric_tau_sym,
        statistic_function = statistic_function,
    )

    push!(
        ews_df,
        (
            ews_metric,
            individual_test_specification,
            ews_metric_specification,
            ews_enddate_type,
            ews_tau,
            getproperty(test_ewsmetrics, ews_metric_tau_sym),
        ),
    )

    return nothing
end
