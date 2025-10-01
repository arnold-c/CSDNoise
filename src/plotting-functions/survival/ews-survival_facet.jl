function ews_survival_facet!(
        gl,
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        facet_title = "Survival",
        ews_aggregation = Dates.Day(7),
        burnin = Dates.Year(5),
        endpoint_aggregation = Dates.Day(30),
        emergent_color = "#C2192D",
        null_color = "#453948",
        alpha = 1.0,
        trim_burnin = true,
    )
    times,
        enddate_times,
        enddate_counts,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec,
        nsims = prepare_survival_facet_params(
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        ews_aggregation = ews_aggregation,
        endpoint_aggregation = endpoint_aggregation,
    )

    survival_plot_histogram!(
        gl,
        enddate_times,
        enddate_counts;
        trim_burnin = trim_burnin,
        burnin = burnin,
    )

    survival_plot_lines!(
        gl,
        enddate_times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        emergent_color = emergent_color,
        null_color = null_color,
        alpha = alpha,
        nsims = nsims,
        facet_title = facet_title,
        trim_burnin = trim_burnin,
        burnin = burnin,
    )

    return nothing
end
