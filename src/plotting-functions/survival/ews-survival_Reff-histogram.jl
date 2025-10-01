function survival_plot_histogram!(
        gl,
        enddate_times,
        enddate_counts;
        tipping_point_color = ("#A4BED5", 0.7),
        trim_burnin = true,
        burnin = Dates.Year(5),
    )
    hist_ax = Axis(
        gl[1, 1];
        limits = (0, maximum(enddate_times), nothing, nothing),
    )

    hidexdecorations!(hist_ax)
    hideydecorations!(hist_ax)

    barplot!(
        hist_ax,
        enddate_times,
        enddate_counts;
        color = tipping_point_color,
    )

    if trim_burnin
        xlims!(hist_ax, Dates.value(burnin), enddate_times[end])
    end

    return nothing
end
