function survival_plot_lines!(
        gl::T1,
        times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        linestyle = :solid,
        emergent_color = "#C2192D",
        null_color = "#453948",
        alpha = 1.0,
        nsims = 100,
        facet_title = "Survival",
        trim_burnin = true,
        burnin = Dates.Year(5),
        facet_fontsize = 20,
        xlabelsize = 22,
        ylabelsize = 22,
    ) where {T1 <: GLMakie.GridLayout}
    surv_ax = Axis(
        gl[2, 1];
        xlabel = "Time (Years)",
        ylabel = "Survival Numbers",
        xticks = 1:1:10,
        limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
    )
    Box(gl[1, 1]; color = :lightgray, strokevisible = false)
    Label(
        gl[1, 1],
        facet_title;
        fontsize = facet_fontsize,
        padding = (0, 0, 0, 0),
        valign = :bottom,
        tellwidth = false,
    )
    rowsize!(gl, 2, Relative(0.9))

    survival_plot_lines!(
        surv_ax,
        times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        linestyle = linestyle,
        emergent_color = emergent_color,
        null_color = null_color,
        alpha = alpha,
        nsims = nsims,
        trim_burnin = trim_burnin,
        burnin = burnin,
    )
    return nothing
end

function survival_plot_lines!(
        surv_ax::T1,
        times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        linestyle = :solid,
        emergent_color = "#C2192D",
        null_color = "#453948",
        alpha = 1.0,
        nsims = 100,
        trim_burnin = true,
        burnin = Dates.Year(5),
    ) where {T1 <: GLMakie.Axis}
    lines!(
        surv_ax,
        detection_survival_times,
        detection_survival_vec;
        color = (emergent_color, alpha),
        linestyle = linestyle,
        label = "Detection",
    )

    lines!(
        surv_ax,
        null_survival_times,
        null_survival_vec;
        color = (null_color, alpha),
        linestyle = linestyle,
        label = "Null",
    )

    if trim_burnin
        xlims!(surv_ax, Dates.value(burnin), times[end])
    else
        vlines!(
            surv_ax, Int64(Dates.value(burnin)); color = :black,
            linestyle = :dash,
        )
    end
    return nothing
end
