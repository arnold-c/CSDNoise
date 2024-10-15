using GLMakie
using ColorSchemes
using UnPack
using DataFrames
using FLoops
using LaTeXStrings
using NaNMath: NaNMath
using Dates
using Debugger: Debugger
using StructArrays: StructArrays

function Reff_plot(
    incidencearr,
    Reffarr,
    Reff_thresholds_vec,
    thresholdsarr,
    timeparams;
    sim = 1,
    outbreak_colormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    threshold = 5,
    kwargs...,
)
    @unpack trange = timeparams
    times = collect(trange) ./ 365
    kwargs_dict = Dict(kwargs)

    period_sum_arr = zeros(Int64, length(times), 2)
    for (lower, upper, periodsum, outbreakstatus) in
        eachrow(thresholdsarr[sim])
        period_sum_arr[lower:upper, 1] .= periodsum
        period_sum_arr[lower:upper, 2] .= outbreakstatus
    end

    Reff_dur_arr = zeros(Int64, length(times), 2)
    for (lower, upper, duration) in
        eachrow(Reff_thresholds_vec[sim])
        Reff_dur_arr[lower:upper, 1] .= duration
        Reff_dur_arr[lower:upper, 2] .= 1
    end

    fig = Figure()
    incax = Axis(
        fig[1, 1]; ylabel = "Incidence"
    )
    reffax = Axis(
        fig[2, 1]; ylabel = "Reff"
    )
    reff_sum_ax = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = "Reff Duration"
    )

    linkxaxes!(incax, reffax, reff_sum_ax)

    lines!(
        incax,
        times,
        incidencearr[:, 1, sim];
        color = period_sum_arr[:, 2],
        colormap = outbreak_colormap,
    )
    hlines!(
        incax,
        threshold;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    )
    lines!(
        reffax,
        times,
        Reffarr[:, sim];
        color = Reff_dur_arr[:, 2],
        colormap = Reff_colormap,
        linewidth = 3,
    )
    hlines!(
        reffax,
        1.0;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    )
    lines!(
        reff_sum_ax,
        times,
        Reff_dur_arr[:, 1];
        color = Reff_dur_arr[:, 2],
        colormap = Reff_colormap,
        linewidth = 3,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [incax, reffax, reff_sum_ax],
        )
    end

    if haskey(kwargs_dict, :ylims_Reffdur)
        ylims!(reff_sum_ax, kwargs_dict[:ylims_Reffdur])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    map(hidexdecorations!, [incax, reffax])

    Legend(
        fig[1, 2],
        [PolyElement(; color = col) for col in outbreak_colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    Legend(
        fig[2:3, 2],
        [PolyElement(; color = col) for col in Reff_colormap],
        ["Reff < 1", "Reff >= 1"],
        "Reff Status",
    )
    return fig
end

function incidence_prevalence_plot(
    incidencearr,
    ensemblearr,
    thresholdsarr,
    timeparams;
    sim = 1,
    colormap = [BASE_COLOR, OUTBREAK_COLOR],
    threshold = 5,
    kwargs...,
)
    @unpack trange = timeparams
    times = collect(trange) ./ 365
    kwargs_dict = Dict(kwargs)

    period_sum_arr = zeros(Int64, length(times), 2)
    for (lower, upper, periodsum, outbreakstatus) in
        eachrow(thresholdsarr[sim])
        period_sum_arr[lower:upper, 1] .= periodsum
        period_sum_arr[lower:upper, 2] .= outbreakstatus
    end

    fig = Figure()
    ax_prev = Axis(fig[1, 1]; ylabel = "Prevalence")
    ax_inc = Axis(fig[2, 1]; ylabel = "Incidence")
    ax_periodsum = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = "Period Sum"
    )

    linkxaxes!(ax_prev, ax_inc, ax_periodsum)

    lines!(
        ax_prev,
        times,
        ensemblearr[:, 2, sim];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )
    lines!(
        ax_inc,
        times,
        incidencearr[:, 1, sim];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )
    hlines!(
        ax_inc,
        threshold;
        color = :black,
        linestyle = :dash,
        linewidth = 2,
    )
    barplot!(
        ax_periodsum,
        times,
        period_sum_arr[:, 1];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )

    map(hidexdecorations!, [ax_prev, ax_inc])

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [ax_prev, ax_inc, ax_periodsum],
        )
    end

    if haskey(kwargs_dict, :ylims_periodsum)
        ylims!(ax_periodsum, kwargs_dict[:ylims_periodsum])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(ax_inc, kwargs_dict[:ylims_inc])
    end

    Legend(
        fig[:, 2],
        [PolyElement(; color = col) for col in colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    return fig
end

function incidence_testing_plot(
    incvec,
    outbreak_status_vec,
    noisevec,
    testvec,
    test_movingavg_vec,
    timeparams;
    aggregation = 1,
    outbreakcolormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    plottitle = "",
    kwargs...,
)
    times = collect(timeparams.trange) ./ 365
    aggregated_vec_length = length(incvec)
    times = times[1:aggregation:(aggregation * aggregated_vec_length)]

    kwargs_dict = Dict(kwargs)

    inc_test_fig = Figure()
    inc_test_ax1 = Axis(inc_test_fig[1, 1]; ylabel = "Incidence")
    inc_test_ax2 = Axis(inc_test_fig[2, 1]; ylabel = "Incidence + Noise")
    inc_test_ax3 = Axis(inc_test_fig[3, 1]; ylabel = "Test Positive")
    inc_test_ax4 = Axis(
        inc_test_fig[4, 1];
        xlabel = "Time (years)",
        ylabel = "7d Avg Test Positive",
    )

    lines!(
        inc_test_ax1, times, incvec;
        color = outbreak_status_vec,
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax2, times, incvec .+ noisevec;
        color = outbreak_status_vec,
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax3, times, testvec
    )
    lines!(
        inc_test_ax4, times, test_movingavg_vec
    )

    linkxaxes!(inc_test_ax1, inc_test_ax2, inc_test_ax3, inc_test_ax4)

    map(hidexdecorations!, [inc_test_ax1, inc_test_ax3])

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [inc_test_ax1, inc_test_ax2, inc_test_ax3, inc_test_ax4],
        )
    end

    if haskey(kwargs_dict, :ylims)
        map(
            ax -> ylims!(ax, kwargs_dict[:ylims]),
            [inc_test_ax1, inc_test_ax2, inc_test_ax3, inc_test_ax4],
        )
    end

    map(
        ax -> hlines!(
            ax,
            5 * aggregation;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax1, inc_test_ax2],
    )

    Label(
        inc_test_fig[0, :, Top()],
        plottitle,
    )

    Legend(
        inc_test_fig[1:2, 2],
        [PolyElement(; color = col) for col in outbreakcolormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    rowsize!(inc_test_fig.layout, 0, 5)
    colsize!(inc_test_fig.layout, 1, Relative(0.92))

    return inc_test_fig
end

function ensemble_incarr_Reff_plot(
    incarr,
    Reffarr,
    Reff_thresholds_vec,
    thresholdsarr;
    aggregation = 1,
    linewidth = 3,
    hlinewidth = 2,
    threshold = 5 * aggregation,
    outbreak_alpha = 0.2,
    Reff_alpha = 0.2,
    outbreak_colormap = [
        (BASE_COLOR, outbreak_alpha),
        (OUTBREAK_COLOR, outbreak_alpha),
    ],
    Reff_colormap = [
        (BASE_COLOR, Reff_alpha),
        (REFF_GT_ONE_COLOR, Reff_alpha),
    ])
    fig = Figure()
    times = collect(1:aggregation:size(incarr, 1)) ./ 365

    fig = Figure()
    reffax = Axis(
        fig[1, 1]; ylabel = "Reff"
    )
    incax = Axis(
        fig[2, 1]; ylabel = "Incidence"
    )

    linkxaxes!(incax, reffax)
    aggregated_vec_length = size(Reffarr, 1) ÷ aggregation * aggregation
    times = times[1:aggregation:aggregated_vec_length]

    for sim in axes(incarr, 3)
        Reff_above_one = zeros(Int64, length(times))

        for (lower, upper, _) in
            eachrow(Reff_thresholds_vec[sim])
            Reff_above_one[lower:upper] .= 1
        end
        Reff_above_one = aggregate_timeseries(Reff_above_one, aggregation)

        filtered_thresholdsarr = thresholdsarr[sim][
            (thresholdsarr[sim][:, 4] .== 1), :,
        ]

        outbreak_status_vec = zeros(Int64, length(times))
        for (lower, upper, periodsum, outbreakstatus) in
            eachrow(filtered_thresholdsarr)
            outbreak_status_vec[(lower ÷ aggregation):1:(upper ÷ aggregation)] .=
                outbreakstatus
        end

        lines!(
            reffax,
            times,
            @view(
                Reffarr[
                    1:aggregation:aggregated_vec_length,
                    sim,
                ]
            );
            linewidth = linewidth,
            color = Reff_above_one,
            colormap = Reff_colormap,
        )

        lines!(
            incax,
            times,
            aggregate_timeseries(@view(incarr[:, 1, sim]), aggregation);
            linewidth = linewidth,
            color = outbreak_status_vec,
            colormap = outbreak_colormap,
        )
    end

    hlines!(reffax, [0.9, 1]; linestyle = :dash)
    hlines!(incax, threshold; linestyle = :dash)

    Legend(
        fig[1, 2],
        [PolyElement(; color = (col[1], 1)) for col in Reff_colormap],
        ["Reff < 1", "Reff >= 1"],
        "Reff Status",
    )

    Legend(
        fig[2, 2],
        [PolyElement(; color = (col[1], 1)) for col in outbreak_colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    return fig
end
