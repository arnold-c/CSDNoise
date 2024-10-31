using GLMakie
using UnPack: @unpack
using NaNMath: NaNMath
using StructArrays: StructArray
using DataFrames

function Reff_ews_plot(
    incvec,
    Reffvec,
    Reff_thresholds,
    ewsmetrics,
    null_ewsmetrics,
    ewsmetric::S,
    thresholdsarr,
    exceeds_thresholds_vec,
    detection_index,
    null_exceeds_thresholds_vec,
    null_detection_index,
    timeparams;
    outbreak_colormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    threshold = 5,
    metric_color = Makie.wong_colors()[1],
    kwargs...,
) where {S<:Symbol}
    @unpack trange = timeparams
    @unpack aggregation, bandwidth, lag, method = ewsmetrics.ews_specification
    method = split(string(method), "::")[1]

    times = collect(trange) ./ 365
    aggregated_vec_length = length(Reffvec)
    times = times[1:aggregation:(aggregation * aggregated_vec_length)]

    kwargs_dict = Dict(kwargs)

    outbreak_status_vec = zeros(Int64, length(times))
    for (lower, upper) in
        eachrow(thresholdsarr)
        outbreak_status_vec[lower:1:upper] .= 1
    end

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds)
        Reff_above_one[lower:1:upper] .= 1
    end

    ewsmetric_vec = getproperty(ewsmetrics, ewsmetric)
    ewsmetric_endpoint = length(ewsmetric_vec)
    ewsmetric_vec = vcat(
        ewsmetric_vec,
        fill(NaN, length(times) - ewsmetric_endpoint),
    )
    ewsmetric_tau = getproperty(ewsmetrics, Symbol(String(ewsmetric) * "_tau"))

    null_ewsmetric_vec = getproperty(null_ewsmetrics, ewsmetric)
    null_ewsmetric_endpoint = length(null_ewsmetric_vec)
    null_ewsmetric_vec = vcat(
        null_ewsmetric_vec,
        fill(NaN, length(times) - null_ewsmetric_endpoint),
    )

    null_ewsmetric_tau = getproperty(
        null_ewsmetrics, Symbol(String(ewsmetric) * "_tau")
    )

    fig = Figure()
    reffax = Axis(
        fig[1, 1]; ylabel = "Reff"
    )
    incax = Axis(
        fig[2, 1]; ylabel = "Incidence"
    )
    metric_ax = Axis(
        fig[3, 1]; xlabel = "Time (years)",
        ylabel = "Detection " * String(ewsmetric),
    )
    null_metric_ax = Axis(
        fig[4, 1]; xlabel = "Time (years)", ylabel = "Null " * String(ewsmetric)
    )

    linkxaxes!(incax, reffax, metric_ax, null_metric_ax)

    map(
        ax ->
            vspan!(
                ax, 0, (ewsmetric_endpoint - 1) * aggregation / 365;
                color = (:lightgrey, 0.5),
            ),
        [incax, reffax, metric_ax, null_metric_ax],
    )

    line_and_hline!(
        reffax,
        times,
        Reffvec,
        1;
        color = Reff_above_one,
        colormap = Reff_colormap,
    )

    outbreak_threshold_hline = aggregation == 1 ? threshold : NaN
    line_and_hline!(
        incax,
        times,
        incvec,
        outbreak_threshold_hline;
        color = outbreak_status_vec,
        colormap = outbreak_colormap,
    )

    lines!(
        metric_ax,
        times,
        ewsmetric_vec;
        color = metric_color,
        linewidth = 3,
    )
    if !isnothing(detection_index)
        exceeds_threshold_indices = findall(
            x -> x == true, exceeds_thresholds_vec
        )

        vlines!(metric_ax, times[detection_index]; color = :black)
        scatter!(
            metric_ax,
            times[exceeds_threshold_indices],
            ewsmetric_vec[exceeds_threshold_indices];
            markersize = 10,
            strokecolor = :grey20,
            strokewidth = 2,
            color = (:grey20, 0.4),
        )
    end

    lines!(
        null_metric_ax,
        times,
        null_ewsmetric_vec;
        color = metric_color,
        linewidth = 3,
    )
    if !isnothing(null_detection_index)
        null_exceeds_threshold_indices = findall(
            x -> x == true, exceeds_thresholds_vec
        )

        vlines!(metric_ax, times[null_detection_index]; color = :black)
        scatter!(
            metric_ax,
            times[null_exceeds_threshold_indices],
            ewsmetric_vec[null_exceeds_threshold_indices];
            markersize = 10,
            strokecolor = :grey20,
            strokewidth = 2,
            color = (:grey20, 0.4),
        )
    end

    ewsmetric_extrema = extrema(
        replace(ewsmetric_vec[1:ewsmetric_endpoint], NaN => 0)
    )
    ewsmetric_range_buffer =
        abs(ewsmetric_extrema[2] - ewsmetric_extrema[1]) / 10

    ewsmetric_tau_yvalue =
        if ewsmetric_vec[ewsmetric_endpoint] >=
            ewsmetric_extrema[2] - ewsmetric_range_buffer
            ewsmetric_extrema[2] - ewsmetric_range_buffer
        elseif ewsmetric_vec[ewsmetric_endpoint] <=
            ewsmetric_extrema[1] + ewsmetric_range_buffer
            ewsmetric_extrema[1] + ewsmetric_range_buffer
        else
            ewsmetric_vec[ewsmetric_endpoint]
        end

    text!(
        metric_ax,
        times[ewsmetric_endpoint] + 0.5,
        ewsmetric_tau_yvalue;
        text = "τ = $(round(ewsmetric_tau; digits = 2))",
        justification = :right,
    )

    null_ewsmetric_extrema = extrema(
        replace(null_ewsmetric_vec[1:null_ewsmetric_endpoint], NaN => 0)
    )
    ewsmetric_range_buffer =
        abs(null_ewsmetric_extrema[2] - null_ewsmetric_extrema[1]) / 10

    null_ewsmetric_tau_yvalue =
        if null_ewsmetric_vec[null_ewsmetric_endpoint] >=
            null_ewsmetric_extrema[2] - ewsmetric_range_buffer
            null_ewsmetric_extrema[2] - ewsmetric_range_buffer
        elseif null_ewsmetric_vec[null_ewsmetric_endpoint] <=
            null_ewsmetric_extrema[1] + ewsmetric_range_buffer
            null_ewsmetric_extrema[1] + ewsmetric_range_buffer
        else
            null_ewsmetric_vec[null_ewsmetric_endpoint]
        end

    text!(
        null_metric_ax,
        times[null_ewsmetric_endpoint] + 0.5,
        null_ewsmetric_tau_yvalue;
        text = "τ = $(round(null_ewsmetric_tau; digits = 2))",
        justification = :right,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [incax, reffax, metric_ax, null_metric_ax],
        )
    end

    if haskey(kwargs_dict, :ylims_metric)
        ylims!(metric_ax, kwargs_dict[:ylims_metric])
        ylims!(null_metric_ax, kwargs_dict[:ylims_metric])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    map(hidexdecorations!, [incax, reffax])

    if haskey(kwargs_dict, :plottitle)
        Label(
            fig[0, :, Top()],
            kwargs_dict[:plottitle],
        )
    end
    Legend(
        fig[1, 2],
        [PolyElement(; color = col) for col in Reff_colormap],
        ["Reff < 1", "Reff >= 1"],
        "Reff Status",
    )

    Legend(
        fig[2, 2],
        [PolyElement(; color = col) for col in outbreak_colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    Label(
        fig[0, :],
        "EWS Metric: $(string(ewsmetric)), $(method), Lag: $(lag), Bandwidth: $(bandwidth), Aggregation: $(aggregation) days\nEWS calculated for shaded region",
    )

    rowsize!(fig.layout, 0, Relative(0.05))

    return fig
end

function Reff_ews_plot(
    incvec,
    Reffvec,
    Reff_thresholds,
    ewsmetrics,
    ewsmetric::S,
    thresholdsarr,
    exceeds_thresholds_vec,
    detection_index,
    timeparams;
    outbreak_colormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    threshold = 5,
    metric_color = Makie.wong_colors()[1],
    kwargs...,
) where {S<:Symbol}
    @unpack trange = timeparams
    @unpack aggregation, bandwidth, lag, method = ewsmetrics.ews_specification
    method = split(string(method), "::")[1]

    times = collect(trange) ./ 365
    aggregated_vec_length = length(Reffvec)
    times = times[1:aggregation:(aggregation * aggregated_vec_length)]

    kwargs_dict = Dict(kwargs)

    outbreak_status_vec = zeros(Int64, length(times))
    for (lower, upper) in
        eachrow(thresholdsarr)
        outbreak_status_vec[lower:1:upper] .= 1
    end

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds)
        Reff_above_one[lower:1:upper] .= 1
    end

    ewsmetric_vec = getproperty(ewsmetrics, ewsmetric)
    ewsmetric_endpoint = length(ewsmetric_vec)
    ewsmetric_vec = vcat(
        ewsmetric_vec,
        fill(NaN, length(times) - ewsmetric_endpoint),
    )

    ewsmetric_tau = getproperty(ewsmetrics, Symbol(String(ewsmetric) * "_tau"))

    fig = Figure()
    reffax = Axis(
        fig[1, 1]; ylabel = "Reff"
    )
    incax = Axis(
        fig[2, 1]; ylabel = "Incidence"
    )
    metric_ax = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = String(ewsmetric)
    )

    linkxaxes!(incax, reffax, metric_ax)

    map(
        ax ->
            vspan!(
                ax, 0, (ewsmetric_endpoint - 1) * aggregation / 365;
                color = (:lightgrey, 0.5),
            ),
        [incax, reffax, metric_ax],
    )

    line_and_hline!(
        reffax,
        times,
        Reffvec,
        1;
        color = Reff_above_one,
        colormap = Reff_colormap,
    )

    outbreak_threshold_hline = aggregation == 1 ? threshold : NaN
    line_and_hline!(
        incax,
        times,
        incvec,
        outbreak_threshold_hline;
        color = outbreak_status_vec,
        colormap = outbreak_colormap,
    )

    lines!(
        metric_ax,
        times,
        ewsmetric_vec;
        color = metric_color,
        linewidth = 3,
    )
    if !isnothing(detection_index)
        exceeds_threshold_indices = findall(
            x -> x == true, exceeds_thresholds_vec
        )

        vlines!(metric_ax, times[detection_index]; color = :black)
        scatter!(
            metric_ax,
            times[exceeds_threshold_indices],
            ewsmetric_vec[exceeds_threshold_indices];
            markersize = 10,
            strokecolor = :grey20,
            strokewidth = 2,
            color = (:grey20, 0.4),
        )
    end

    ewsmetric_extrema = extrema(
        replace(ewsmetric_vec[1:ewsmetric_endpoint], NaN => 0)
    )
    ewsmetric_range_buffer =
        abs(ewsmetric_extrema[2] - ewsmetric_extrema[1]) / 10

    ewsmetric_tau_yvalue =
        if ewsmetric_vec[ewsmetric_endpoint] >=
            ewsmetric_extrema[2] - ewsmetric_range_buffer
            ewsmetric_extrema[2] - ewsmetric_range_buffer
        elseif ewsmetric_vec[ewsmetric_endpoint] <=
            ewsmetric_extrema[1] + ewsmetric_range_buffer
            ewsmetric_extrema[1] + ewsmetric_range_buffer
        else
            ewsmetric_vec[ewsmetric_endpoint]
        end

    text!(
        metric_ax,
        times[ewsmetric_endpoint] + 0.5,
        ewsmetric_tau_yvalue;
        text = "τ = $(round(ewsmetric_tau; digits = 2))",
        justification = :right,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [incax, reffax, metric_ax],
        )
    end

    if haskey(kwargs_dict, :ylims_metric)
        ylims!(metric_ax, kwargs_dict[:ylims_metric])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    map(hidexdecorations!, [incax, reffax])

    if haskey(kwargs_dict, :plottitle)
        Label(
            fig[0, :, Top()],
            kwargs_dict[:plottitle],
        )
    end
    Legend(
        fig[1, 2],
        [PolyElement(; color = col) for col in Reff_colormap],
        ["Reff < 1", "Reff >= 1"],
        "Reff Status",
    )

    Legend(
        fig[2, 2],
        [PolyElement(; color = col) for col in outbreak_colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    Label(
        fig[0, :],
        "EWS Metric: $(string(ewsmetric)), $(method), Lag: $(lag), Bandwidth: $(bandwidth), Aggregation: $(aggregation) days\nEWS calculated for shaded region",
    )

    rowsize!(fig.layout, 0, Relative(0.05))

    return fig
end

function Reff_ews_plot(
    incidencearr,
    Reffarr,
    Reff_thresholds_vec,
    ewsmetric_sa::T1,
    spaero_ewsmetric_sa::T2,
    ewsmetric::S,
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
    metric_color = Makie.wong_colors()[1],
    spaero_metric_color = Makie.wong_colors()[3],
    linewidth = 3,
    hlinewidth = 2,
    kwargs...,
) where {T1<:StructArray,T2<:StructArray,S<:Symbol}
    @unpack trange = timeparams
    times = collect(trange) ./ 365

    kwargs_dict = Dict(kwargs)

    outbreak_status_vec = zeros(Int64, length(times))
    for (lower, upper, periodsum, outbreakstatus) in
        eachrow(thresholdsarr[sim])
        outbreak_status_vec[lower:upper] .= outbreakstatus
    end

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds_vec[sim])
        Reff_above_one[lower:upper] .= 1
    end

    ewsmetric_vec = getproperty(ewsmetric_sa, ewsmetric)[sim]
    replace!(ewsmetric_vec, Inf => NaN)
    spaero_ewsmetric_vec = getproperty(spaero_ewsmetric_sa, ewsmetric)[sim]
    replace!(spaero_ewsmetric_vec, Inf => NaN)

    fig = Figure()
    reffax = Axis(
        fig[1, 1]; ylabel = "Reff"
    )
    incax = Axis(
        fig[2, 1]; ylabel = "Incidence"
    )
    metric_ax = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = String(ewsmetric)
    )
    spaero_ax = Axis(
        fig[4, 1]; xlabel = "Time (years)",
        ylabel = "SPAERO " * String(ewsmetric),
    )

    linkxaxes!(incax, reffax, metric_ax, spaero_ax)

    line_and_hline!(
        reffax,
        times,
        @view(Reffarr[:, sim]),
        1;
        color = Reff_above_one,
        colormap = Reff_colormap,
        linewidth = linewidth,
        hlinewidth = hlinewidth,
    )

    line_and_hline!(
        incax,
        times,
        @view(incidencearr[:, 1, sim]),
        threshold;
        color = @view(outbreak_status_vec[:, 2]),
        colormap = outbreak_colormap,
        linewidth = linewidth,
        hlinewidth = hlinewidth,
    )

    lines!(
        metric_ax,
        times,
        ewsmetric_vec;
        color = metric_color,
        linewidth = linewidth,
    )

    lines!(
        spaero_ax,
        times,
        spaero_ewsmetric_vec;
        color = spaero_metric_color,
        linewidth = linewidth,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [incax, reffax, metric_ax],
        )
    end

    if haskey(kwargs_dict, :ylims_metric)
        ylims!(metric_ax, kwargs_dict[:ylims_metric])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    map(hidexdecorations!, [incax, reffax])

    if haskey(kwargs_dict, :plottitle)
        Label(
            fig[0, :, Top()],
            kwargs_dict[:plottitle],
        )
    end
    Legend(
        fig[1, 2],
        [PolyElement(; color = col) for col in Reff_colormap],
        ["Reff < 1", "Reff >= 1"],
        "Reff Status",
    )

    Legend(
        fig[2, 2],
        [PolyElement(; color = col) for col in outbreak_colormap],
        ["Not Outbreak", "Outbreak"],
        "True\nOutbreak Status",
    )

    rowsize!(fig.layout, 0, Relative(0.05))

    return fig
end

function Reff_ews_plot(
    incidencearr,
    Reffarr,
    Reff_thresholds_vec,
    ewsmetric_sa::T1,
    spaero_ewsmetric_sa::T2,
    ewsmetric::S,
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
    metric_color = Makie.wong_colors()[1],
    spaero_metric_color = Makie.wong_colors()[3],
    linewidth = 3,
    hlinewidth = 2,
    aggregation = 1,
    kwargs...,
) where {T1<:StructArray,T2<:DataFrame,S<:Symbol}
    @unpack trange = timeparams
    times = collect(trange) ./ 365

    kwargs_dict = Dict(kwargs)

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds_vec[sim])
        Reff_above_one[lower:upper] .= 1
    end
    Reff_above_one = aggregate_timeseries(Reff_above_one, aggregation)

    ewsmetric_vec = getproperty(ewsmetric_sa, ewsmetric)[sim]
    replace!(ewsmetric_vec, Inf => NaN)
    spaero_ewsmetric_vec = getproperty(spaero_ewsmetric_sa, ewsmetric)
    replace!(spaero_ewsmetric_vec, Inf => NaN)

    fig = Figure()
    reffax = Axis(
        fig[1, 1]; ylabel = "Reff"
    )
    incax = Axis(
        fig[2, 1]; ylabel = "Incidence"
    )
    metric_ax = Axis(
        fig[3, 1]; xlabel = "Time (years)", ylabel = String(ewsmetric)
    )
    diff_ax = Axis(
        fig[4, 1]; xlabel = "Time (years)", ylabel = "Diff between CSD & SPAERO"
    )

    linkxaxes!(incax, reffax, metric_ax, diff_ax)
    aggregated_vec_length = size(Reffarr, 1) ÷ aggregation * aggregation
    times = times[1:aggregation:aggregated_vec_length]

    line_and_hline!(
        reffax,
        times,
        @view(
            Reffarr[
                1:aggregation:aggregated_vec_length,
                sim,
            ]
        ),
        1;
        linewidth = linewidth,
        hlinewidth = hlinewidth,
    )

    line_and_hline!(
        incax,
        times,
        aggregate_timeseries(@view(incidencearr[:, 1, sim]), aggregation),
        threshold;
        linewidth = linewidth,
        hlinewidth = hlinewidth,
    )

    lines!(
        metric_ax,
        times,
        ewsmetric_vec;
        color = metric_color,
        linewidth = linewidth,
    )

    lines!(
        metric_ax,
        times,
        spaero_ewsmetric_vec;
        color = spaero_metric_color,
        linewidth = linewidth,
    )

    lines!(
        diff_ax,
        times,
        ewsmetric_vec .- spaero_ewsmetric_vec;
        color = :black,
        linewidth = linewidth,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [incax, reffax, metric_ax],
        )
    end

    if haskey(kwargs_dict, :ylims_metric)
        ylims!(metric_ax, kwargs_dict[:ylims_metric])
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    map(hidexdecorations!, [incax, reffax])

    if haskey(kwargs_dict, :plottitle)
        Label(
            fig[0, :, Top()],
            kwargs_dict[:plottitle],
        )
    end

    Legend(
        fig[3, 2],
        [
            PolyElement(; color = col) for
            col in [metric_color, spaero_metric_color]
        ],
        ["CSD Noise Metric", "SPAERO CSD Noise Metric"],
        "CSD Noise Metric",
    )

    Legend(
        fig[4, 2],
        [PolyElement(; color = :black)],
        ["Diff between CSD & SPAERO"],
        "Diff between CSD & SPAERO",
    )

    return fig
end

function simulation_tau_distribution(
    tested_ews_vec::T1,
    actual_ews_vec::T1,
    tau_metric;
    plottitle = string(tau_metric) * " Distribution",
) where {T1<:Vector{<:Union{<:Missing,<:EWSMetrics}}}
    if split(tau_metric, "_")[end] != "tau"
        tau_metric *= "_tau"
    end
    return simulation_tau_distribution(
        convert(Vector{EWSMetrics}, tested_ews_vec),
        convert(Vector{EWSMetrics}, actual_ews_vec),
        tau_metric;
        plottitle = plottitle,
    )
end

function simulation_tau_distribution(
    tested_ews_vec::T1,
    actual_ews_vec::T1,
    tau_metric;
    plottitle = string(tau_metric) * " Distribution",
) where {T1<:Vector{<:EWSMetrics}}
    tested_ews_sa = StructArray(tested_ews_vec)
    actual_ews_sa = StructArray(actual_ews_vec)

    if split(tau_metric, "_")[end] != "tau"
        tau_metric *= "_tau"
    end

    return simulation_tau_distribution(
        tested_ews_sa,
        actual_ews_sa,
        tau_metric;
        plottitle = plottitle,
    )
end

function simulation_tau_distribution(
    tested_ews::T1,
    actual_ews::T1,
    tau_metric;
    plottitle = string(tau_metric) * " Distribution",
) where {T1<:StructArray{EWSMetrics}}
    if split(tau_metric, "_")[end] != "tau"
        tau_metric *= "_tau"
    end

    tau = filter(!isnan, getproperty(tested_ews, Symbol(tau_metric)))
    actual_tau = filter(!isnan, getproperty(actual_ews, Symbol(tau_metric)))

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Kendall's Tau",
        ylabel = "Count",
    )

    hist!(
        ax,
        tau;
        color = (Makie.wong_colors()[1], 0.5),
        bins = 20,
        label = "Tested EWS",
    )

    hist!(
        ax,
        actual_tau;
        color = (Makie.wong_colors()[2], 0.5),
        bins = 20,
        label = "Actual EWS",
    )

    vlines!(
        ax,
        mean(tau);
        color = Makie.wong_colors()[1],
        linestyle = :dash,
        linewidth = 4,
    )

    vlines!(
        ax,
        mean(actual_tau);
        color = :grey20,
        linewidth = 4,
    )

    text!(
        ax,
        mean(actual_tau) + 0.005,
        length(actual_tau) / 10;
        text = "Actual Tau (Mean)",
        justification = :left,
    )

    lplots = filter(Makie.get_plots(ax)) do plot
        haskey(plot.attributes, :label)
    end

    if !isempty(lplots)
        Legend(fig[1, 2], ax, "Aggregation (wks)")
    end
    return fig
end
