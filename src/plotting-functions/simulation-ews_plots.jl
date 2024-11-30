using GLMakie
using UnPack: @unpack
using NaNMath: NaNMath
using StructArrays: StructArray
using DataFrames

function Reff_ews_plot(
    incvec,
    null_incvec,
    testvec,
    null_testvec,
    noisevec,
    Reffvec,
    null_Reffvec,
    Reff_thresholds,
    null_Reff_thresholds,
    ewsmetrics,
    null_ewsmetrics,
    ewsmetric::S,
    thresholdsarr,
    null_thresholdsarr,
    exceeds_thresholds_vec,
    null_exceeds_thresholds_vec,
    thresholds_percentile_vec,
    null_thresholds_percentile_vec,
    detection_index,
    null_detection_index,
    timeparams;
    plottitle = "",
    outbreak_colormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    outbreak_color_labels = ["Not Outbreak", "Outbreak"],
    threshold = 5,
    metric_color = Makie.wong_colors()[1],
    rowsize = Makie.Relative(0.03),
    kwargs...,
) where {S<:Symbol}
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @unpack ews_specification = ewsmetrics
    @unpack aggregation, bandwidth, lag, method = ews_specification
    method = split(string(method), "::")[1]

    min_Reff, max_Reff = extrema(vcat(Reffvec, null_Reffvec))
    ylims_Reff = (min_Reff - 0.1, max_Reff + 0.1)

    function calculate_ylims(
        ylims, kwargs_dict, detection_vec, null_detection_vec; link = true
    )
        if !link
            return (nothing, nothing)
        end
        if !haskey(kwargs_dict, ylims)
            min_metric, max_metric = extrema(
                vcat(
                    filter(!isnan, detection_vec),
                    filter(!isnan, null_detection_vec),
                ),
            )
            return (min_metric - 0.1, max_metric + 0.1)
        end
        return kwargs_dict[ylims]
    end

    ylims_metric = calculate_ylims(
        :ylims_metric,
        kwargs_dict,
        getproperty(ewsmetrics, ewsmetric),
        getproperty(null_ewsmetrics, ewsmetric),
    )

    ylims_Reff = calculate_ylims(
        :ylims_Reff,
        kwargs_dict,
        Reffvec,
        null_Reffvec,
    )

    ylims_inc = calculate_ylims(
        :ylims_inc,
        kwargs_dict,
        incvec,
        null_incvec;
        link = false,
    )

    fig = Figure()
    gl_detection = fig[1, 1] = GridLayout()
    gl_null = fig[1, 2] = GridLayout()

    Reff_ews_plot_facet!(
        gl_detection,
        incvec,
        testvec,
        noisevec,
        Reffvec,
        Reff_thresholds,
        ewsmetrics,
        ewsmetric,
        thresholdsarr,
        exceeds_thresholds_vec,
        detection_index,
        timeparams;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        outbreak_colormap = outbreak_colormap,
        outbreak_color_labels = outbreak_color_labels,
        threshold = threshold,
        metric_color = metric_color,
        ylims_Reff = ylims_Reff,
        ylims_inc = ylims_inc,
        ylims_metric = ylims_metric,
        thresholds_percentile_vec = thresholds_percentile_vec,
        kwargs...,
    )

    Reff_ews_plot_facet!(
        gl_null,
        null_incvec,
        null_testvec,
        noisevec,
        null_Reffvec,
        null_Reff_thresholds,
        null_ewsmetrics,
        ewsmetric,
        null_thresholdsarr,
        null_exceeds_thresholds_vec,
        null_detection_index,
        timeparams;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        outbreak_colormap = outbreak_colormap,
        outbreak_color_labels = outbreak_color_labels,
        threshold = threshold,
        metric_color = metric_color,
        ylims_Reff = ylims_Reff,
        ylims_inc = ylims_inc,
        ylims_metric = ylims_metric,
        thresholds_percentile_vec = null_thresholds_percentile_vec,
        kwargs...,
    )
    if plottitle == ""
        plottitle = "EWS Metric: $(string(ewsmetric)), $(method), Lag: $(lag), Bandwidth: $(bandwidth), Aggregation: $(aggregation)\nEWS calculated for shaded region"
    end

    Label(
        fig[0, :],
        plottitle,
    )

    rowsize!(fig.layout, 0, rowsize)
    # colsize!(fig.layout, 1, Relative(1))

    return fig
end

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
    rowsize = Makie.Relative(0.03),
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

    ewsmetric_endpoint = length(
        getproperty(ewsmetrics, propertynames(ewsmetrics)[2])
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

    add_ews_lines!(
        metric_ax,
        ewsmetrics,
        ewsmetric,
        ewsmetric_endpoint,
        detection_index,
        exceeds_thresholds_vec,
        times;
        metric_color = Makie.wong_colors()[1],
    )

    add_ews_lines!(
        null_metric_ax,
        null_ewsmetrics,
        ewsmetric,
        ewsmetric_endpoint,
        null_detection_index,
        null_exceeds_thresholds_vec,
        times;
        metric_color = Makie.wong_colors()[1],
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
        "EWS Metric: $(string(ewsmetric)), $(method), Lag: $(lag), Bandwidth: $(bandwidth), Aggregation: $(aggregation)\nEWS calculated for shaded region",
    )

    rowsize!(fig.layout, 0, rowsize)

    return fig
end

function add_ews_lines!(
    ax,
    ewsmetrics,
    ewsmetric,
    ewsmetric_endpoint,
    detection_index,
    exceeds_thresholds_vec,
    times;
    metric_color = Makie.wong_colors()[1],
    plot_tau = true,
    plot_detection_index = true,
    kwargs...,
)
    ewsmetric_vec = getproperty(ewsmetrics, ewsmetric)
    ewsmetric_vec = vcat(
        ewsmetric_vec,
        fill(NaN, length(times) - ewsmetric_endpoint),
    )
    ewsmetric_tau = getproperty(ewsmetrics, Symbol(String(ewsmetric) * "_tau"))

    lines!(
        ax,
        times,
        ewsmetric_vec;
        color = metric_color,
        linewidth = 3,
    )
    if !isnothing(detection_index)
        exceeds_threshold_indices = findall(
            x -> x == true, exceeds_thresholds_vec
        )

        if plot_detection_index
            vlines!(ax, times[detection_index]; color = :black)
        end
        scatter!(
            ax,
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

    if plot_tau
        text!(
            ax,
            times[ewsmetric_endpoint] + 0.5,
            ewsmetric_tau_yvalue;
            text = "τ = $(round(ewsmetric_tau; digits = 2))",
            justification = :right,
        )
    end

    return nothing
end

function Reff_ews_plot(
    incvec,
    testvec,
    noisevec,
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
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    outbreak_color_labels = ["Not Outbreak", "Outbreak"],
    threshold = 5,
    metric_color = Makie.wong_colors()[1],
    rowsize = Makie.Relative(0.03),
    kwargs...,
) where {S<:Symbol}
    @unpack trange = timeparams
    @unpack aggregation, bandwidth, lag, method = ewsmetrics.ews_specification
    method = split(string(method), "::")[1]

    fig = Figure()
    gl = fig[1, 1] = GridLayout()

    Reff_ews_plot_facet!(
        gl,
        incvec,
        testvec,
        noisevec,
        Reffvec,
        Reff_thresholds,
        ewsmetrics,
        ewsmetric,
        thresholdsarr,
        exceeds_thresholds_vec,
        detection_index,
        timeparams;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        outbreak_colormap = outbreak_colormap,
        outbreak_color_labels = outbreak_color_labels,
        threshold = threshold,
        metric_color = metric_color,
        kwargs...,
    )

    return fig
end

function Reff_noise_inc_test_plot(
    testvec,
    incvec,
    thresholdsarr,
    Reffvec,
    Reff_thresholds,
    noisevec,
    timeparams,
    ewsmetrics;
    outbreak_colormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    legends = true,
    # size = (1300, 800 / 2),
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @unpack trange = timeparams
    @unpack aggregation = ewsmetrics.ews_specification
    ewsmetric_endpoint = length(
        getproperty(ewsmetrics, propertynames(ewsmetrics)[2])
    )

    aggregation_int = Dates.days(aggregation)

    aggregated_vec_length = length(Reffvec)
    times =
        collect(1:aggregation_int:(aggregated_vec_length * aggregation_int)) ./
        365

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

    fig = Figure()

    reffax = Axis(
        fig[1, 1]; ylabel = L"R_{E}"
    )

    Reff_facet!(
        reffax,
        Reffvec,
        Reff_thresholds,
        Reff_above_one,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        legends = legends,
        kwargs...,
    )

    hidexdecorations!(reffax)

    incax = Axis(
        fig[2, 1]; xlabel = "Time (Years)", ylabel = "Incidence"
    )

    noise_facet!(
        incax,
        noisevec,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        kwargs...,
    )

    lines!(
        incax,
        times,
        incvec;
        color = outbreak_status_vec,
        colormap = outbreak_colormap,
    )

    lines!(
        incax,
        times,
        testvec;
        color = :black,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [reffax, incax],
        )
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    if haskey(kwargs_dict, :ylims_Reff)
        ylims!(reffax, kwargs_dict[:ylims_Reff])
    end

    return fig
end

function Reff_noise_inc_plot(
    incvec,
    thresholdsarr,
    Reffvec,
    Reff_thresholds,
    noisevec,
    timeparams,
    ewsmetrics;
    outbreak_colormap = [
        BASE_COLOR, OUTBREAK_COLOR
    ],
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    legends = true,
    # size = (1300, 800 / 2),
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @unpack trange = timeparams
    @unpack aggregation = ewsmetrics.ews_specification
    ewsmetric_endpoint = length(
        getproperty(ewsmetrics, propertynames(ewsmetrics)[2])
    )

    aggregation_int = Dates.days(aggregation)

    aggregated_vec_length = length(Reffvec)
    times =
        collect(1:aggregation_int:(aggregated_vec_length * aggregation_int)) ./
        365

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

    fig = Figure()

    reffax = Axis(
        fig[1, 1]; ylabel = L"R_{E}"
    )

    Reff_facet!(
        reffax,
        Reffvec,
        Reff_thresholds,
        Reff_above_one,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        legends = legends,
        kwargs...,
    )

    hidexdecorations!(reffax)

    incax = Axis(
        fig[2, 1]; xlabel = "Time (Years)", ylabel = "Incidence"
    )

    noise_facet!(
        incax,
        noisevec,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        kwargs...,
    )

    lines!(
        incax,
        times,
        incvec;
        color = outbreak_status_vec,
        colormap = outbreak_colormap,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [reffax, incax],
        )
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    if haskey(kwargs_dict, :ylims_Reff)
        ylims!(reffax, kwargs_dict[:ylims_Reff])
    end

    return fig
end

function Reff_noise_plot(
    Reffvec,
    Reff_thresholds,
    noisevec,
    timeparams,
    ewsmetrics;
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    legends = true,
    # size = (1300, 800 / 2),
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @unpack trange = timeparams
    @unpack aggregation = ewsmetrics.ews_specification
    ewsmetric_endpoint = length(
        getproperty(ewsmetrics, propertynames(ewsmetrics)[2])
    )

    aggregation_int = Dates.days(aggregation)

    aggregated_vec_length = length(Reffvec)
    times =
        collect(1:aggregation_int:(aggregated_vec_length * aggregation_int)) ./
        365

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds)
        Reff_above_one[lower:1:upper] .= 1
    end

    fig = Figure()

    reffax = Axis(
        fig[1, 1]; ylabel = L"R_{E}"
    )

    Reff_facet!(
        reffax,
        Reffvec,
        Reff_thresholds,
        Reff_above_one,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        legends = legends,
        kwargs...,
    )

    hidexdecorations!(reffax)

    incax = Axis(
        fig[2, 1]; xlabel = "Time (Years)", ylabel = "Incidence"
    )

    noise_facet!(
        incax,
        noisevec,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        kwargs...,
    )

    if haskey(kwargs_dict, :xlims)
        map(
            ax -> xlims!(ax, kwargs_dict[:xlims]),
            [reffax, incax],
        )
    end

    if haskey(kwargs_dict, :ylims_inc)
        ylims!(incax, kwargs_dict[:ylims_inc])
    end

    if haskey(kwargs_dict, :ylims_Reff)
        ylims!(reffax, kwargs_dict[:ylims_Reff])
    end

    return fig
end

function noise_facet!(
    noiseax,
    noisevec,
    times,
    aggregation_int,
    ewsmetric_endpoint;
    ewscolor_vspan = :lightgrey,
    ewscolor_vspan_alpha = 0.5,
    kwargs...,
)
    add_ews_vspan!(
        noiseax,
        aggregation_int,
        ewsmetric_endpoint;
        color = ewscolor_vspan,
        alpha = ewscolor_vspan_alpha,
    )

    band!(
        noiseax,
        times,
        repeat([0], length(times)),
        noisevec;
        color = (:darkred, 0.4),
    )

    return nothing
end

function add_ews_vspan!(
    ax,
    aggregation_int,
    ewsmetric_endpoint;
    color = :lightgrey,
    alpha = 0.5,
)
    vspan!(
        ax, 0, (ewsmetric_endpoint - 1) * aggregation_int / 365;
        color = (color, alpha),
    )

    return nothing
end

function Reff_plot(
    Reffvec,
    Reff_thresholds,
    timeparams,
    ewsmetrics;
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    legends = true,
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @unpack trange = timeparams
    @unpack aggregation = ewsmetrics.ews_specification
    ewsmetric_endpoint = length(
        getproperty(ewsmetrics, propertynames(ewsmetrics)[2])
    )

    aggregation_int = Dates.days(aggregation)

    aggregated_vec_length = length(Reffvec)
    times =
        collect(1:aggregation_int:(aggregated_vec_length * aggregation_int)) ./
        365

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds)
        Reff_above_one[lower:1:upper] .= 1
    end

    fig = Figure()

    reffax = Axis(
        fig[1, 1]; xlabel = "Time (Years)", ylabel = L"R_{E}"
    )

    Reff_facet!(
        reffax,
        Reffvec,
        Reff_thresholds,
        Reff_above_one,
        times,
        aggregation_int,
        ewsmetric_endpoint;
        Reff_colormap = Reff_colormap,
        Reff_color_labels = Reff_color_labels,
        legends = legends,
        kwargs...,
    )

    if haskey(kwargs_dict, :xlims)
        xlims!(reffax, kwargs_dict[:xlims])
    end

    if haskey(kwargs_dict, :ylims_Reff)
        ylims!(reffax, kwargs_dict[:ylims_Reff])
    end

    return fig
end

function Reff_facet!(
    reffax,
    Reffvec,
    Reff_thresholds,
    Reff_above_one,
    times,
    aggregation_int,
    ewsmetric_endpoint;
    Reff_colormap = [
        BASE_COLOR, REFF_GT_ONE_COLOR
    ],
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    legends = true,
    ewscolor_vspan = :lightgrey,
    ewscolor_vspan_alpha = 0.5,
    kwargs...,
)
    add_ews_vspan!(
        reffax,
        aggregation_int,
        ewsmetric_endpoint;
        color = ewscolor_vspan,
        alpha = ewscolor_vspan_alpha,
    )

    if legends
        if length(unique(Reff_above_one)) > 1
            Legend(
                reffax[1, 2],
                [PolyElement(; color = col) for col in Reff_colormap],
                Reff_color_labels,
                "Reff Status",
            )
        else
            Reff_colormap = [Reff_colormap[1]]
            Reff_color_labels = [Reff_color_labels[1]]
        end
    end

    line_and_hline!(
        reffax,
        times,
        Reffvec,
        1;
        color = Reff_above_one,
        colormap = Reff_colormap,
        linewidth = 3,
    )

    return nothing
end

function Reff_ews_plot_facet!(
    gl,
    incvec,
    testvec,
    noisevec,
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
    Reff_color_labels = ["Reff < 1", "Reff >= 1"],
    outbreak_color_labels = ["Not Outbreak", "Outbreak"],
    threshold = 5,
    metric_color = Makie.wong_colors()[1],
    legends = true,
    plot_tau = true,
    plot_detection_index = true,
    kwargs...,
) where {S<:Symbol}
    @unpack trange = timeparams
    @unpack aggregation, method =
        ewsmetrics.ews_specification
    method = split(string(method), "::")[1]

    aggregation_int = Dates.days(aggregation)

    aggregated_vec_length = length(Reffvec)
    times =
        collect(1:aggregation_int:(aggregated_vec_length * aggregation_int)) ./
        365

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

    ewsmetric_endpoint = length(
        getproperty(ewsmetrics, propertynames(ewsmetrics)[2])
    )

    reffax = Axis(
        gl[1, 1]; ylabel = L"R_{E}"
    )
    incax = Axis(
        gl[2, 1]; ylabel = "Incidence"
    )
    metric_ax = Axis(
        gl[3, 1]; xlabel = "Time (years)", ylabel = titlecase(String(ewsmetric))
    )

    linkxaxes!(incax, reffax, metric_ax)

    map(
        ax ->
            vspan!(
                ax, 0, (ewsmetric_endpoint - 1) * aggregation_int / 365;
                color = (:lightgrey, 0.5),
            ),
        [incax, reffax, metric_ax],
    )

    if legends
        if length(unique(Reff_above_one)) > 1
            Legend(
                gl[1, 2],
                [PolyElement(; color = col) for col in Reff_colormap],
                Reff_color_labels,
                "Reff Status",
            )
        else
            Reff_colormap = [Reff_colormap[1]]
            Reff_color_labels = [Reff_color_labels[1]]
        end

        if length(unique(outbreak_status_vec)) > 1
            Legend(
                gl[2, 2],
                [PolyElement(; color = col) for col in outbreak_colormap],
                outbreak_color_labels,
                "True\nOutbreak Status",
            )
        else
            outbreak_colormap = [outbreak_colormap[1]]
            outbreak_color_labels = [outbreak_color_labels[1]]
        end
    end

    line_and_hline!(
        reffax,
        times,
        Reffvec,
        1;
        color = Reff_above_one,
        colormap = Reff_colormap,
        linewidth = 3,
    )

    band!(
        incax,
        times,
        repeat([0], length(times)),
        noisevec;
        color = (:darkred, 0.4),
        # linewidth = 5,
    )

    outbreak_threshold_hline = aggregation_int == 1 ? threshold : NaN
    line_and_hline!(
        incax,
        times,
        incvec,
        outbreak_threshold_hline;
        color = outbreak_status_vec,
        colormap = outbreak_colormap,
        linewidth = 3,
    )

    lines!(
        incax,
        times,
        testvec;
        linewidth = 3,
        color = :black,
    )

    add_ews_lines!(
        metric_ax,
        ewsmetrics,
        ewsmetric,
        ewsmetric_endpoint,
        detection_index,
        exceeds_thresholds_vec,
        times;
        metric_color = metric_color,
        plot_tau = plot_tau,
        plot_detection_index = plot_detection_index,
    )

    if haskey(kwargs_dict, :burnin_vline)
        burnin_vline = if typeof(kwargs_dict[:burnin_vline]) == Dates.Year
            Dates.value(kwargs_dict[:burnin_vline])
        else
            Dates.days(kwargs_dict[:burnin_vline]) / 365
        end
        map(
            ax -> vlines!(
                ax,
                burnin_vline;
                color = :darkred,
                linestyle = :dash,
            ),
            [incax, reffax, metric_ax],
        )
    end

    if haskey(kwargs_dict, :thresholds_percentile_vec)
        lines!(
            metric_ax,
            times[1:ewsmetric_endpoint],
            kwargs_dict[:thresholds_percentile_vec];
            color = Makie.wong_colors()[2],
            linewidth = 3,
        )
    end

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

    if haskey(kwargs_dict, :ylims_Reff)
        ylims!(reffax, kwargs_dict[:ylims_Reff])
    end

    map(hidexdecorations!, [incax, reffax])

    if haskey(kwargs_dict, :plottitle)
        Label(
            gl[0, :, Top()],
            kwargs_dict[:plottitle],
        )
    end

    return nothing
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
    rowsize = Makie.Relative(0.03),
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

    rowsize!(fig.layout, 0, rowsize)

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
