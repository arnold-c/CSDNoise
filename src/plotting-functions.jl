# module PlottingFunctions
#
# export seircolors, seir_state_labels, create_sir_plot, draw_sir_plot,
#     sir_quantiles_array_base_plot, create_sir_quantiles_plot, outbreakcols,
#     detect_outbreak_plot, visualize_ensemble_noise, incidence_testing_plot,
#     testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot

using GLMakie
using AlgebraOfGraphics
using ColorSchemes
using UnPack
using DataFrames
using FLoops
using LaTeXStrings
using NaNMath: NaNMath
using Dates
using Debugger: Debugger

seircolors = ["dodgerblue4", "green", "firebrick3", "chocolate2", "purple"]
seir_state_labels = ["S", "E", "I", "R", "N"]

ACCURACY_COLOR = "#004643"
DAILY_SENSITIVITY_COLOR = "#2A3965"
DAILY_SPECIFICITY_COLOR = "#C31D60"
DAILY_PPV_COLOR = "#22D37D"
DAILY_NPV_COLOR = "#6f366bff"
DETECTION_DELAY_COLOR = "#AE560A"
N_ALERTS_PER_OUTBREAK_COLOR = "#86B1A3"
N_FALSE_ALERTS_COLOR = "#D06778"
N_ALERTS_COLOR = "#00857E"
N_OUTBREAKS_COLOR = "#F4A157"
N_MISSED_OUTBREAKS_COLOR = "#5E5C6C"
PERC_OUTBREAKS_DETECTED_COLOR = "#F0780F"
PERC_OUTBREAKS_MISSED_COLOR = "#3A3842"
PERC_ALERTS_CORRECT_COLOR = "#004643"
PERC_ALERTS_FALSE_COLOR = "#852938"

function create_sir_plot(sol_df; labels = ["S", "I", "R", "N"], annual = annual)
    time_function(t) = annual == true ? t / 365.0 : t
    if annual == true
        time_label = "Time (years)"
    else
        time_label = "Time (days)"
    end

    return data(sol_df) *
           mapping(
               :time => time_function => time_label, :Number;
               color = :State => sorter(labels...),
           ) *
           visual(Lines; linewidth = 4)
end

function draw_sir_plot(
    sir_plot; colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    xlims = xlims, ylims = ylims,
)
    return draw(
        sir_plot;
        palettes = (; color = colors),
        axis = (; limits = (xlims, ylims)),
    )
end

function draw_sir_plot(
    sol_df::DataFrame;
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"],
    annual = false,
    xlims = nothing,
    ylims = nothing,
)
    return draw_sir_plot(
        create_sir_plot(sol_df; labels = labels, annual = annual);
        colors = colors,
        xlims = xlims,
        ylims = ylims,
    )
end

function sir_quantiles_array_base_plot(
    sim_quantiles, lower_index, med_index, upper_index, timeparams, colors,
    labels,
    annual; xlims, ylims, caption,
)
    times = timeparams.trange
    xlab = "Time (days)"
    if annual == true
        times = times ./ 365
        xlab = "Time (years)"
    end

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = xlab, ylabel = "Number")

    # Medians
    map(
        state -> lines!(
            ax,
            times,
            sim_quantiles[med_index, :, state];
            color = colors[state],
            linewidth = 2,
            label = labels[state],
        ),
        eachindex(labels),
    )

    map(
        state -> band!(
            ax,
            times,
            sim_quantiles[lower_index, :, state],
            sim_quantiles[upper_index, :, state];
            color = (colors[state], 0.5),
        ),
        eachindex(labels),
    )

    if xlims != false
        xlims!(ax, xlims)
    end

    if ylims != false
        ylims!(ax, ylims)
    end

    Legend(fig[1, 2], ax, "State")

    if caption != false
        Label(fig[2, :, Bottom()], caption)
        rowsize!(fig.layout, 1, Relative(0.98))
    end

    return fig
end

function create_sir_quantiles_plot(
    sim_quantiles = sim_quantiles;
    timeparams,
    colors = ["dodgerblue4", "firebrick3", "chocolate2", "purple"],
    labels = ["S", "I", "R", "N"],
    annual = false,
    xlims = false,
    ylims = false,
    caption = false,
)
    med_index = 2
    lower_index = 1
    upper_index = 3

    return sir_quantiles_array_base_plot(
        sim_quantiles, lower_index, med_index, upper_index, timeparams,
        colors,
        labels,
        annual;
        xlims = xlims,
        ylims = ylims,
        caption = caption,
    )
end

function Reff_plot(
    incidencearr,
    Reffarr,
    Reff_thresholds_vec,
    thresholdsarr,
    timeparams;
    sim = 1,
    outbreak_colormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    Reff_colormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
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
    colormap = [N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR],
    threshold = 5,
    kwargs...,
)
    @unpack trange = timeparams
    times = collect(trange) ./ 365
    kwargs_dict = Dict(kwargs)

    period_sum_arr = zeros(Int64, length(times), 2)
    for (lower, upper, periodsum, outbreakstatus) in
        eachrow(thresholdsarr[1])
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
        ensemblearr[:, 2, 1];
        color = period_sum_arr[:, 2],
        colormap = colormap,
    )
    lines!(
        ax_inc,
        times,
        incidencearr[:, 1, 1];
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

function visualize_ensemble_noise(
    ensemble_noise_arr, poisson_noise_prop_vec, timespecification, noisedir;
    xlabel = "Time (years)", ylabel = "Noise Incidence",
)
    times = collect(timespecification.trange) ./ 365
    meanline = vec(mean(ensemble_noise_arr; dims = 2))
    dailymean = NaNMath.mean(meanline)
    poisson_noise_prop = NaNMath.mean(poisson_noise_prop_vec)

    poisson_noise_arr = ensemble_noise_arr .* poisson_noise_prop_vec
    poisson_noise_daily_mean = NaNMath.mean(poisson_noise_arr)
    dynamic_noise_arr = ensemble_noise_arr .- poisson_noise_arr
    dynamic_noise_daily_mean = NaNMath.mean(dynamic_noise_arr)

    fig = Figure()
    ax = Axis(fig[2, 1]; xlabel = xlabel, ylabel = ylabel)

    for noise_sim in eachcol(ensemble_noise_arr)
        lines!(
            ax,
            times,
            noise_sim;
            color = (:gray, 0.2),
        )
    end

    lines!(
        ax,
        times,
        meanline;
        color = :black,
    )

    Label(
        fig[1, :],
        "Noise: $(noisedir), Total Daily Mean: $(round(dailymean, digits = 2)), Poisson Noise Proportion: $(round(poisson_noise_prop, digits = 2))\nDynamic Daily Mean: $(round(dynamic_noise_daily_mean, digits = 2)) Poisson Daily Mean: $(round(poisson_noise_daily_mean, digits = 2))",
    )

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, Relative(1))

    return fig
end

function incidence_testing_plot(
    incarr,
    noisearr,
    testingarr,
    test_movingvg_arr,
    detection_specification,
    timeparams;
    sim = 1,
    outbreakcolormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    alertcolormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    plottitle = "",
    kwargs...,
)
    times = collect(timeparams.trange) ./ 365
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
        inc_test_ax1, times, incarr[:, 1, sim];
        color = incarr[:, 3, sim],
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax2, times, incarr[:, 1, sim] .+ noisearr[:, sim];
        color = incarr[:, 3, sim],
        colormap = outbreakcolormap,
    )
    lines!(
        inc_test_ax3, times, testingarr[:, 5, sim];
        color = testingarr[:, 6, sim],
        colormap = alertcolormap,
    )
    lines!(
        inc_test_ax4, times, test_movingvg_arr[:, sim];
        color = testingarr[:, 6, sim],
        colormap = alertcolormap,
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
            5;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax1, inc_test_ax2],
    )

    map(
        ax -> hlines!(
            ax,
            detection_specification.alert_threshold;
            color = :black,
            linestyle = :dash,
            linewidth = 2,
        ),
        [inc_test_ax3, inc_test_ax4],
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

    Legend(
        inc_test_fig[3:4, 2],
        [PolyElement(; color = col) for col in alertcolormap],
        ["Not Outbreak", "Outbreak"],
        "Alert Status",
    )

    rowsize!(inc_test_fig.layout, 0, 5)
    colsize!(inc_test_fig.layout, 1, Relative(0.92))

    return inc_test_fig
end

function testing_plot(
    testingarr, timeparams; plottitle = "", sim1_num = 1, sim2_num = 25
)
    times = collect(timeparams.trange) ./ 365

    testing_fig = Figure()
    testing_grid = testing_fig[1, 1] = GridLayout()
    sim1_ax = Axis(
        testing_grid[1, 1];
        title = "Simulation $sim1_num",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )
    sim2_ax = Axis(
        testing_grid[2, 1];
        title = "Simulation $sim2_num",
        xlabel = "Time (years)",
        ylabel = "Number tested",
    )

    for (ind, label, col) in zip(
        (1, 2, 5),
        ("Infectious", "Noise", "Total Positive"),
        (:red, :blue, :black),
    )
        lines!(
            sim1_ax, times, testingarr[:, ind, sim1_num]; color = col,
            label = label,
        )
        lines!(
            sim2_ax, times, testingarr[:, ind, sim2_num]; color = col,
            label = label,
        )
    end

    linkxaxes!(sim1_ax, sim2_ax)

    hidexdecorations!(sim1_ax)

    Legend(
        testing_fig[:, 2],
        sim1_ax,
        "Type of Individual";
        orientation = :vertical,
    )

    Label(
        testing_fig[0, :, Top()],
        plottitle,
    )

    rowsize!(testing_fig.layout, 0, 5)
    colsize!(testing_fig.layout, 1, Relative(0.92))

    return testing_fig
end

function ensemble_outbreak_distribution_plot(testarr, infecarr; plottitle = "")
    outbreak_dist_fig = Figure()
    outbreak_dist_ax = Axis(
        outbreak_dist_fig[1, 1];
        xlabel = "Proportion of Time Series with Outbreak",
    )

    hist!(
        outbreak_dist_ax,
        vec(sum(@view(infecarr[:, 3, :]); dims = 1)) ./ size(infecarr, 1);
        bins = 0.0:0.01:0.7,
        color = (:blue, 0.5),
        strokecolor = :black,
        strokewidth = 1,
        normalization = :pdf,
        label = "True Outbreaks",
    )

    hist!(
        outbreak_dist_ax,
        vec(sum(@view(testarr[:, 6, :]); dims = 1)) ./ size(testarr, 1);
        bins = 0.0:0.01:0.7,
        color = (:red, 0.5),
        strokecolor = :black,
        strokewidth = 1,
        normalization = :pdf,
        label = "Tested Outbreaks",
    )

    Legend(outbreak_dist_fig[1, 2], outbreak_dist_ax, "Outbreak Proportion")

    Label(
        outbreak_dist_fig[0, :, Top()],
        plottitle,
    )

    rowsize!(outbreak_dist_fig.layout, 0, 5)
    colsize!(outbreak_dist_fig.layout, 1, Relative(0.92))

    return outbreak_dist_fig
end

function ensemble_OTChars_plot(
    OTChars,
    testspec,
    detectspec,
    plottingchars;
    plottitle = "",
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    binwidth = 10.0,
    legend = true,
    legendlabel = "# Outbreaks",
    meanlines = true,
    meanlabels = true,
    normalization = :none,
    kwargs...,
)
    OT_char_tuple = (;
        OT_chars = OTChars,
        ind_test_spec = testspec,
        outbreak_detect_spec = detectspec,
    )
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @pack! kwargs_dict = binwidth, xlabel, ylabel, legendlabel

    charvecs = map(
        chartuple -> reduce(
            vcat, getproperty(OT_char_tuple.OT_chars, chartuple.char)
        ),
        plottingchars,
    )

    fig = Figure()
    ax = Axis(fig[2, 1]; xlabel = xlabel, ylabel = ylabel)

    construct_single_OTchars_facet!(
        ax,
        charvecs,
        plottingchars,
        kwargs_dict;
        meanlines = meanlines,
        meanlabels = meanlabels,
        normalization = normalization,
    )

    if legend
        Legend(
            fig[2, 2],
            [
                PolyElement(; color = col) for
                col in map(chartuple -> chartuple.color, plottingchars)
            ],
            [chartuple.label for chartuple in plottingchars];
            label = legendlabel,
        )
    end

    Label(
        fig[1, :, Top()],
        plottitle,
    )

    rowsize!(fig.layout, 1, 5)

    return fig
end

function ensemble_outbreak_detect_diff_plot(OT_chars; binwidth = 1)
    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel = "Difference Between Actual - Detected Outbreaks"
    )

    difference = OT_chars.noutbreaks .- OT_chars.nalerts

    bins = minimum(difference):binwidth:maximum(difference)

    hist!(
        ax,
        difference;
        color = (:purple, 0.5),
        bins = bins,
    )

    vlines!(
        ax,
        mean(difference);
        color = :black,
        linestyle = :dash,
        linewidth = 4,
    )

    return fig
end

function singlescenario_test_positivity_plot(
    test_positivity_struct_vec; agg = :seven_day
)
    posoddsmatrix = reduce(
        hcat,
        map(array -> array[:, 1], getfield.(test_positivity_struct_vec, agg)),
    )
    avgpositivity = vec(mapslices(NaNMath.mean, posoddsmatrix; dims = 2))

    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel = "Time steps by $(agg)", ylabel = "Test Positivity"
    )
    lines!(ax, 1:length(avgpositivity), avgpositivity)

    return fig
end

function test_positivity_distribution_plot(
    test_positivity_struct_vec; agg = :seven_day, kwargs...
)
    df = @chain test_positivity_struct_vec begin
        getfield.(agg)
        reduce(vcat, _)
        DataFrame(Tables.table(_), [:positivity, :outbreak])
    end

    df[!, :outbreak] = string.(df[:, :outbreak])

    return draw(
        data(df) *
        mapping(
            :positivity => "Test Positivity aggregated by $(agg)"; kwargs...
        ) *
        histogram(; bins = 0.0:0.05:1.05), ;
        axis = (ylabel = "Count",),
    )
end

function save_compare_ensemble_OTchars_plot(
    char_struct_vec,
    columnfacetchar::Symbol,
    plottingchars,
    percent_clinic_tested;
    plotname,
    plotdirpath,
    plotformat = "png",
    size = (2200, 1200),
    columnfacetchar_label = "Alert Threshold",
    binwidth = 1.0,
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    legend = true,
    legendlabel = "Outbreak Chacteristic",
    meanlines = false,
    meanlabels = false,
    normalization = :none,
    force = false,
    kwargs...,
)
    mkpath(plotdirpath)

    plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

    if !isfile(plotpath) || force
        plot = compare_ensemble_OTchars_plots(
            char_struct_vec,
            columnfacetchar,
            plottingchars,
            percent_clinic_tested;
            columnfacetchar_label = columnfacetchar_label,
            binwidth = binwidth,
            xlabel = xlabel,
            ylabel = ylabel,
            legend = legend,
            legendlabel = legendlabel,
            meanlines = meanlines,
            meanlabels = meanlabels,
            normalization = normalization,
            kwargs...,
        )

        save(plotpath, plot; size = size)

        Makie.empty!(plot)
    end

    return nothing
end

function compare_ensemble_OTchars_plots(
    char_struct_vec,
    columnfacetchar::Symbol,
    plottingchars,
    percent_clinic_tested;
    columnfacetchar_label = "Alert Threshold",
    binwidth = 1.0,
    xlabel = "Alert Characteristic Value",
    ylabel = "Density",
    legend = true,
    legendlabel = "Outbreak Chacteristic",
    meanlines = false,
    meanlabels = false,
    normalization = :none,
    kwargs...,
)
    xs, ys = calculate_comparison_plot_facet_dims(
        char_struct_vec, columnfacetchar
    )
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    @pack! kwargs_dict = columnfacetchar_label,
    binwidth,
    xlabel,
    ylabel,
    legendlabel

    fig = Figure()

    construct_OTchars_facets!(
        fig,
        char_struct_vec,
        plottingchars,
        xs,
        ys,
        kwargs_dict;
        meanlines = meanlines,
        meanlabels = meanlabels,
        normalization = normalization,
    )

    if legend
        Legend(
            fig[:, end + 1],
            [
                PolyElement(; color = col) for
                col in map(chartuple -> chartuple.color, plottingchars)
            ],
            map(chartuple -> chartuple.label, plottingchars);
            label = legendlabel,
        )
    end

    Label(
        fig[1, :, Top()],
        "Perc Clinic Tested: $(percent_clinic_tested), Noise: $(getdirpath(char_struct_vec[end].noise_specification))",
    )

    unique_thresholds = unique(
        getfield.(
            getfield.(char_struct_vec, :outbreak_detect_spec), :alert_threshold
        ),
    )

    for (j, threshold) in pairs(unique_thresholds)
        Label(
            fig[2, j + 1, Top()],
            "Alert Threshold: $(threshold)",
        )
    end

    unique_test_specs = unique(getfield.(char_struct_vec, :ind_test_spec))

    for (i, test_spec) in pairs(unique_test_specs)
        Label(
            fig[i + 2, 1, Left()],
            "Sens: $(test_spec.sensitivity), Spec: $(test_spec.specificity),\nLag: $(test_spec.test_result_lag)";
            rotation = pi / 2,
        )
    end

    rowsize!(fig.layout, 1, 5)
    rowsize!(fig.layout, 2, 7)
    colsize!(fig.layout, 1, 7)

    return fig
end

function calculate_comparison_plot_facet_dims(
    char_struct_vec, facetchar
)
    xlength = length(
        Set(
            getfield.(
                getfield.(char_struct_vec, :outbreak_detect_spec),
                facetchar,
            ),
        ),
    )
    ylength = length(
        Set(getfield.(char_struct_vec, :ind_test_spec))
    )

    ys = repeat(1:ylength, xlength)
    xs = repeat(1:xlength; inner = ylength)

    return xs, ys
end

function construct_single_OTchars_facet!(
    ax,
    charvecs,
    plottingchars,
    kwargs_dict;
    meanlines = false,
    meanlabels = false,
    normalization = :pdf,
)
    @unpack binwidth = kwargs_dict

    if !haskey(kwargs_dict, :bins)
        bins = calculate_bins(charvecs, binwidth)
    else
        bins = kwargs_dict[:bins]
    end

    for charnumber in eachindex(plottingchars)
        if isempty(charvecs[charnumber])
            break
        end
        hist!(
            ax,
            charvecs[charnumber];
            bins = bins,
            color = plottingchars[charnumber].color,
            normalization = normalization,
        )

        if meanlines || meanlabels
            charmean = mean(charvecs[charnumber])
        end
        if meanlines
            vlines!(
                ax,
                charmean;
                color = :black,
                linestyle = :dash,
                linewidth = 4,
            )
        end
        if meanlabels
            hjust = 0
            vjust = 0
            if haskey(plottingchars[charnumber], :hjust)
                hjust = plottingchars[charnumber].hjust
            end
            if haskey(plottingchars[charnumber], :vjust)
                vjust = plottingchars[charnumber].vjust
            end
            text!(
                Point(charmean + hjust, 0 + vjust);
                text = "Mean ($(plottingchars[charnumber].label)):\n$(round(charmean, digits = 2))",
            )
        end
    end
end

function construct_OTchars_facets!(
    fig,
    char_struct_vec,
    plottingchars,
    xs,
    ys,
    kwargs_dict;
    meanlines = false,
    meanlabels = false,
    normalization = :pdf,
)
    number_tests = length(
        Set(getfield.(char_struct_vec, :ind_test_spec))
    )

    for (OT_char_tuple, x, y) in zip(char_struct_vec, xs, ys)
        skipped_plottingchar = 0

        charvecs = Vector{Vector{Union{Int64,Float64}}}(
            undef, length(plottingchars)
        )

        for (i, chartuple) in pairs(plottingchars)
            charvec = getproperty(OT_char_tuple.OT_chars, chartuple.char)

            if sum(isempty.(charvec)) == length(charvec)
                skipped_plottingchar += 1
                charvecs[i] = eltype(charvec)[]
            else
                charvecs[i] = reduce(vcat, charvec)
            end
        end

        @unpack xlabel, ylabel = kwargs_dict

        gl = fig[y + 2, x + 1] = GridLayout()
        ax = Axis(
            gl[1, 1];
            xlabel = xlabel,
            ylabel = ylabel,
        )

        if y < number_tests
            hidexdecorations!(ax; ticklabels = false, ticks = false)
        end

        hideydecorations!(ax)

        if skipped_plottingchar == length(plottingchars)
            continue
        end

        construct_single_OTchars_facet!(
            ax,
            charvecs,
            plottingchars,
            kwargs_dict;
            meanlines = meanlines,
            meanlabels = meanlabels,
            normalization = normalization,
        )
    end
end

function calculate_bins(charvec, binwidth)
    filteredcharvec = filter(
        !isempty, charvec
    )
    minbinvec = minimum.(filteredcharvec)
    maxbinvec = maximum.(filteredcharvec)
    minbin = minimum(minbinvec)
    maxbin = maximum(maxbinvec)
    minbin -= 3 * binwidth / 2
    maxbin += 3 * binwidth / 3
    if minbin == maxbin
        minbin -= binwidth
        maxbin += binwidth
    end
    return minbin:binwidth:maxbin
end

function compare_optimal_thresholds_chars_plot(
    optimal_thresholds_vec,
    plottingchars;
    plotdirpath = plotsdir("ensemble/optimal-thresholds"),
    plotformat = "png",
    force = false,
    kwargs...,
)
    unique_percent_clinic_tested = unique(
        optimal_thresholds_vec.percent_clinic_tested
    )

    mkpath(plotdirpath)

    for percent_clinic_tested in unique_percent_clinic_tested
        plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_best-thresholds"
        plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

        if !isfile(plotpath) || force
            optimal_thresholds_chars = filter(
                optimal_thresholds ->
                    optimal_thresholds.percent_clinic_tested ==
                    percent_clinic_tested ||
                        (
                            optimal_thresholds.percent_clinic_tested ==
                            1.0 &&
                            optimal_thresholds.individual_test_specification ==
                            CLINICAL_CASE_TEST_SPEC
                        ),
                optimal_thresholds_vec,
            )

            plot = create_optimal_thresholds_chars_plot(
                optimal_thresholds_chars,
                plottingchars;
                kwargs...,
            )

            save(
                joinpath(
                    plotdirpath,
                    "compare-outbreak_clinic-tested-$(percent_clinic_tested)_best-thresholds.png",
                ),
                plot;
                size = (2200, 1600),
            )

            Makie.empty!(plot)

            @info "Created optimal thresholds plot for % clinic tested $(percent_clinic_tested)"

            continue
        end

        GC.gc(true)

        @info "Optimal thresholds plot for % clinic tested $(percent_clinic_tested) already exists. Skipping."
    end

    return nothing
end

function create_optimal_thresholds_chars_plot(
    optimal_thresholds_chars,
    plottingchars;
    kwargs...,
)
    number_tests = length(optimal_thresholds_chars)
    number_plotting_chars = length(plottingchars)
    midpoint_plotting_chars = Int64(
        round(number_plotting_chars / 2; digits = 0)
    )

    sort!(
        optimal_thresholds_chars;
        by = threshold ->
            threshold.individual_test_specification.specificity,
    )

    fig = Figure()

    Label(
        fig[
            1,
            (midpoint_plotting_chars + 1):(midpoint_plotting_chars + 2),
            Top(),
        ],
        "Perc Clinic Tested: $(optimal_thresholds_chars[end].percent_clinic_tested), Noise: $(getdirpath(optimal_thresholds_chars[end].noise_specification))",
    )

    thresholdschars_structarr =
        optimal_thresholds_chars.outbreak_threshold_chars

    for (column, chartuple) in pairs(plottingchars)
        x = column + 1

        thresholdschars_vec =
            getproperty.(thresholdschars_structarr, chartuple.char)

        charvecs = reduce(vcat, thresholdschars_vec)

        if !haskey(chartuple, :bins)
            if !haskey(chartuple, :binwidth)
                @error "The metric $(chartuple.char) wasn't provided with bins or a binwidth"
                break
            end
            bins = calculate_bins(charvecs, chartuple.binwidth)
        else
            bins .= chartuple.bins
        end

        if !haskey(chartuple, :label)
            label = :none
        else
            label = chartuple.label
        end

        for (row, optimal_thresholds) in pairs(optimal_thresholds_chars)
            y = row + 1

            gl = fig[y, x] = GridLayout()
            ax = Axis(gl[1, 1]; xlabel = label)

            hist!(
                ax,
                reduce(vcat, thresholdschars_vec[row]);
                bins = bins,
                color = chartuple.color,
            )

            Label(
                fig[y, 1, Left()],
                "Sens: $(optimal_thresholds.individual_test_specification.sensitivity), Spec: $(optimal_thresholds.individual_test_specification.specificity),\nLag: $(optimal_thresholds.individual_test_specification.test_result_lag), Threshold: $(optimal_thresholds.alert_threshold)";
                rotation = pi / 2,
            )

            if row < number_tests
                hidexdecorations!(ax; ticklabels = false, ticks = false)
            end
        end
    end

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, 7)

    return fig
end

function compare_optimal_thresholds_test_chars_plot(
    optimal_thresholds_vec,
    plottingchars;
    plotdirpath = plotsdir("ensemble/optimal-thresholds"),
    plotformat = "png",
    force = false,
    kwargs...,
)
    unique_tests = unique(
        optimal_thresholds_vec.individual_test_specification
    )

    filter!(x -> !in(x, CLINICAL_TEST_SPECS), unique_tests)

    mkpath(plotdirpath)

    for test_specification in unique_tests
        plotname = "compare-outbreak_clinic-test-specification_sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)_best-thresholds"
        plotpath = joinpath(plotdirpath, "$plotname.$plotformat")

        if !isfile(plotpath) || force
            optimal_thresholds_chars = filter(
                optimal_thresholds ->
                    optimal_thresholds.individual_test_specification ==
                    test_specification,
                optimal_thresholds_vec,
            )

            plot = create_optimal_thresholds_test_chars_plot(
                optimal_thresholds_chars,
                plottingchars;
                kwargs...,
            )

            save(
                joinpath(
                    plotdirpath,
                    "compare-outbreak_clinic-test-specification_sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)_best-thresholds.png",
                ),
                plot;
                size = (2200, 1600),
            )

            Makie.empty!(plot)

            @info "Created optimal thresholds plot for test specification $(test_specification.sensitivity)_$(test_specification.specificity)_$(test_specification.test_result_lag)"
            continue
        end

        GC.gc(true)

        @info "Optimal thresholds plot for test specification $(test_specification.sensitivity)_$(test_specification.specificity)_$(test_specification.test_result_lag) already exists. Skipping."
    end

    return nothing
end

function create_optimal_thresholds_test_chars_plot(
    optimal_thresholds_chars,
    plottingchars;
    kwargs...,
)
    number_clinic_testing_rates = length(optimal_thresholds_chars)
    number_plotting_chars = length(plottingchars)
    midpoint_plotting_chars = Int64(
        round(number_plotting_chars / 2; digits = 0)
    )

    sort!(
        optimal_thresholds_chars;
        by = threshold ->
            threshold.percent_clinic_tested,
    )

    fig = Figure()

    Label(
        fig[
            1,
            (midpoint_plotting_chars + 1):(midpoint_plotting_chars + 2),
            Top(),
        ],
        "Test characteristics: Sensitivity $(optimal_thresholds_chars[end].individual_test_specification.sensitivity), Specificity $(optimal_thresholds_chars[end].individual_test_specification.specificity), Lag $(optimal_thresholds_chars[end].individual_test_specification.test_result_lag), Noise: $(getdirpath(optimal_thresholds_chars[end].noise_specification))",
    )

    thresholdschars_structarr =
        optimal_thresholds_chars.outbreak_threshold_chars

    for (column, chartuple) in pairs(plottingchars)
        x = column + 1

        thresholdschars_vec =
            getproperty.(thresholdschars_structarr, chartuple.char)

        charvecs = reduce(vcat, thresholdschars_vec)

        if !haskey(chartuple, :bins)
            if !haskey(chartuple, :binwidth)
                @error "The metric $(chartuple.char) wasn't provided with bins or a binwidth"
                break
            end
            bins = calculate_bins(charvecs, chartuple.binwidth)
        else
            bins .= chartuple.bins
        end

        if !haskey(chartuple, :label)
            label = :none
        else
            label = chartuple.label
        end

        for (row, optimal_thresholds) in pairs(optimal_thresholds_chars)
            y = row + 1

            gl = fig[y, x] = GridLayout()
            ax = Axis(gl[1, 1]; xlabel = label)

            hist!(
                ax,
                reduce(vcat, thresholdschars_vec[row]);
                bins = bins,
                color = chartuple.color,
            )

            Label(
                fig[y, 1, Left()],
                "% Clinic Tested: $(optimal_thresholds.percent_clinic_tested), Threshold: $(optimal_thresholds.alert_threshold)";
                rotation = pi / 2,
            )

            if row < number_clinic_testing_rates
                hidexdecorations!(ax; ticklabels = false, ticks = false)
            end
        end
    end

    rowsize!(fig.layout, 1, 5)
    colsize!(fig.layout, 1, 7)

    return fig
end

function Reff_ews_plot(
    incidencearr,
    Reffarr,
    Reff_thresholds_vec,
    ewsmetric_sa::T,
    ewsmetric::S,
    thresholdsarr,
    timeparams;
    sim = 1,
    outbreak_colormap = [
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    Reff_colormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    threshold = 5,
    aggregation = 1,
    metric_color = Makie.wong_colors()[1],
    kwargs...,
) where {T<:StructArray,S<:Symbol}
    @unpack trange = timeparams
    times = collect(trange) ./ 365
    aggregated_vec_length = size(Reffarr, 1) ÷ aggregation * aggregation
    times = times[1:aggregation:aggregated_vec_length]

    kwargs_dict = Dict(kwargs)

    filtered_thresholdsarr = thresholdsarr[sim][
        (thresholdsarr[sim][:, 4] .== 1), :,
    ]

    outbreak_status_vec = zeros(Int64, length(times))
    for (lower, upper, periodsum, outbreakstatus) in
        eachrow(filtered_thresholdsarr)
        outbreak_status_vec[(lower ÷ aggregation):1:(upper ÷ aggregation)] .=
            outbreakstatus
    end

    Reff_above_one = zeros(Int64, length(times))
    for (lower, upper, _) in
        eachrow(Reff_thresholds_vec[sim])
        Reff_above_one[(lower ÷ aggregation):1:(upper ÷ aggregation)] .= 1
    end

    ewsmetric_vec = getproperty(ewsmetric_sa, ewsmetric)[sim]

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
        color = Reff_above_one,
        colormap = Reff_colormap,
    )

    line_and_hline!(
        incax,
        times,
        aggregate_timeseries(@view(incidencearr[:, 1, sim]), aggregation),
        threshold * aggregation;
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

    # rowsize!(fig.layout, 0, Relative(0.05))

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
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    Reff_colormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
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
        N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR
    ],
    Reff_colormap = [
        N_MISSED_OUTBREAKS_COLOR, N_ALERTS_COLOR
    ],
    threshold = 5,
    metric_color = Makie.wong_colors()[1],
    spaero_metric_color = Makie.wong_colors()[3],
    linewidth = 3,
    hlinewidth = 2,
    aggregation = 1,
    kwargs...,
) where {T1<:StructArray,T2<:DataFrames.DataFrame,S<:Symbol}
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

function line_and_hline!(
    ax,
    times,
    yvec,
    hlinevec;
    linewidth = 3,
    hlinecolor = :black,
    hlinewidth = 2,
    hlinestyle = :dash,
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    color_keys = [:color, :colormap]

    if !(
        sum(k -> haskey(kwargs_dict, k), color_keys) in (0, length(color_keys))
    )
        error(
            "Must provide both $color_keys for the `lines!` function, or none. You have specified: $(keys(kwargs_dict))"
        )
    end

    lines!(
        ax,
        times,
        yvec;
        linewidth = linewidth,
        kwargs...,
    )
    hlines!(
        ax,
        hlinevec;
        color = hlinecolor,
        linewidth = hlinewidth,
        linestyle = hlinestyle,
    )
    return nothing
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
        (N_MISSED_OUTBREAKS_COLOR, outbreak_alpha),
        (PERC_OUTBREAKS_DETECTED_COLOR, outbreak_alpha),
    ],
    Reff_colormap = [
        (N_MISSED_OUTBREAKS_COLOR, Reff_alpha),
        (N_ALERTS_COLOR, Reff_alpha),
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

    hlines!(reffax, 1; linestyle = :dash)
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

function tycho_epicurve(
    plot_dates,
    weekly_cases::T1,
    biweekly_cases::T1,
    monthly_cases::T1;
    plottitle = "",
    subtitle = "",
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
) where {T1<:AbstractVector{<:Integer}}
    xticks = calculate_xticks(plot_dates)

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        subtitle = subtitle,
        xlabel = "Date",
        ylabel = "Number of Cases",
        xticks = xticks,
    )

    tycho_epicurve!(
        ax,
        plot_dates,
        weekly_cases,
        biweekly_cases,
        monthly_cases;
        obsdate = obsdate,
    )

    Legend(fig[1, 2], ax, "Aggregation (wks)")

    return fig
end

function tycho_epicurve(
    plot_dates,
    cases_tuple::T1,
    ews_tuple::T2;
    plottitle = "",
    subtitle = "",
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
    ews_ylabel = "EWS",
    cases_ylabel = "Test Positives",
) where {T1<:Tuple,T2<:Tuple}
    xticks = calculate_xticks(plot_dates)

    weekly_cases, biweekly_cases, monthly_cases = cases_tuple
    weekly_ews, biweekly_ews, monthly_ews = ews_tuple

    fig = Figure()
    ews_ax = Axis(
        fig[1, 1];
        xlabel = "Date",
        ylabel = ews_ylabel,
    )
    cases_ax = Axis(
        fig[2, 1];
        xlabel = "Date",
        ylabel = cases_ylabel,
        xticks = xticks,
    )

    tycho_epicurve!(
        cases_ax,
        plot_dates,
        weekly_cases,
        biweekly_cases,
        monthly_cases;
        obsdate = obsdate,
    )

    ews_timeseries!(
        ews_ax,
        plot_dates,
        weekly_ews,
        biweekly_ews,
        monthly_ews;
        obsdate = obsdate,
    )

    linkxaxes!(cases_ax, ews_ax)
    hidexdecorations!(ews_ax)

    Legend(fig[:, 2], cases_ax, "Aggregation (wks)")

    Label(fig[1, 1, Top()], plottitle * "\n" * subtitle; font = :bold)

    return fig
end

function ews_timeseries!(
    ax,
    plot_dates,
    weekly_ews,
    biweekly_ews,
    monthly_ews;
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
)
    weekly_xaxes = (1:length(weekly_ews)) .* 7
    biweekly_xaxes = (1:length(biweekly_ews)) .* (2 * 7) .- 7
    monthly_xaxes = (1:length(monthly_ews)) .* (4 * 7) .- 21

    obsdate_indices = findfirst(x -> x == obsdate, plot_dates)
    vspan!(
        ax,
        obsdate_indices,
        length(plot_dates);
        color = (:black, 0.1),
    )

    lines!(
        ax, monthly_xaxes, monthly_ews; color = :darkred, label = "Monthly"
    )
    lines!(
        ax, biweekly_xaxes, biweekly_ews; color = :blue, label = "Biweekly"
    )
    lines!(ax, weekly_xaxes, weekly_ews; color = :black, label = "Weekly")

    return nothing
end

function tycho_epicurve(
    plot_dates,
    cases_tuple::T1,
    ews_tuple::T2,
    ews_threshold_tuple::T3,
    ews_threshold_indices::T4;
    plottitle = "",
    subtitle = "",
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
    ews_ylabel = "EWS",
    cases_ylabel = "Test Positives",
    threshold_percentile = 0.95,
    consecutive_thresholds = 2,
) where {
    S1<:AbstractVector,
    S2<:AbstractMatrix,
    S3<:Union{Nothing,Integer},
    T1<:Tuple{<:S1,<:S1,<:S1},
    T2<:Tuple{<:S1,<:S1,<:S1},
    T3<:Tuple{<:S2,<:S2,<:S2},
    T4<:Tuple{<:S3,<:S3,<:S3},
}
    reshaped_ews_threshold_tuple = map(x -> reshape(x, :), ews_threshold_tuple)

    for i in eachindex(reshaped_ews_threshold_tuple)
        @assert length(reshaped_ews_threshold_tuple[i]) ==
            length(ews_threshold_tuple[i][:, 1])
    end

    return tycho_epicurve(
        plot_dates,
        cases_tuple,
        ews_tuple,
        reshaped_ews_threshold_tuple,
        ews_threshold_indices;
        plottitle = plottitle,
        subtitle = subtitle,
        obsdate = obsdate,
        ews_ylabel = ews_ylabel,
        cases_ylabel = cases_ylabel,
        threshold_percentile = threshold_percentile,
        consecutive_thresholds = consecutive_thresholds,
    )
end

function tycho_epicurve(
    plot_dates,
    cases_tuple::T1,
    ews_tuple::T2,
    ews_threshold_tuple::T3,
    ews_threshold_indices::T4;
    plottitle = "",
    subtitle = "",
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
    ews_ylabel = "EWS",
    cases_ylabel = "Test Positives",
    threshold_percentile = 0.95,
    consecutive_thresholds = 2,
) where {
    S1<:AbstractVector,
    S2<:Union{Nothing,Integer},
    T1<:Tuple{<:S1,<:S1,<:S1},
    T2<:Tuple{<:S1,<:S1,<:S1},
    T3<:Tuple{<:S1,<:S1,<:S1},
    T4<:Tuple{<:S2,<:S2,<:S2},
}
    xticks = calculate_xticks(plot_dates)

    weekly_cases, biweekly_cases, monthly_cases = cases_tuple

    fig = Figure()
    ews_ax = Axis(
        fig[1, 1];
        xlabel = "Date",
        ylabel = ews_ylabel,
    )
    cases_ax = Axis(
        fig[2, 1];
        xlabel = "Date",
        ylabel = cases_ylabel,
        xticks = xticks,
    )

    tycho_epicurve!(
        cases_ax,
        plot_dates,
        weekly_cases,
        biweekly_cases,
        monthly_cases;
        obsdate = obsdate,
    )

    ews_timeseries!(
        ews_ax,
        plot_dates,
        ews_tuple,
        ews_threshold_tuple,
        ews_threshold_indices;
        obsdate = obsdate,
    )

    linkxaxes!(cases_ax, ews_ax)
    hidexdecorations!(ews_ax)

    Legend(fig[:, 2], cases_ax, "Aggregation (wks)")

    Label(
        fig[1, 2],
        "Circles Represent\nExceeds Threshold\nof $threshold_percentile Percentile\n\nVertical Lines Represent\n$consecutive_thresholds Consecutive EWS\nExceeds Threshold",
        ;
        valign = :center,
        font = :bold,
    )
    Label(fig[1, 1, Top()], plottitle * "\n" * subtitle; font = :bold)

    rowsize!(fig.layout, 1, Relative(0.5))

    return fig
end

function ews_timeseries!(
    ax,
    plot_dates,
    ews_tuple,
    ews_threshold_tuple,
    ews_threshold_indices;
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
)
    obsdate_indices = findfirst(x -> x == obsdate, plot_dates)
    vspan!(
        ax,
        obsdate_indices,
        length(plot_dates);
        color = (:black, 0.1),
    )

    for (ews, threshold, ind, aggregation, color, label) in zip(
        ews_tuple,
        ews_threshold_tuple,
        ews_threshold_indices,
        [1, 2, 4],
        [:grey20, :blue, :darkred],
        ["Weekly", "Biweekly", "Monthly"],
    )
        multiplier = aggregation * 7

        xaxes = (1:length(ews)) .* multiplier

        threshold = ews .* threshold
        replace!(
            threshold,
            0.0 => NaN,
            -0.0 => NaN,
        )

        lines!(
            ax, xaxes, ews; color = color, label = label
        )
        scatter!(
            ax,
            xaxes,
            threshold;
            markersize = 10,
            strokecolor = color,
            strokewidth = 2,
            color = (color, 0.4),
        )
        if !isnothing(ind)
            ind = ind .* multiplier
            vlines!(ax, ind; color = color)
        end
    end

    return nothing
end

function calculate_xticks(plot_dates)
    unique_years = unique(Dates.year.(plot_dates))
    year_indices = map(
        x -> findfirst(==(x), Dates.year.(plot_dates)), unique_years
    )

    return (year_indices, string.(unique_years))
end

function tycho_epicurve!(
    ax,
    plot_dates,
    weekly_cases,
    biweekly_cases,
    monthly_cases;
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
)
    xaxes = 1:length(plot_dates)

    obsdate_indices = findfirst(x -> x == obsdate, plot_dates)
    vspan!(ax, obsdate_indices, maximum(xaxes); color = (:black, 0.1))

    lines!(ax, xaxes, monthly_cases; color = :darkred, label = "Monthly")
    band!(
        ax, xaxes, 0, biweekly_cases; color = (:blue, 0.2), label = "Biweekly"
    )
    lines!(ax, xaxes, biweekly_cases; color = :blue)
    band!(ax, xaxes, 0, weekly_cases; color = :black, label = "Weekly")

    return nothing
end

function tycho_noise_components_epicurve(
    plot_dates,
    cases_tuple,
    noise_tuple;
    plottitle = "",
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
)
    weekly_cases, biweekly_cases, monthly_cases = cases_tuple
    weekly_noise, biweekly_noise, monthly_noise = noise_tuple

    xticks = calculate_xticks(plot_dates)

    fig = Figure()
    gl = fig[1, 1] = GridLayout()

    obsax = Axis(gl[1, 1]; title = "Observed", xlabel = "Date")
    tycho_epicurve!(
        obsax,
        plot_dates,
        weekly_cases .+ weekly_noise,
        biweekly_cases .+ biweekly_noise,
        monthly_cases .+ monthly_noise;
        obsdate = obsdate,
    )

    noiseax = Axis(gl[2, 1]; title = "Noise", xlabel = "Date")
    tycho_epicurve!(
        noiseax,
        plot_dates,
        weekly_noise,
        biweekly_noise,
        monthly_noise;
        obsdate = obsdate,
    )

    incax = Axis(
        gl[3, 1];
        title = "Incidence",
        xlabel = "Date",
        xticks = xticks,
    )

    tycho_epicurve!(
        incax,
        plot_dates,
        weekly_cases,
        biweekly_cases,
        monthly_cases;
        obsdate = obsdate,
    )

    map(ax -> hidexdecorations!(ax), (obsax, noiseax))
    linkxaxes!(incax, noiseax, obsax)

    Legend(fig[1, 2], incax, "Aggregation (wks)")

    Label(fig[0, :], plottitle)
    rowsize!(fig.layout, 0, Relative(0.03))

    return fig
end

function tycho_test_positive_components_epicurve(
    plot_dates,
    inc_vec_tuple,
    test_arr_tuple;
    plottitle = "",
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
)
    weekly_cases, biweekly_cases, monthly_cases = inc_vec_tuple
    test_weekly_cases_arr, test_biweekly_cases_arr, test_monthly_cases_arr =
        test_arr_tuple

    xticks = calculate_xticks(plot_dates)

    fig = Figure()
    gl = fig[1, 1] = GridLayout()

    trueposax = Axis(gl[1, 1]; title = "True Positives", xlabel = "Date")
    tycho_epicurve!(
        trueposax,
        plot_dates,
        test_weekly_cases_arr[:, 3],
        test_biweekly_cases_arr[:, 3],
        test_monthly_cases_arr[:, 3];
    )

    falseposax = Axis(gl[2, 1]; title = "False Positives", xlabel = "Date")
    tycho_epicurve!(
        falseposax,
        plot_dates,
        test_weekly_cases_arr[:, 4],
        test_biweekly_cases_arr[:, 4],
        test_monthly_cases_arr[:, 4];
    )

    testposax = Axis(gl[3, 1]; title = "Test Positive", xlabel = "Date")
    tycho_epicurve!(
        testposax,
        plot_dates,
        test_weekly_cases_arr[:, 5],
        test_biweekly_cases_arr[:, 5],
        test_monthly_cases_arr[:, 5];
    )

    ntestedax = Axis(gl[4, 1]; title = "Total Tested", xlabel = "Date")
    tycho_epicurve!(
        ntestedax,
        plot_dates,
        test_weekly_cases_arr[:, 1] .+ test_weekly_cases_arr[:, 2],
        test_biweekly_cases_arr[:, 1] .+ test_biweekly_cases_arr[:, 2],
        test_monthly_cases_arr[:, 1] .+ test_monthly_cases_arr[:, 2];
    )

    incax = Axis(
        gl[5, 1];
        title = "Incidence",
        xlabel = "Date",
        xticks = xticks,
    )

    tycho_epicurve!(
        incax,
        plot_dates,
        weekly_cases,
        biweekly_cases,
        monthly_cases;
        obsdate = obsdate,
    )

    map(
        ax -> hidexdecorations!(ax),
        (trueposax, falseposax, testposax, ntestedax),
    )
    linkxaxes!(incax, trueposax, falseposax, testposax, ntestedax)

    Legend(fig[1, 2], incax, "Aggregation (wks)")

    Label(fig[0, :], plottitle)
    rowsize!(fig.layout, 0, Relative(0.03))

    return fig
end

function tycho_tau_distribution(
    tested_ews_tuple,
    actual_ews_sa,
    tau_metric;
    plottitle = string(tau_metric) * " Distribution",
)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Kendall's Tau",
        ylabel = "Count",
    )
    for (i, (aggregation, ews)) in
        zip(1:length(tested_ews_tuple), pairs(tested_ews_tuple))
        aggregation_tau = getproperty(ews, tau_metric)
        if sum(isnan.(aggregation_tau)) == length(aggregation_tau)
            continue
        end
        hist!(
            ax,
            aggregation_tau;
            color = (Makie.wong_colors()[i], 0.5),
            bins = 20,
            label = string(aggregation),
        )
        vlines!(
            ax,
            mean(aggregation_tau);
            color = Makie.wong_colors()[i],
            linestyle = :dash,
            linewidth = 4,
        )
    end

    vlines!(
        ax,
        mean(getproperty(actual_ews_sa, tau_metric));
        color = :grey20,
        linewidth = 4,
    )

    text!(
        ax,
        mean(getproperty(actual_ews_sa, tau_metric)) + 0.005,
        length(getproperty(tested_ews_tuple[1], tau_metric)) / 10;
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

function tycho_tau_heatmap_plot(
    df;
    baseline_test = IndividualTestSpecification(1.0, 1.0, 0),
    colormap = :RdBu,
    textcolorthreshold = 0.6,
    statistic_function = "mean",
    plottitle = "Kendall's Tau Heatmap: " * titlecase(statistic_function),
)
    df[!, :test_sens] = getproperty.(df.test_specification, :sensitivity)
    df[!, :test_spec] = getproperty.(df.test_specification, :sensitivity)
    df[!, :test_result_lag] =
        getproperty.(df.test_specification, :test_result_lag)
    sort!(
        df,
        [:test_sens, :test_spec, :test_result_lag];
        rev = [true, true, false],
    )

    unique_tests = unique(df.test_specification)
    if mapreduce(x -> x == baseline_test, +, unique_tests) == 0
        error("default_test $baseline_test must be in unique_tests")
    end

    ordered_df =
        unstack(
            select(df, [:ews_metric, :test_specification, :ews_metric_value]),
            :test_specification,
            :ews_metric_value,
        ) |>
        df -> sort(df, order(2; rev = false))

    default_test_metric_order = ordered_df.ews_metric

    mat = Matrix(ordered_df[:, 2:end])'

    function test_axis_label(test)
        return "($(test.sensitivity), $(test.specificity), $(test.test_result_lag))"
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Test Specification (Sensitivity, Specificity, Lag)",
        ylabel = "EWS Metric",
        xticks = (1:length(unique_tests), test_axis_label.(unique_tests)),
        yticks = (
            1:length(default_test_metric_order), default_test_metric_order
        ),
    )

    hmap = heatmap!(
        ax,
        mat;
        colormap = colormap,
        colorrange = (-1, 1),
    )

    for j in axes(mat, 2), i in axes(mat, 1)
        val = mat[i, j]
        textcolor = abs(val) < textcolorthreshold ? :black : :white
        text!(
            ax,
            "$(round(mat[i,j], digits = 5))";
            position = (i, j),
            color = textcolor,
            align = (:center, :center),
        )
    end

    limits!(
        ax,
        (0, length(unique_tests) + 1),
        (0, length(default_test_metric_order) + 1),
    )
    #
    Colorbar(fig[1, 2], hmap; label = "values", width = 15, ticksize = 15)
    return fig
end

function ews_lead_time_plot(
    lead_time_df;
    week_aggregation = 1,
    ews_method = Main.Backward,
    lead_time_units = :days,
)
    unique_tests = unique(lead_time_df.test_specification)
    lead_time_units_string = string(lead_time_units)

    unique_noise = unique(lead_time_df.noise_magnitude)

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Lead Time By Test and Noise Magnitude",
        subtitle = "$week_aggregation Week Aggregation and $(titlecase(method_string(ews_method))) EWS",
        xlabel = "Noise Magnitude",
        ylabel = "Lead Time ($(lead_time_units))",
        xticks = minimum(unique_noise):maximum(unique_noise),
    )

    for (i, test) in pairs(unique_tests)
        filtered_df = subset(
            lead_time_df,
            :test_specification => ByRow(==(test)),
            :week_aggregation => ByRow(==(week_aggregation)),
            :ews_method => ByRow(==(ews_method)),
            :lead_time_units => ByRow(==(lead_time_units_string)),
        )

        band!(
            ax,
            filtered_df.noise_magnitude,
            filtered_df.lead_time_lower,
            filtered_df.lead_time_upper;
            color = (Makie.wong_colors()[i], 0.1),
            # label = get_test_description(test),
        )

        scatterlines!(
            ax,
            filtered_df.noise_magnitude,
            filtered_df.lead_time_median;
            color = Makie.wong_colors()[i],
            label = get_test_description(test),
        )
    end

    Legend(fig[1, 2], ax, "Test")

    return fig
end
