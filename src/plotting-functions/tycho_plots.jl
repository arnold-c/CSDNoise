using GLMakie
using DataFrames
using Dates: Dates

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
