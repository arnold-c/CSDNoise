function ews_survival_plot(
    survival_df::T1;
    test_specification_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
    ],
    noise_specification_vec = [
        PoissonNoiseSpecification(1.0),
        PoissonNoiseSpecification(7.0),
        DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.8734),
        DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.102),
    ],
    endpoint_aggregation = Dates.Day(30),
    linestyle_vec = [:solid, :dot],
    alpha = 1.0,
    trim_burnin = true,
    plottitle = "Survival",
    subtitle = "",
    nbanks = length(test_specification_vec) * length(linestyle_vec),
    legend_rowsize = Makie.Relative(0.05),
) where {T1<:DataFrames.DataFrame}
    if sum(
        filter(
            !isnothing,
            indexin(
                unique(survival_df.test_specification),
                test_specification_vec,
            ),
        ),
    ) != sum(eachindex(test_specification_vec))
        error(
            "Not all test specified are present in the dataframe.\nTrying to plot survival curves for:\n$(test_specification_vec).\nFound these in the dataframe:\n$(unique(survival_df.test_specification))"
        )
    end

    if length(test_specification_vec) > length(linestyle_vec)
        error(
            "There are more tests specified ($(length(test_specification_vec))) than line styles ($(length(linestyle_vec)))."
        )
    end

    if sum(
        filter(
            !isnothing,
            indexin(
                unique(survival_df.noise_specification),
                noise_specification_vec,
            ),
        ),
    ) != sum(eachindex(noise_specification_vec))
        error(
            "Not all noise structures specified are present in the dataframe.\nTrying to plot survival curves for:\n$(noise_specification_vec).\nFound these in the dataframe:\n$(unique(survival_df.noise_specification))"
        )
    end

    @assert length(unique(survival_df.ews_metric_specification)) == 1
    @assert length(unique(survival_df.ews_threshold_burnin)) == 1
    ews_aggregation = survival_df.ews_metric_specification[1].aggregation
    burnin = survival_df.ews_threshold_burnin[1]

    subset_survival_df = subset(
        survival_df,
        :test_specification => ByRow(in(test_specification_vec)),
        :noise_specification => ByRow(in(noise_specification_vec)),
    )

    noise_grouped_dfs = groupby(subset_survival_df, :noise_specification)

    num_noise = length(noise_specification_vec)

    if num_noise > 4
        error(
            "Trying to plot too many metric facets $(length(ewsmetric_grouped_dfs)). Max allowed is 4"
        )
    end

    fig = Figure()

    for (noise_num, noise_gdf) in enumerate(noise_grouped_dfs)
        noise_description = noise_table_description(
            noise_gdf.noise_specification[1]
        )
        ax_position = Match.@match noise_num begin
            1 => (1, 1)
            2 => (1, 2)
            3 => (2, 1)
            4 => (2, 2)
        end
        gl = fig[ax_position...] = GridLayout()
        surv_ax = nothing
        for (i, test_specification) in pairs(test_specification_vec)
            detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
                subset(
                    noise_gdf,
                    :test_specification => ByRow(==(test_specification)),
                ),
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
                noise_gdf.enddate;
                ews_aggregation = ews_aggregation,
                endpoint_aggregation = endpoint_aggregation,
            )

            if i == 1
                survival_plot_histogram!(
                    gl,
                    enddate_times,
                    enddate_counts;
                    trim_burnin = trim_burnin,
                    burnin = burnin,
                )
                surv_ax = Axis(
                    gl[1, 1];
                    title = noise_description,
                    xlabel = "Time (Years)",
                    ylabel = "Survival Numbers",
                    xticks = 1:1:10,
                    limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
                )
            end

            survival_plot_lines!(
                surv_ax,
                times,
                detection_survival_times,
                detection_survival_vec,
                null_survival_times,
                null_survival_vec;
                alpha = alpha,
                nsims = nsims,
                trim_burnin = trim_burnin,
                burnin = burnin,
                linestyle = linestyle_vec[i],
            )
        end
        if ax_position[2] != 1
            hideydecorations!(surv_ax)
        end

        nrows = ceil(num_noise / 2)
        if ax_position[1] != nrows && ax_position[2] + 2 <= num_noise
            hidexdecorations!(surv_ax)
        end
    end

    Legend(
        fig[0, :],
        vcat(
            [
                LineElement(; linestyle = style) for
                style in linestyle_vec[1:length(test_specification_vec)]
            ],
            [
                PolyElement(; color = col) for
                col in [:red, :blue]
            ],
        ),
        vcat(
            get_test_description.(test_specification_vec),
            ["Null", "Emergent"],
        ),
        "";
        nbanks = nbanks,
    )

    rowsize!(fig.layout, 0, legend_rowsize)
    if num_noise == 1
        colsize!(fig.layout, 1, Makie.Relative(1))
    end

    return fig
end

function ews_survival_plot(
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    facet_title = "Survival",
    ews_aggregation = Day(7),
    burnin = Year(5),
    endpoint_aggregation = Day(30),
    alpha = 1.0,
    trim_burnin = true,
)
    fig = Figure()

    gl = fig[1, 1] = GridLayout()

    ews_survival_facet!(
        gl,
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        facet_title = facet_title,
        ews_aggregation = ews_aggregation,
        burnin = burnin,
        endpoint_aggregation = endpoint_aggregation,
        alpha = alpha,
        trim_burnin = trim_burnin,
    )

    Legend(
        fig[1, 2],
        contents(gl)[2];
        orientation = :vertical,
    )

    return fig
end

function ews_survival_facet!(
    gl,
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    facet_title = "Survival",
    ews_aggregation = Day(7),
    burnin = Year(5),
    endpoint_aggregation = Day(30),
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
        alpha = alpha,
        nsims = nsims,
        facet_title = facet_title,
        trim_burnin = trim_burnin,
        burnin = burnin,
    )

    return nothing
end

function prepare_survival_facet_params(
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    ews_aggregation = Day(7),
    endpoint_aggregation = Day(30),
)
    @unpack detection_survival_vec, detection_indices_vec =
        detection_survival_vecs
    @unpack null_survival_vec, null_indices_vec = null_survival_vecs

    filtered_enddate_vec = filter(isinteger, enddate_vec)

    times =
        collect(1:Dates.days(ews_aggregation):maximum(filtered_enddate_vec)) ./
        365

    detection_survival_times = vcat(
        0,
        times[detection_indices_vec]...,
        times[end],
    )
    null_survival_times = vcat(
        0,
        times[null_indices_vec]...,
        times[end],
    )

    nsims = detection_survival_vec[1]

    enddate_vec = div.(enddate_vec, Dates.days(endpoint_aggregation))
    unique_enddate_vec = sort(unique(enddate_vec))
    enddate_times =
        (unique_enddate_vec .* Dates.days(endpoint_aggregation)) ./ 365

    enddate_counts = map(unique_enddate_vec) do enddate
        sum(enddate_vec .== enddate)
    end

    return times, enddate_times,
    enddate_counts,
    detection_survival_times,
    detection_survival_vec,
    null_survival_times,
    null_survival_vec,
    nsims
end

function survival_plot_histogram!(
    gl,
    enddate_times,
    enddate_counts;
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
        color = (:darkgray, 1.0),
    )

    if trim_burnin
        xlims!(hist_ax, Dates.value(burnin), enddate_times[end])
    end

    return nothing
end

function survival_plot_lines!(
    gl::T1,
    times,
    detection_survival_times,
    detection_survival_vec,
    null_survival_times,
    null_survival_vec;
    linestyle = :solid,
    alpha = 1.0,
    nsims = 100,
    facet_title = "Survival",
    trim_burnin = true,
    burnin = Dates.Year(5),
) where {T1<:Makie.GridLayout}
    surv_ax = Axis(
        gl[1, 1];
        title = facet_title,
        xlabel = "Time (Years)",
        ylabel = "Survival Numbers",
        xticks = 1:1:10,
        limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
    )
    survival_plot_lines!(
        surv_ax,
        times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        linestyle = linestyle,
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
    alpha = 1.0,
    nsims = 100,
    facet_title = "Survival",
    trim_burnin = true,
    burnin = Dates.Year(5),
) where {T1<:Makie.Axis}
    lines!(
        surv_ax,
        detection_survival_times,
        detection_survival_vec;
        color = (:blue, alpha),
        linestyle = linestyle,
        label = "Detection",
    )

    lines!(
        surv_ax,
        null_survival_times,
        null_survival_vec;
        color = (:red, alpha),
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
