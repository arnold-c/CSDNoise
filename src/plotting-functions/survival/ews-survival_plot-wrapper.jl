export ews_survival_plot

function ews_survival_plot(
        survival_df::T1;
        test_specification_vec = [
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
        ],
        noise_specification_vec = [
            NoiseSpecification(PoissonNoise(1.0)),
            NoiseSpecification(PoissonNoise(7.0)),
            NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.8734)),
            NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.102)),
        ],
        endpoint_aggregation = Dates.Day(30),
        linestyle_vec = [:solid, :dot],
        emergent_color = "#C2192D",
        null_color = "#453948",
        tipping_point_color = ("#A4BED5", 0.7),
        alpha = 1.0,
        trim_burnin = true,
        plottitle = "Survival",
        subtitle = "",
        nbanks = length(test_specification_vec) + length(linestyle_vec) + 1,
        xlabel = "Time After Burn-In (Years)",
        ylabel = "Survival Numbers",
        legend_rowsize = GLMakie.Relative(0.05),
        xlabel_rowsize = GLMakie.Relative(0.03),
        ylabel_rowsize = GLMakie.Relative(0.02),
        facet_fontsize = 20,
        legendsize = 22,
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 22,
        yticklabelsize = 22,
        kwargs...,
    ) where {T1 <: DF.DataFrame}
    kwargs_dict = Dict{Symbol, Any}(kwargs)
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

    subset_survival_df = DF.subset(
        survival_df,
        :test_specification => DF.ByRow(in(test_specification_vec)),
        :noise_specification => DF.ByRow(in(noise_specification_vec)),
    )

    noise_grouped_dfs = DF.groupby(subset_survival_df, :noise_specification)

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

        ax_position = if noise_num == 1
            (1, 1)
        elseif noise_num == 2
            (1, 2)
        elseif noise_num == 3
            (2, 1)
        elseif noise_num == 4
            (2, 2)
        else
            @error "Too many noise types"
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
                    gl[2, 1],
                    enddate_times,
                    enddate_counts;
                    tipping_point_color = tipping_point_color,
                    trim_burnin = trim_burnin,
                    burnin = burnin,
                )
                surv_ax = Axis(
                    gl[2, 1];
                    xticks = (1:1:10, string.(collect(-4:1:5))),
                    xticklabelsize = xticklabelsize,
                    yticklabelsize = yticklabelsize,
                    limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
                )
                Box(gl[1, 1]; color = :lightgray, strokevisible = false)
                Label(
                    gl[1, 1],
                    noise_description;
                    fontsize = facet_fontsize,
                    padding = (0, 0, 0, 0),
                    valign = :bottom,
                    tellwidth = false,
                )
                rowsize!(gl, 2, Relative(0.9))
            end

            survival_plot_lines!(
                surv_ax,
                times,
                detection_survival_times,
                detection_survival_vec,
                null_survival_times,
                null_survival_vec;
                emergent_color = emergent_color,
                null_color = null_color,
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
                    col in [emergent_color, null_color, tipping_point_color]
            ],
        ),
        vcat(
            get_test_description.(test_specification_vec),
            ["Emergent", "Null", "Tipping Point"],
        ),
        "";
        nbanks = nbanks,
        labelsize = legendsize,
    )

    rowsize!(fig.layout, 0, legend_rowsize)
    if num_noise == 1
        colsize!(fig.layout, 1, GLMakie.Relative(1))
    end

    Label(
        fig[:, 0],
        ylabel;
        fontsize = ylabelsize,
        font = :bold,
        rotation = pi / 2,
    )
    colsize!(fig.layout, 0, ylabel_rowsize)

    xlabel_position = num_noise > 2 ? 3 : 2
    Label(
        fig[xlabel_position, :],
        xlabel;
        fontsize = xlabelsize,
        font = :bold,
    )
    rowsize!(fig.layout, xlabel_position, xlabel_rowsize)

    return fig
end

function ews_survival_plot(
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        facet_title = "Survival",
        ews_aggregation = Dates.Day(7),
        burnin = Year(5),
        endpoint_aggregation = Dates.Day(30),
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
