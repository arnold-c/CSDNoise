export ews_reff_histogram_plot

function ews_reff_histogram_plot(
        survival_df::T1;
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
        nbanks = 1,
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
        hist_ax = nothing
        test_specification_vec = unique(noise_gdf.test_specification)
        for (i, test_specification) in pairs(test_specification_vec)
            detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
                DF.subset(
                    noise_gdf,
                    :test_specification => DF.ByRow(==(test_specification)),
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
                # max_count = maximum(enddate_counts)
                # ytick_step = max_count / 2
                hist_ax = Axis(
                    gl[2, 1];
                    limits = (
                        5, maximum(enddate_times),
                        nothing, nothing,
                        # nothing, ceil(1.1 * max_count),
                    ),
                    xticks = (1:1:10, string.(collect(-4:1:5))),
                    # yticks = (
                    #     0:ytick_step:max_count, string.(collect(0:50:100))
                    # ),
                    xticklabelsize = xticklabelsize,
                    yticklabelsize = yticklabelsize,
                )

                # hideydecorations!(hist_ax)

                barplot!(
                    hist_ax,
                    enddate_times,
                    enddate_counts;
                    color = tipping_point_color,
                )

                # survival_plot_histogram!(
                #     gl[2, 1],
                #     enddate_times,
                #     enddate_counts;
                #     tipping_point_color = tipping_point_color,
                #     trim_burnin = trim_burnin,
                #     burnin = burnin,
                # )
                # surv_ax = Axis(
                #     gl[2, 1];
                #     xticks = (
                #         0:1:10,
                #         mapreduce(
                #             i -> if i < 5
                #                 [string(i), ""]
                #             else
                #                 [string(i)]
                #             end,
                #             vcat, 0:5,
                #         ),
                #     ),
                #     xticklabelsize = xticklabelsize,
                #     yticklabelsize = yticklabelsize,
                #     limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
                # )
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
        end
        if ax_position[2] != 1
            hideydecorations!(hist_ax)
        end

        nrows = ceil(num_noise / 2)
        if ax_position[1] != nrows && ax_position[2] + 2 <= num_noise
            hidexdecorations!(hist_ax)
        end
    end

    Legend(
        fig[0, :],
        [
            PolyElement(; color = tipping_point_color),
        ],
        ["Tipping Point"],
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
