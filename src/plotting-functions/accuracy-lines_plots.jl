using DataFrames
using DrWatson: DrWatson
using StatsBase: StatsBase
using Match: Match
using StructArrays

function prepare_line_plot_df!(
    output_df,
    gdf::T1,
    ewsmetric = "mean",
    outcome = :accuracy,
    tests = [
        IndividualTestSpecification(0.8, 0.8, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.95, 0.95, 0),
        IndividualTestSpecification(0.97, 0.97, 0),
        IndividualTestSpecification(0.98, 0.98, 0),
        IndividualTestSpecification(0.99, 0.99, 0),
        IndividualTestSpecification(1.0, 1.0, 0),
    ];
    tiebreaker_preference = "specificity",
) where {T1<:Union{<:DataFrame,<:SubDataFrame}}
    tiebreaker_args = Match.@match tiebreaker_preference begin
        "speed" => (:ews_consecutive_thresholds, false)
        "specificity" => (:specificity, true)
        _ => error(
            "Invalid preference: $tiebreaker_preference. Please choose either \"speed\" or \"specificity\"."
        )
    end

    metric_df = subset(
        gdf,
        :ews_metric => ByRow(==(ewsmetric)),
        :test_specification => ByRow(in(tests)),
    )

    filtered_test_df =
        map(collect(groupby(metric_df, :test_specification))) do df
            sort(df, order(tiebreaker_args[1]; rev = tiebreaker_args[2]))[
                1, :,
            ]
        end |>
        x -> vcat(DataFrame.(x)...; cols = :union)

    filtered_test_df[!, :test_sens] =
        getproperty.(filtered_test_df.test_specification, :sensitivity)
    filtered_test_df[!, :test_spec] =
        getproperty.(filtered_test_df.test_specification, :sensitivity)
    filtered_test_df[!, :test_result_lag] =
        getproperty.(filtered_test_df.test_specification, :test_result_lag)
    sort!(
        filtered_test_df,
        [:test_sens, :test_spec, :test_result_lag];
        rev = [true, true, false],
    )

    append!(output_df, filtered_test_df; cols = :union)

    return nothing
end

function line_plot(
    df;
    legend_rowsize = Relative(0.05),
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    if !haskey(kwargs_dict, :ylims)
        ylims = (
            maximum([0, minimum(df.accuracy) - 0.05]),
            minimum([1, maximum(df.accuracy) + 0.05]),
        )
    else
        ylims = kwargs_dict[:ylims]
    end

    noise_descriptions =
        heatmap_noise_description.(unique(df.noise_specification))
    num_noise_descriptions = length(noise_descriptions)

    if !haskey(kwargs_dict, :nbanks)
        nbanks = num_noise_descriptions
    else
        nbanks = kwargs_dict[:nbanks]
    end

    ewsmetric_grouped_dfs = groupby(df, :ews_metric)

    num_metrics = length(ewsmetric_grouped_dfs)
    if num_metrics > 4
        error(
            "Trying to plot too many metric facets $(length(ewsmetric_grouped_dfs)). Max allowed is 4",
        )
    end

    fig = Figure()

    for (metric_num, metric_gdf) in enumerate(ewsmetric_grouped_dfs)
        ax_position = Match.@match metric_num begin
            1 => (1, 1)
            2 => (1, 2)
            3 => (2, 1)
            4 => (2, 2)
        end
        gl = fig[ax_position...] = GridLayout()
        line_plot_facet!(
            gl,
            metric_gdf;
            ax_position = ax_position,
            num_metrics = num_metrics,
            facet_title = sentencecase(replace(metric_gdf[1, :ews_metric], "_" => " ")),
            ylims = ylims,
        )
    end

    Legend(
        fig[0, :],
        [
            PolyElement(; color = col) for
            col in Makie.wong_colors()[1:num_noise_descriptions]
        ],
        noise_descriptions,
        "";
        nbanks = nbanks,
    )

    rowsize!(fig.layout, 0, legend_rowsize)

    return fig
end

function line_plot_facet!(
    gl,
    df;
    ax_position = (1, 1),
    num_metrics = 1,
    facet_title = "",
    ylims = nothing,
)
    grouped_dfs = groupby(df, :noise_specification)
    unique_tests = unique(df[!, :test_specification])
    test_labels = map(
        test ->
            "[$(test.sensitivity), $(test.test_result_lag)]",
        unique_tests,
    )

    ax = Axis(
        gl[1, 1];
        title = facet_title,
        xlabel = "Test Specification ([Sensitivity/Specificity, Lag])",
        ylabel = "EWS Accuracy",
        xticks = (collect(1:length(test_labels)), test_labels),
    )

    for gdf in grouped_dfs
        lines!(
            ax,
            1:nrow(gdf),
            gdf.accuracy,
        )
    end

    if !(isnothing(ylims))
        ylims!(ax, ylims)
    end

    if ax_position[2] != 1
        hideydecorations!(ax)
    end

    nrows = ceil(num_metrics / 2)
    if ax_position[1] != nrows && ax_position[2] + 2 <= num_metrics
        hidexdecorations!(ax)
    end
    return nothing
end

function heatmap_noise_description(
    noise_specification::T
) where {T<:PoissonNoiseSpecification}
    return string(
        "Poisson noise: ", noise_specification.noise_mean_scaling, "x Measles"
    )
end

function heatmap_noise_description(
    noise_specification::T
) where {T<:DynamicalNoiseSpecification}
    mean_vaccination_coverage = mean([
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
    ])

    return string(
        "Dynamical noise: mean vaccination $(mean_vaccination_coverage)"
    )
end
