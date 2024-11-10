#%%
using DrWatson
using CSDNoise
using Try: Try
using DataFrames
using Dates

include(srcdir("cairomakie-plotting-setup.jl"))
CairoMakie.activate!(; type = "svg")

#%%
optimization_filepath = get_most_recent_hyperparam_filepath(
    "ews-hyperparam-optimization.jld2",
    outdir("ensemble", "ews-hyperparam-optimization"),
)

previous_optimal_ews_df = load(Try.unwrap(optimization_filepath))["optimal_ews_df"]

subset_optimal_df = subset(
    previous_optimal_ews_df,
    :ews_threshold_burnin => ByRow(==(Year(5))),
    :ews_metric_specification =>
        ByRow(==(EWSMetricSpecification(Backward, Day(7), Day(364), 1))),
)

#%%
create_optimal_ews_plots(
    subset_optimal_df,
    ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec,
    ["specificity"];
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
    ],
    ews_metrics = [
        "autocorrelation",
        "autocovariance",
        "coefficient_of_variation",
        "index_of_dispersion",
        "kurtosis",
        "mean",
        "skewness",
        "variance",
    ],
    base_plotpath = joinpath(plotsdir(), "ensemble"),
    output_format = "svg",
)

#%%
gdfs = groupby(
    subset_optimal_df,
    [
        :noise_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_burnin,
    ],
)

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

function heatmap_noise_filename(
    noise_specification::T
) where {T<:PoissonNoiseSpecification}
    return "poisson_$(noise_specification.noise_mean_scaling)x"
end

function heatmap_noise_filename(
    noise_specification::T
) where {T<:DynamicalNoiseSpecification}
    mean_vaccination_coverage = mean([
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
    ])
    return "dynamical_$(mean_vaccination_coverage)"
end

for gdf in gdfs
    optimal_heatmap_df = optimal_ews_heatmap_df(
        gdf;
        tiebreaker_preference = "specificity",
        optimal_grouping_parameters = [
            :noise_specification,
            :test_specification,
            :percent_tested,
            :ews_metric_specification,
            :ews_enddate_type,
            :ews_threshold_window,
            :ews_threshold_burnin,
            :ews_metric,
        ],
    )

    noise_specification = optimal_heatmap_df.noise_specification[1]
    noise_description = heatmap_noise_description(noise_specification)

    optimal_heatmap_plot = optimal_ews_heatmap_plot(
        optimal_heatmap_df; subtitle = noise_description
    )

    optimal_heatmap_plot_name =
        "optimal_heatmap_" * heatmap_noise_filename(noise_specification) *
        ".svg"

    heatmap_plotdir = projectdir("manuscript", "manuscript_files", "plots")

    save(
        joinpath(heatmap_plotdir, optimal_heatmap_plot_name),
        optimal_heatmap_plot,
    )
end
