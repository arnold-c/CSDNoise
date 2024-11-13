#%%
using DrWatson
@quickactivate "CSDNoise"
using CSDNoise
using Try: Try
using DataFrames
using Dates
using StructArrays
using StatsBase: StatsBase

include(srcdir("cairomakie-plotting-setup.jl"))
CairoMakie.activate!(; type = "svg")

#%%
ensemble_model_type = ("seasonal-infectivity-import", "tau-leaping")

burnin_years = 5
nyears = 20
burnin_time = 365.0 * burnin_years
ensemble_time_specification = SimTimeParameters(;
    burnin = 365.0 * burnin_years, tmin = 0.0, tmax = 365.0 * nyears,
    tstep = 1.0,
)

ensemble_state_specification = StateParameters(
    500_000,
    Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
)

mu = calculate_mu(27)
beta_mean = calculate_beta(
    R0, GAMMA, mu, 1, ensemble_state_specification.init_states.N
)
epsilon = calculate_import_rate(
    mu, R0, ensemble_state_specification.init_states.N
)

min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
    0.9,
    burnin_years * 2,
    ensemble_state_specification.init_states.S,
    ensemble_state_specification.init_states.N,
    R0,
    mu,
)

max_burnin_vaccination_coverage = 1.0
min_vaccination_coverage = 0.0
max_vaccination_coverage = 0.8

ensemble_dynamics_specification = DynamicsParameterSpecification(
    beta_mean,
    0.0,
    cos,
    SIGMA,
    GAMMA,
    mu,
    27,
    epsilon,
    R0,
    min_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage,
    0.6,
    0.8,
)

null_dynamics_specification = DynamicsParameterSpecification(
    map(
        pn -> getproperty(ensemble_dynamics_specification, pn),
        filter(
            name ->
                name != :min_vaccination_coverage &&
                    name != :max_vaccination_coverage,
            propertynames(ensemble_dynamics_specification),
        ),
    )...,
    nothing,
    nothing,
)

#%%
ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

null_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    null_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

ensemble_outbreak_specification = OutbreakSpecification(
    5, 30, 500
)

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

null_single_scenario_inc_file = get_ensemble_file(
    null_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
null_single_incarr = null_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

#%%
ews_metrics = [
    "autocorrelation",
    "autocovariance",
    "coefficient_of_variation",
    "index_of_dispersion",
    "kurtosis",
    "mean",
    "skewness",
    "variance",
]

ews_metric_specification = EWSMetricSpecification(
    Backward, Day(28), Day(364), 1
)
ews_threshold_burnin = Year(5)
ews_enddate_type = Reff_start
ews_threshold_window = ExpandingThresholdWindow
percent_tested = 1.0

optimization_filepath = get_most_recent_hyperparam_filepath(
    "ews-hyperparam-optimization.jld2",
    outdir("ensemble", "ews-hyperparam-optimization"),
)

previous_optimal_ews_df = load(Try.unwrap(optimization_filepath))["optimal_ews_df"]

subset_optimal_df = subset(
    previous_optimal_ews_df,
    :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
    :ews_metric_specification =>
        ByRow(==(ews_metric_specification)),
    :percent_tested => ByRow(==(percent_tested)),
    :ews_metric_specification => ByRow(==(ews_metric_specification)),
    :ews_enddate_type => ByRow(==(ews_enddate_type)),
)

#%%
function plot_noise_filename(
    noise_specification::T
) where {T<:PoissonNoiseSpecification}
    return "poisson_$(noise_specification.noise_mean_scaling)x"
end

function plot_noise_filename(
    noise_specification::T
) where {T<:DynamicalNoiseSpecification}
    mean_vaccination_coverage = mean([
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
    ])
    return "dynamical_$(mean_vaccination_coverage)"
end

function plot_test_filename(test_specification)
    return "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)"
end

#%%
gdfs = groupby(
    subset_optimal_df,
    [
        :noise_specification
    ],
)

for gdf in gdfs
    @assert length(unique(gdf.noise_specification)) == 1
    noise_specification = gdf.noise_specification[1]

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

    noise_description = heatmap_noise_description(noise_specification)

    optimal_heatmap_plot = optimal_ews_heatmap_plot(
        optimal_heatmap_df;
        subtitle = noise_description,
        colormap = :Blues,
        colorrange = (0.5, 0.8),
        textcolorthreshold = 0.68,
    )

    optimal_heatmap_plot_name =
        "optimal_heatmap_" * plot_noise_filename(noise_specification) *
        ".svg"

    heatmap_plotdir = projectdir(
        "manuscript", "manuscript_files", "plots", "optimal-threshold-heatmaps"
    )

    mkpath(heatmap_plotdir)

    save(
        joinpath(heatmap_plotdir, optimal_heatmap_plot_name),
        optimal_heatmap_plot;
        size = (1700, 1000),
    )

    test_df = subset(
        optimal_heatmap_df,
        :ews_metric_specification =>
            ByRow(==(ews_metric_specification)),
        :ews_enddate_type => ByRow(==(ews_enddate_type)),
        :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
        :ews_threshold_window => ByRow(==(ews_threshold_window)),
        :noise_specification => ByRow(==(noise_specification)),
        :percent_tested => ByRow(==(percent_tested)),
    )

    for (i, survival_ews_metric) in pairs(["mean", "variance"])
        survival_df, (
        vec_of_testarr,
        vec_of_null_testarr,
        vec_of_ews_vals_vec,
        vec_of_null_ews_vals_vec,
        vec_of_exceed_thresholds,
        vec_of_null_exceed_thresholds,
        vec_of_threshold_percentiles,
        vec_of_null_threshold_percentiles,
        vec_of_detection_index_vec,
        vec_of_null_detection_index_vec
), noisearr = simulate_ews_survival_data(
            test_df,
            ensemble_specification,
            ensemble_single_incarr,
            null_single_incarr,
            ensemble_single_Reff_thresholds_vec;
            ews_metric = survival_ews_metric,
        )

        for test_specification in [
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(0.90, 0.90, 0),
        ]
            detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
                subset(
                    survival_df,
                    :test_specification => ByRow(==(test_specification)),
                ),
            )

            survival_plot = ews_survival_plot(
                detection_survival_vecs,
                null_survival_vecs,
                survival_df.enddate;
                ews_aggregation = ews_metric_specification.aggregation,
                burnin = ews_threshold_burnin,
            )

            survival_plot_name =
                "survival_" *
                "test-" * plot_test_filename(test_specification) *
                "_ews-" * survival_ews_metric * "_" *
                plot_noise_filename(noise_specification) *
                ".svg"

            survival_plotdir = projectdir(
                "manuscript",
                "manuscript_files",
                "plots",
                "survival",
                "test-$(plot_test_filename(test_specification))",
            )

            mkpath(survival_plotdir)

            save(
                joinpath(survival_plotdir, survival_plot_name),
                survival_plot,
            )
        end

        if i == 1
            ews_df = DataFrame(
                "ews_metric" => String[],
                "test_specification" => IndividualTestSpecification[],
                "ews_metric_specification" => EWSMetricSpecification[],
                "ews_enddate_type" => EWSEndDateType[],
                "ews_metric_value" => Float64[],
                "ews_metric_vector" => Vector{Float64}[],
            )

            unique_tests = unique(survival_df.test_specification)

            for test_ind in eachindex(vec_of_ews_vals_vec)
                test_specification = unique_tests[test_ind]
                sv = StructVector(
                    convert(Vector{EWSMetrics}, vec_of_ews_vals_vec[test_ind])
                )
                for metric in ews_metrics
                    simulation_tau_heatmap_df!(
                        ews_df,
                        sv,
                        metric;
                        individual_test_specification = test_specification,
                        ews_metric_specification = ews_metric_specification,
                        ews_enddate_type = ews_enddate_type,
                        statistic_function = StatsBase.mean,
                    )
                end
            end

            tau_heatmap = tycho_tau_heatmap_plot(
                subset(
                    ews_df,
                    :ews_enddate_type => ByRow(==(ews_enddate_type)),
                    :ews_metric_specification =>
                        ByRow(==(ews_metric_specification)),
                ),
            )

            tau_heatmap_plot_name =
                "tau-heatmap_" *
                plot_noise_filename(noise_specification) *
                ".svg"

            tau_heatmap_plotdir = projectdir(
                "manuscript",
                "manuscript_files",
                "plots",
                "tau-heatmaps",
                "full-length",
            )

            mkpath(tau_heatmap_plotdir)

            save(
                joinpath(tau_heatmap_plotdir, tau_heatmap_plot_name),
                tau_heatmap,
            )
        end

        if i == 1
            ews_df = DataFrame(
                "ews_metric" => String[],
                "test_specification" => IndividualTestSpecification[],
                "ews_metric_specification" => EWSMetricSpecification[],
                "ews_enddate_type" => EWSEndDateType[],
                "ews_metric_value" => Float64[],
                "ews_metric_vector" => Vector{Float64}[],
            )

            unique_tests = unique(survival_df.test_specification)
            subset_start_ind = Int64(round(Dates.days(Year(5))))

            for test_ind in eachindex(vec_of_ews_vals_vec)
                test_specification = unique_tests[test_ind]
                subset_testarr = vec_of_testarr[test_ind][
                    :, 5, :,
                ]
                sv =
                    map(axes(subset_testarr, 3)) do sim
                        enddate = Try.unwrap(
                            calculate_ews_enddate(
                                ensemble_single_Reff_thresholds_vec[sim],
                                ews_enddate_type,
                            ),
                        )

                        EWSMetrics(
                            ews_metric_specification,
                            subset_testarr[subset_start_ind:enddate, sim],
                        )
                    end |>
                    v -> StructVector(v)

                for metric in ews_metrics
                    simulation_tau_heatmap_df!(
                        ews_df,
                        sv,
                        metric;
                        individual_test_specification = test_specification,
                        ews_metric_specification = ews_metric_specification,
                        ews_enddate_type = ews_enddate_type,
                        statistic_function = StatsBase.mean,
                    )
                end
            end

            tau_heatmap = tycho_tau_heatmap_plot(
                subset(
                    ews_df,
                    :ews_enddate_type => ByRow(==(ews_enddate_type)),
                    :ews_metric_specification =>
                        ByRow(==(ews_metric_specification)),
                ),
            )

            tau_heatmap_plot_name =
                "tau-heatmap_" *
                plot_noise_filename(noise_specification) *
                ".svg"

            tau_heatmap_plotdir = projectdir(
                "manuscript",
                "manuscript_files",
                "plots",
                "tau-heatmaps",
                "after-burnin",
            )

            mkpath(tau_heatmap_plotdir)

            save(
                joinpath(tau_heatmap_plotdir, tau_heatmap_plot_name),
                tau_heatmap,
            )
        end
    end
end

#%%
lineplot_df = similar(gdfs[1], 0)
for gdf in gdfs, ewsmetric in ["mean", "variance", "autocovariance"]
    prepare_line_plot_df!(lineplot_df, gdf, ewsmetric)
end

line_plot(
    lineplot_df;
    plottitle = "Mean",
)
