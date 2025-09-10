#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack

using CSDNoise
using StructArrays
using SumTypes
using Try
using DataFrames
using StatsBase: StatsBase
using Random: Random
using Distributions: Distributions
using StyledStrings
using Dates

using ProgressMeter

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

#%%
ensemble_dynamics_parameters_sa = StructVector(
    get_ensemble_file(
        ensemble_specification
    )["dynamics_parameters"],
)

null_dynamics_parameters_sa = StructVector(
    get_ensemble_file(
        null_specification
    )["dynamics_parameters"],
)

@assert ensemble_dynamics_parameters_sa.burnin_vaccination_coverage ==
    null_dynamics_parameters_sa.burnin_vaccination_coverage ==
    null_dynamics_parameters_sa.vaccination_coverage

#%%
ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

null_single_scenario_inc_file = get_ensemble_file(
    null_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

null_single_incarr = null_single_scenario_inc_file["ensemble_inc_arr"]
null_single_periodsum_vecs = null_single_scenario_inc_file["ensemble_thresholds_vec"]

ensemble_single_Reff_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_arr"]

ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

null_single_Reff_arr = get_ensemble_file(
    null_specification
)["ensemble_Reff_arr"]

null_single_Reff_thresholds_vec = get_ensemble_file(
    null_specification
)["ensemble_Reff_thresholds_vec"]

#%%
logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt")

noise_specification_vec = [
    NoiseSpecification(PoissonNoise(1.0)),
    NoiseSpecification(PoissonNoise(7.0)),
    NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.8734)),
    NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.102)),
]

test_specification_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(0.95, 0.95, 0),
    IndividualTestSpecification(0.96, 0.96, 0),
    IndividualTestSpecification(0.97, 0.97, 0),
    IndividualTestSpecification(0.98, 0.98, 0),
    IndividualTestSpecification(0.99, 0.99, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

percent_tested_vec = [1.0]

#%%
ews_method_vec = [
    # Centered,
    Backward,
]
ews_aggregation_vec = [
    Day(7),
    Day(14),
    Day(28),
]
ews_bandwidth_vec = [Week(52)]
ews_lag_days_vec = [1]

ews_metric_specification_vec = create_combinations_vec(
    EWSMetricSpecification,
    (
        ews_method_vec,
        ews_aggregation_vec,
        ews_bandwidth_vec,
        ews_lag_days_vec,
    ),
)

ews_metric_vec = [
    "autocorrelation",
    "autocovariance",
    "coefficient_of_variation",
    "index_of_dispersion",
    "kurtosis",
    "mean",
    "skewness",
    "variance",
]

#%%
ews_enddate_type_vec = [
    Reff_start,
    # Reff_end
]
ews_threshold_window_vec = [ExpandingThresholdWindow]
ews_threshold_percentile_vec = collect(0.5:0.01:0.99)
ews_consecutive_thresholds_vec = [collect(2:1:30)...]
ews_threshold_burnin_vec = [
    # Day(50),
    Year(5),
]

#%%
specification_vecs = (;
    noise_specification_vec,
    test_specification_vec,
    percent_tested_vec,
    ews_metric_specification_vec,
    ews_enddate_type_vec,
    ews_threshold_window_vec,
    ews_threshold_burnin_vec,
    ews_threshold_percentile_vec,
    ews_consecutive_thresholds_vec,
    ews_metric_vec,
)

specification_vec_tuples = (
    noise_specification = NoiseSpecification[],
    test_specification = IndividualTestSpecification[],
    percent_tested = Float64[],
    ews_metric_specification = EWSMetricSpecification[],
    ews_enddate_type = EWSEndDateType[],
    ews_threshold_window = Union{
        Type{ExpandingThresholdWindow}, Type{RollingThresholdWindow},
    }[],
    ews_threshold_burnin = Union{Dates.Day, Dates.Year}[],
    ews_threshold_percentile = Float64[],
    ews_consecutive_thresholds = Int[],
    ews_metric = String[],
)

#%%
optimal_ews_df = ews_hyperparam_optimization(
    specification_vecs,
    (;
        ensemble_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs,
    );
    filedir = outdir("ensemble", "ews-hyperparam-optimization"),
    gridsearch_filename_base = "ews-hyperparam-gridsearch.jld2",
    optimization_filename_base = "ews-hyperparam-optimization.jld2",
    logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
    force = false,
    return_df = true,
    specification_vec_tuples = specification_vec_tuples,
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
        :ews_threshold_window,
        :ews_threshold_burnin,
    ],
    disable_time_check = false,
)

#%%
create_optimal_ews_plots(
    optimal_ews_df,
    ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec,
    ["specificity", "speed"];
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
    ],
    output_format = "png",
    force_heatmap = true,
    force_survival = false,
)

#%%
debug_Reff_plots = true

if debug_Reff_plots
    # test_noise_specification = DynamicalNoiseSpecification(
    #     5.0,
    #     7,
    #     14,
    #     "in-phase",
    #     0.15,
    #     0.8734,
    #     # 0.102,
    # )
    test_noise_specification = NoiseSpecification(PoissonNoise(1.0))
    # test_ews_metric = "mean"
    test_ews_metric = "autocovariance"
    # test_ews_metric = "skewness"
    test_ews_metric_specification = EWSMetricSpecification(
        # Backward, Day(7), Week(52), 1
        Backward, Day(28), Week(52), 1
    )
    test_ews_enddate_type = Reff_start
    test_ews_threshold_burnin = Dates.Year(5)
    test_ews_threshold_window = ExpandingThresholdWindow
    percent_tested = 1.0
    test_tiebreaker_preference = "specificity"

    optimal_heatmap_df = optimal_ews_heatmap_df(
        optimal_ews_df;
        tiebreaker_preference = test_tiebreaker_preference,
    )

    test_df = subset(
        optimal_heatmap_df,
        :ews_metric_specification =>
            ByRow(==(test_ews_metric_specification)),
        :ews_enddate_type => ByRow(==(test_ews_enddate_type)),
        :ews_threshold_burnin => ByRow(==(test_ews_threshold_burnin)),
        :ews_threshold_window => ByRow(==(test_ews_threshold_window)),
        :noise_specification => ByRow(==(test_noise_specification)),
        # :test_specification => ByRow(==(test_specification)),
        :percent_tested => ByRow(==(percent_tested)),
    )

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
            vec_of_null_detection_index_vec,
        ), noisearr = simulate_ews_survival_data(
        test_df,
        ensemble_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = test_ews_metric,
    )
end;

#%%
if debug_Reff_plots
    # include(srcdir("makie-plotting-setup.jl"))
    include(srcdir("cairomakie-plotting-setup.jl"))

    # test_specification = IndividualTestSpecification(0.8, 0.8, 0)
    test_specification = IndividualTestSpecification(0.9, 0.9, 0)
    # test_specification = IndividualTestSpecification(0.95, 0.95, 0)
    # test_specification = IndividualTestSpecification(0.97, 0.97, 0)
    # test_specification = IndividualTestSpecification(0.98, 0.98, 0)
    # test_specification = IndividualTestSpecification(0.99, 0.99, 0)
    # test_specification = IndividualTestSpecification(1.0, 1.0, 0)

    test_ind = findfirst(
        t -> t == test_specification, unique(survival_df.test_specification)
    )

    selected_sim = 90

    plottitle =
        "Noise: $(get_noise_magnitude_description(test_noise_specification)), Percent Tested: $(percent_tested), $(split(string(test_ews_enddate_type), "::")[1])" *
        "\nEWS Metric: $(test_ews_metric), $(get_ews_metric_specification_description(test_ews_metric_specification)), Threshold Burnin: $(test_ews_threshold_burnin), Tiebreaker: $(test_tiebreaker_preference)" *
        "\nEWS calculated for shaded region"

    hyperparam_debugging_Reff_plot(
        ensemble_single_incarr[:, 1, selected_sim],
        null_single_incarr[:, 1, selected_sim],
        ensemble_single_Reff_arr[:, selected_sim],
        null_single_Reff_arr[:, selected_sim],
        ensemble_single_Reff_thresholds_vec[selected_sim],
        null_single_Reff_thresholds_vec[selected_sim],
        ensemble_single_periodsum_vecs[selected_sim],
        null_single_periodsum_vecs[selected_sim],
        Symbol(test_ews_metric),
        vec_of_testarr[test_ind][:, 5, selected_sim],
        vec_of_null_testarr[test_ind][:, 5, selected_sim],
        noisearr[:, selected_sim],
        vec_of_ews_vals_vec[test_ind][selected_sim],
        vec_of_null_ews_vals_vec[test_ind][selected_sim],
        vec(vec_of_exceed_thresholds[test_ind][selected_sim, 1]),
        vec(vec_of_null_exceed_thresholds[test_ind][selected_sim, 1]),
        vec(vec_of_threshold_percentiles[test_ind][selected_sim, 1]),
        vec(vec_of_null_threshold_percentiles[test_ind][selected_sim, 1]),
        vec_of_detection_index_vec[test_ind][selected_sim],
        vec_of_null_detection_index_vec[test_ind][selected_sim],
        ensemble_time_specification;
        xlims = (0, 12),
        burnin_vline = test_ews_threshold_burnin,
        plottitle = plottitle,
    )
end

#%%
selected_sim = 20

ews_specification =
    vec_of_ews_vals_vec[test_ind][selected_sim].ews_specification
@unpack aggregation = ews_specification
aggregation_int = Dates.days(aggregation)

aggregated_inc_vec = aggregate_timeseries(
    ensemble_single_incarr[:, 1, selected_sim],
    aggregation,
)

aggregated_Reff_vec = aggregate_Reff_vec(
    ensemble_single_Reff_arr[:, selected_sim],
    aggregation,
)

aggregated_Reff_thresholds =
    ensemble_single_Reff_thresholds_vec[selected_sim] .÷ aggregation_int

aggregated_null_Reff_vec = aggregate_Reff_vec(
    null_single_Reff_arr[:, selected_sim],
    aggregation,
)

aggregated_null_Reff_thresholds =
    null_single_Reff_thresholds_vec[selected_sim] .÷ aggregation_int

aggregated_outbreak_bounds =
    ensemble_single_periodsum_vecs[selected_sim][
    (
        ensemble_single_periodsum_vecs[selected_sim][:, 4] .== 1
    ), [1, 2],
] .÷
    aggregation_int

aggregated_test_vec = aggregate_timeseries(
    vec_of_testarr[test_ind][:, 5, selected_sim],
    aggregation,
)

aggregated_noise_vec = aggregate_timeseries(
    noisearr[:, selected_sim],
    aggregation,
)

#%%
null_reff_schematic = Reff_plot(
    aggregated_null_Reff_vec,
    aggregated_null_Reff_thresholds,
    ensemble_time_specification,
    vec_of_ews_vals_vec[test_ind][selected_sim];
    xlims = (0, 8),
    ylims_Reff = (0.6, 1.2),
    legends = false,
    Reff_colormap = [BASE_COLOR],
)

save(
    plotsdir("null-reff-schematic.png"),
    null_reff_schematic;
    size = (600, 800 / 3),
)

#%%
reff_schematic = Reff_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    ensemble_time_specification,
    vec_of_ews_vals_vec[test_ind][selected_sim];
    xlims = (0, 8),
    ylims_Reff = (0.6, 1.2),
    legends = false,
)

save(
    plotsdir("autocovariance-reff-schematic.png"),
    reff_schematic;
    size = (600, 800 / 3),
)

#%%
reff_inc_schematic = Reff_inc_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_inc_vec,
    aggregated_outbreak_bounds,
    ensemble_time_specification,
    vec_of_ews_vals_vec[test_ind][selected_sim];
    xlims = (0, 12),
    ylims_inc = (0, 950),
    legends = false,
)

save(
    plotsdir("autocovariance-reff-inc-schematic.png"),
    reff_inc_schematic;
    size = (1300, (2 * 800 / 3) * 1.038),
)

#%%
reff_noise_schematic = Reff_noise_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_noise_vec,
    ensemble_time_specification,
    vec_of_ews_vals_vec[test_ind][selected_sim];
    xlims = (0, 12),
    ylims_inc = (0, 950),
    legends = false,
)

save(
    plotsdir("autocovariance-reff-noise-schematic.png"),
    reff_noise_schematic;
    size = (1300, (2 * 800 / 3) * 1.038),
)

#%%
reff_noise_inc_schematic = Reff_noise_inc_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_inc_vec,
    aggregated_outbreak_bounds,
    aggregated_noise_vec,
    ensemble_time_specification,
    vec_of_ews_vals_vec[test_ind][selected_sim];
    xlims = (0, 12),
    legends = false,
    ylims_inc = (0, 950),
)

save(
    plotsdir("autocovariance-reff-noise-inc-schematic.png"),
    reff_noise_inc_schematic;
    size = (1300, (2 * 800 / 3) * 1.038),
)

#%%
reff_noise_inc_test_schematic = Reff_noise_inc_test_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_test_vec,
    aggregated_inc_vec,
    aggregated_outbreak_bounds,
    aggregated_noise_vec,
    ensemble_time_specification,
    vec_of_ews_vals_vec[test_ind][selected_sim];
    xlims = (0, 12),
    legends = false,
    ylims_inc = (0, 950),
    test_positive_color = "#4B2EDF",
)

save(
    plotsdir("autocovariance-reff-noise-inc-test-schematic.png"),
    reff_noise_inc_test_schematic;
    size = (1300, (2 * 800 / 3) * 1.038),
)

#%%
test_enddate = Try.unwrap(
    calculate_ews_enddate(
        ensemble_single_Reff_thresholds_vec[selected_sim],
        test_ews_enddate_type,
    ),
)

inc_ews = EWSMetrics(
    test_ews_metric_specification,
    @view(ensemble_single_incarr[1:test_enddate, 1, selected_sim])
)

inc_ews_no_exceedance_schematic = Reff_ews_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_inc_vec,
    inc_ews,
    Symbol(test_ews_metric),
    aggregated_outbreak_bounds,
    ensemble_time_specification;
    xlims = (0, 12),
    ylims_inc = (0, 950),
    legends = false,
    plot_tau = false,
    plot_detection_index = false,
)

save(
    plotsdir("inc-autocovariance-schematic_no-exceedance.png"),
    inc_ews_no_exceedance_schematic;
    size = (1300, 800),
)

#%%
ews_no_exceedance_schematic = Reff_ews_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_inc_vec,
    aggregated_test_vec,
    aggregated_noise_vec,
    vec_of_ews_vals_vec[test_ind][selected_sim],
    Symbol(test_ews_metric),
    aggregated_outbreak_bounds,
    ensemble_time_specification;
    xlims = (0, 12),
    ylims_inc = (0, 950),
    legends = false,
    plot_tau = false,
    plot_detection_index = false,
    test_positive_color = "#4B2EDF",
)

save(
    plotsdir("autocovariance-schematic_no-exceedance.png"),
    ews_no_exceedance_schematic;
    size = (1300, 800),
)

#%%
ews_schematic = Reff_ews_plot(
    aggregated_Reff_vec,
    aggregated_Reff_thresholds,
    aggregated_inc_vec,
    aggregated_test_vec,
    aggregated_noise_vec,
    vec_of_ews_vals_vec[test_ind][selected_sim],
    Symbol(test_ews_metric),
    aggregated_outbreak_bounds,
    vec(vec_of_exceed_thresholds[test_ind][selected_sim, 1]),
    vec_of_detection_index_vec[test_ind][selected_sim],
    ensemble_time_specification;
    xlims = (0, 12),
    ylims_inc = (0, 950),
    legends = false,
    plot_tau = false,
    plot_detection_index = false,
    test_positive_color = "#4B2EDF",
)

save(
    plotsdir("autocovariance-schematic.png"),
    ews_schematic;
    size = (1300, 800),
)

#%%
if debug_Reff_plots
    ews_survival_plot(
        survival_df;
        noise_specification_vec = [
            NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.8734)),
        ],
        test_specification_vec = [
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
        ],
    )
end

#%%
if debug_Reff_plots
    auc_df = DataFrame(
        "ews_metric" => String[],
        "test_specification" => IndividualTestSpecification[],
        # "ews_metric_specification" => EWSMetricSpecification[],
        # "ews_enddate_type" => EWSEndDateType[],
        "auc" => Float64[],
        "auc_magnitude" => Float64[],
        "emergent_tau" => Vector{Float64}[],
        "null_tau" => Vector{Float64}[],
    )

    unique_tests = unique(survival_df.test_specification)
    subset_start_ind = Int64(round(Dates.days(Year(5))))

    for test_ind in eachindex(vec_of_ews_vals_vec)
        test_specification = unique_tests[test_ind]
        subset_testarr = vec_of_testarr[test_ind][
            :, 5, :,
        ]
        subset_null_testarr = vec_of_null_testarr[test_ind][
            :, 5, :,
        ]

        enddate_vec = map(axes(subset_testarr, 2)) do sim
            Try.unwrap(
                calculate_ews_enddate(
                    ensemble_single_Reff_thresholds_vec[sim],
                    ews_enddate_type,
                ),
            )
        end

        emergent_sv =
            map(axes(subset_testarr, 2)) do sim
            enddate = enddate_vec[sim]
            EWSMetrics(
                ews_metric_specification,
                subset_testarr[subset_start_ind:enddate, sim],
            )
        end |>
            v -> StructVector(v)

        null_sv =
            map(axes(subset_null_testarr, 2)) do sim
            enddate = enddate_vec[sim]
            EWSMetrics(
                ews_metric_specification,
                subset_null_testarr[subset_start_ind:enddate, sim],
            )
        end |>
            v -> StructVector(v)

        for metric in ews_metrics
            metric_tau = Symbol(metric * "_tau")
            emergent_tau = getproperty(emergent_sv, metric_tau)
            null_tau = getproperty(null_sv, metric_tau)
            auc = calculate_auc(
                emergent_tau,
                null_tau,
            )
            auc_magnitude = abs(auc - 0.5)
            push!(
                auc_df,
                (
                    ews_metric = metric,
                    test_specification = test_specification,
                    # ews_metric_specification = ews_metric_specification,
                    # ews_enddate_type = ews_enddate_type,
                    auc = auc,
                    auc_magnitude,
                    emergent_tau = emergent_tau,
                    null_tau = null_tau,
                ),
            )
        end
    end

    tau_auc_heatmap(
        auc_df,
        :auc,
        :auc;
        baseline_test = test_specification,
    )
end

#%%
tau_auc_heatmap(
    auc_df,
    :auc,
    :auc;
    baseline_test = test_specification,
    plottitle = "",
    ylabel = "",
    legendwidth = 20,
    xticklabelsize = 22,
    legendticklabelsize = 22,
    legendsize = 22,
    fontsize = 22,
    digits = 3,
)
