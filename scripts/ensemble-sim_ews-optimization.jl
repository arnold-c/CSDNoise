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
    PoissonNoiseSpecification(1.0),
    PoissonNoiseSpecification(7.0),
    DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.8734),
    DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.102),
]

test_specification_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(0.95, 0.95, 0),
    IndividualTestSpecification(0.97, 0.97, 0),
    IndividualTestSpecification(0.98, 0.98, 0),
    IndividualTestSpecification(0.99, 0.99, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

percent_tested_vec = [1.0]

#%%
ews_method_vec = [
    # Centered,
    Backward
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
    Reff_start
    # Reff_end
]
ews_threshold_window_vec = [ExpandingThresholdWindow]
ews_threshold_percentile_vec = [
    collect(0.5:0.05:0.85)..., collect(0.9:0.02:0.98)..., 0.99
]
ews_consecutive_thresholds_vec = [collect(2:1:30)...]
ews_threshold_burnin_vec = [
    # Day(50),
    Year(5)
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
        Type{ExpandingThresholdWindow},Type{RollingThresholdWindow}
    }[],
    ews_threshold_burnin = Union{Dates.Day,Dates.Year}[],
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
    test_noise_specification = DynamicalNoiseSpecification(
        5.0,
        7,
        14,
        "in-phase",
        0.15,
        # 0.8734
        0.102,
    )
    # test_noise_specification = PoissonNoiseSpecification(1.0)
    # test_specification = IndividualTestSpecification(0.8, 0.8, 0)
    # test_specification = IndividualTestSpecification(0.9, 0.9, 0)
    # test_specification = IndividualTestSpecification(0.95, 0.95, 0)
    # test_specification = IndividualTestSpecification(0.97, 0.97, 0)
    test_specification = IndividualTestSpecification(0.98, 0.98, 0)
    # test_specification = IndividualTestSpecification(0.99, 0.99, 0)
    # test_specification = IndividualTestSpecification(1.0, 1.0, 0)
    # test_ews_metric = "mean"
    test_ews_metric = "variance"
    # test_ews_metric = "skewness"
    test_ews_metric_specification = EWSMetricSpecification(
        Backward, Day(7), Week(52), 1
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
        :test_specification => ByRow(==(test_specification)),
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
    vec_of_null_detection_index_vec
), noisearr = simulate_ews_survival_data(
        test_df,
        ensemble_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = test_ews_metric,
    )
end

#%%
if debug_Reff_plots
    selected_sim = 96
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
        vec_of_testarr[1][:, 5, selected_sim],
        vec_of_null_testarr[1][:, 5, selected_sim],
        noisearr[:, selected_sim],
        vec_of_ews_vals_vec[1][selected_sim],
        vec_of_null_ews_vals_vec[1][selected_sim],
        vec(vec_of_exceed_thresholds[1][selected_sim, 1]),
        vec(vec_of_null_exceed_thresholds[1][selected_sim, 1]),
        vec(vec_of_threshold_percentiles[1][selected_sim, 1]),
        vec(vec_of_null_threshold_percentiles[1][selected_sim, 1]),
        vec_of_detection_index_vec[1][selected_sim],
        vec_of_null_detection_index_vec[1][selected_sim],
        ensemble_time_specification;
        xlims = (0, 12),
        burnin_vline = test_ews_threshold_burnin,
        plottitle = plottitle,
    )
end
