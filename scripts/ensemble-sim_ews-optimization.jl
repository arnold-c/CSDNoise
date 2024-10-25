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

include(srcdir("makie-plotting-setup.jl"))

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

#%%
logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt")

noise_specification_vec = [
    PoissonNoiseSpecification(1.0),
    PoissonNoiseSpecification(8.0),
]

test_specification_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

percent_tested_vec = [1.0]

#%%
ews_method_vec = [
    # Centered,
    Backward
]
ews_aggregation_days_vec = [
    7,
    14,
    28,
]
ews_bandwidth_days_vec = [52 * 7]
ews_lag_days_vec = [1]

ews_metric_specification_vec = create_combinations_vec(
    calculate_bandwidth_and_return_ews_metric_spec,
    (
        ews_method_vec,
        ews_aggregation_days_vec,
        ews_bandwidth_days_vec,
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
ews_threshold_window_vec = [Main.Expanding]
ews_threshold_percentile_vec = [collect(0.9:0.02:0.98)..., 0.99]
ews_consecutive_thresholds_vec = [collect(2:1:30)...]
ews_threshold_burnin_vec = [50]

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
    ews_threshold_window = EWSThresholdWindow[],
    ews_threshold_burnin = Int[],
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
    subset_optimal_parameters = [:ews_threshold_burnin => ByRow(==(50))],
)

#%%
@unpack ews_df = load_most_recent_hyperparam_file(
    "ews-hyperparam-gridsearch.jld2",
    outdir("ensemble", "ews-hyperparam-optimization"),
)

#%%
optimal_heatmap_df = optimal_ews_heatmap_df(
    optimal_ews_df;
    tiebreaker_preference = "specificity",
)

#%%
test_df = subset(
    optimal_heatmap_df,
    :ews_metric_specification =>
        ByRow(==(EWSMetricSpecification(Backward, 7, 52, 1))),
    :ews_enddate_type => ByRow(==(Reff_start)),
    :ews_threshold_burnin => ByRow(==(50)),
    :ews_threshold_window => ByRow(==(Main.Expanding)),
    :noise_specification => ByRow(==(PoissonNoiseSpecification(1.0))),
)

#%%
optimal_ews_heatmap_plot(
    subset(
        optimal_heatmap_df,
        :ews_metric_specification =>
            ByRow(==(EWSMetricSpecification(Backward, 28, 52 ÷ 4, 1))),
        :ews_enddate_type => ByRow(==(Reff_start)),
        :ews_threshold_burnin => ByRow(==(50)),
        :ews_threshold_window => ByRow(==(Main.Expanding)),
        :noise_specification => ByRow(==(PoissonNoiseSpecification(1.0))),
    ),
)

#%%
survival_df = simulate_ews_survival_data(
    test_df,
    ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec;
    ews_metric = "mean",
)

#%%

#%%
subset_survival_df = subset(
    survival_df,
    :test_specification => ByRow(==(IndividualTestSpecification(1.0, 1.0, 0))),
)

#%%
detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
    subset_survival_df
)

#%%
ews_survival_plot(
    detection_survival_vecs,
    null_survival_vecs,
    subset_survival_df.enddate;
)

#%%
hist(enddate_vec)
