#%%
using CSDNoise
using Dates

#%%
time_parameters = SimTimeParameters(
    burnin = Day(5.0 * 365.0),
    tmin = 0.0,
    tmax = 365.0 * 20.0,
    tstep = 1.0
)

population_state_parameters = StateParameters(
    500_000,
    Dict(
        :s_prop => 0.05,
        :e_prop => 0.0,
        :i_prop => 0.0,
        :r_prop => 0.95
    )
)

measles_dynamics_parameters = TargetDiseaseDynamicsParameters(;
    R_0 = 16.0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.0,
    seasonality = SeasonalityFunction(CosineSeasonality()),
    min_post_burnin_vaccination_coverage = 0.6,
    max_post_burnin_vaccination_coverage = 0.8,
    max_burnin_vaccination_coverage = 1.0
)

common_disease_dynamics_parameters = CommonDiseaseDynamicsParameters(
    births_per_k_pop = 27.0,
    nsims = 1000,
    burnin_target_Reff = 0.9
)

measles_dynamical_noise_spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 0.15,
)

measles_ensemble_specification = create_ensemble_specifications(
    time_parameters,
    population_state_parameters,
    measles_dynamics_parameters,
    common_disease_dynamics_parameters,
    measles_dynamical_noise_spec
)

ensemble_specification_vec = [measles_ensemble_specification]

#%%
noise_level_vec = [1.0]

test_specification_vec = [
    IndividualTestSpecification(0.9, 0.9, 0),
]
percent_tested_vec = [1.0]

ews_method_vec = [EWSMethod(Backward())]
ews_aggregation_vec = [Day(28)]
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

ews_enddate_type_vec = [EWSEndDateType(ReffStart())]
ews_threshold_window_vec = [EWSThresholdWindowType(ExpandingThresholdWindow())]
ews_threshold_quantile_vec = collect(0.5:0.1:0.99)
ews_consecutive_thresholds_vec = collect(2:1:20)

#%%
specification_vecs = GridSearchSpecificationVecs(
    ensemble_specification_vec = ensemble_specification_vec,
    noise_level_vec = noise_level_vec,
    test_specification_vec = test_specification_vec,
    percent_tested_vec = percent_tested_vec,
    ews_metric_specification_vec = ews_metric_specification_vec,
    ews_enddate_type_vec = ews_enddate_type_vec,
    ews_threshold_window_vec = ews_threshold_window_vec,
    ews_threshold_quantile_vec = ews_threshold_quantile_vec,
    ews_consecutive_thresholds_vec = ews_consecutive_thresholds_vec,
    ews_metric_vec = ews_metric_vec,
)

#%%
scenarios = create_gridsearch_scenarios_structvector(specification_vecs)

missing_scenarios = find_missing_scenarios(
    scenarios,
    StructVector(OptimizationResult[])
)

evaluate_gridsearch_scenarios(missing_scenarios)

generate_single_ensemble(
    scenarios[1].ensemble_specification
)

#%%
ews_hyperparam_gridsearch_structvector(
    specification_vecs;
    # File management
    filedir = outdir("ensemble", "ews-hyperparam-gridsearch"),
    gridsearch_filename_base = "ews-hyperparam-gridsearch-structvector.jld2",
    gridsearch_output_filepath = joinpath(
        filedir,
        string(Dates.now()) * "_" * gridsearch_filename_base,
    ),
    # Execution configuration
    executor = FLoops.SequentialEx(),
    # Control options
    force = false,
    return_results = true,
    save_results = true,
    save_checkpoints = false,
    verbose = true,
    disable_time_check = false,
)


#%%
ensemble_outbreak_specification = OutbreakSpecification(
    5, 30, 500
)
