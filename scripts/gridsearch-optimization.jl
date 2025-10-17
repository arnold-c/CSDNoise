#%%
using CSDNoiseCore
using Dates
using StructArrays

#%%
time_parameters = SimTimeParameters(
    burnin = Day(5.0 * 365.0),
    tmin = 0.0,
    tmax = 365.0 * 20.0,
    tstep = 1.0
)

births_per_k_pop = 27.0

burnin_target_Reff = 0.9

common_disease_dynamics_parameters = CommonDiseaseDynamicsParameters(
    births_per_k_pop = births_per_k_pop,
    nsims = 1000,
    burnin_target_Reff = burnin_target_Reff
)

initial_Reff = 0.8

#%%
measles_R0 = 16.0

measles_initial_s_prop = round(initial_Reff / measles_R0; digits = 2)

# Set initial population proportions so that initial Reff = 0.8
measles_population_state_parameters = StateParameters(
    500_000,
    Dict(
        :s_prop => measles_initial_s_prop,
        :e_prop => 0.0,
        :i_prop => 0.0,
        :r_prop => 1.0 - measles_initial_s_prop
    )
)

# Calculate the maximum post-burnin vaccination coverage to reach Reff = 1 within
# 3 quarters of the simulation time after the end of the burnin to ensure the tipping
# point is reached (assuming the evaluation period is ReffStart and the value of Reff
# at the end of the burnin is a bit less than the maximum target Reff - burnin_target_Reff).
measles_min_post_burnin_vaccination_coverage,
    measles_max_post_burnin_vaccination_coverage = calculate_post_burnin_vaccination_range(
    time_parameters,
    measles_population_state_parameters,
    measles_R0;
    births_per_k_pop = births_per_k_pop,
    target_Reff = 1.0,
    estimated_post_burnin_Reff = initial_Reff,
    adjust_vaccination_coverage = -0.1,
    vaccination_bounds_spread = 0.2,
    reach_target_prop_of_remaining_simulation = 0.75,
    digits = 1
)

# @guerraBasicReproductionNumber2017 @gastanaduyMeasles2019
measles_dynamics_parameters = TargetDiseaseDynamicsParameters(;
    R_0 = measles_R0,
    latent_period_days = 10.0,
    infectious_duration_days = 8.0,
    beta_force = 0.0,
    seasonality = SeasonalityFunction(CosineSeasonality()),
    min_post_burnin_vaccination_coverage = measles_min_post_burnin_vaccination_coverage,
    max_post_burnin_vaccination_coverage = measles_max_post_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage = 1.0
)

# Choose Rubella-like parameters for the dynamic noise in the measles simulations
# @papadopoulosEstimatesBasicReproduction2022 @RubellaCDCYellow.
rubella_dynamical_noise_spec = DynamicalNoiseSpecification(
    R_0 = 5.0,
    latent_period = 7.0,
    duration_infection = 14.0,
    correlation = "in-phase",
    poisson_component = 0.15,
)

measles_ensemble_specification = create_ensemble_specifications(
    "measles",
    time_parameters,
    measles_population_state_parameters,
    measles_dynamics_parameters,
    common_disease_dynamics_parameters,
    rubella_dynamical_noise_spec
)

#%%
covid_R0 = 3.3

covid_initial_s_prop = round(initial_Reff / covid_R0; digits = 2)

# Set initial population proportions so that initial Reff = 0.8
covid_population_state_parameters = StateParameters(
    500_000,
    Dict(
        :s_prop => covid_initial_s_prop,
        :e_prop => 0.0,
        :i_prop => 0.0,
        :r_prop => 1.0 - covid_initial_s_prop
    )
)

covid_min_post_burnin_vaccination_coverage,
    covid_max_post_burnin_vaccination_coverage = calculate_post_burnin_vaccination_range(
    time_parameters,
    covid_population_state_parameters,
    covid_R0;
    births_per_k_pop = births_per_k_pop,
    target_Reff = 1.0,
    estimated_post_burnin_Reff = initial_Reff,
    adjust_vaccination_coverage = -0.1,
    vaccination_bounds_spread = 0.2,
    reach_target_prop_of_remaining_simulation = 0.75,
    digits = 1
)

# R0 = 3.3 (https://pmc.ncbi.nlm.nih.gov/articles/PMC7280807/)
# Duration of infection = 10 days (https://pmc.ncbi.nlm.nih.gov/articles/PMC7547320/)
# Latent period = 5.5 days (https://academic.oup.com/cid/article/74/9/1678/6359063?login=false)
covid_dynamics_parameters = TargetDiseaseDynamicsParameters(;
    R_0 = covid_R0,
    latent_period_days = 5.5,
    infectious_duration_days = 10.0,
    beta_force = 0.0,
    seasonality = SeasonalityFunction(CosineSeasonality()),
    min_post_burnin_vaccination_coverage = covid_min_post_burnin_vaccination_coverage,
    max_post_burnin_vaccination_coverage = covid_max_post_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage = 1.0
)

# Choose seasonal influenza-like parameters for the dynamic noise in the COVID-19
# and moderate disease simulations
# R0 = 1.28 (https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-14-480)
# Duration of infection = 4.8 days (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
# Latent period = 1 day (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
influenza_dynamical_noise_spec = DynamicalNoiseSpecification(
    R_0 = 1.28,
    latent_period = 1.0,
    duration_infection = 4.8,
    correlation = "in-phase",
    poisson_component = 0.15,
)

covid_ensemble_specification = create_ensemble_specifications(
    "covid-19",
    time_parameters,
    covid_population_state_parameters,
    covid_dynamics_parameters,
    common_disease_dynamics_parameters,
    influenza_dynamical_noise_spec
)

#%%
# Moderate disease with R0 = 5
## Common values
moderate_disease_R0 = 5.0
moderate_disease_latent_period = 5.0

moderate_disease_initial_s_prop = round(initial_Reff / moderate_disease_R0; digits = 2)

# Set initial population proportions so that initial Reff = 0.8
moderate_disease_population_state_parameters = StateParameters(
    500_000,
    Dict(
        :s_prop => moderate_disease_initial_s_prop,
        :e_prop => 0.0,
        :i_prop => 0.0,
        :r_prop => 1.0 - moderate_disease_initial_s_prop
    )
)

moderate_disease_min_post_burnin_vaccination_coverage,
    moderate_disease_max_post_burnin_vaccination_coverage = calculate_post_burnin_vaccination_range(
    time_parameters,
    moderate_disease_population_state_parameters,
    moderate_disease_R0;
    births_per_k_pop = births_per_k_pop,
    target_Reff = 1.0,
    estimated_post_burnin_Reff = initial_Reff,
    adjust_vaccination_coverage = -0.1,
    vaccination_bounds_spread = 0.2,
    reach_target_prop_of_remaining_simulation = 0.75,
    digits = 1
)

#%%
## Medium R0 with measles beta
moderate_disease_measles_like_infectious_duration = calculate_infectious_duration(
    moderate_disease_R0,
    measles_ensemble_specification.emergent_dynamics_parameter_specification.beta_mean,
    moderate_disease_latent_period,
    calculate_mu(births_per_k_pop)
)

moderate_disease_measles_like_dynamics_parameters = TargetDiseaseDynamicsParameters(;
    R_0 = moderate_disease_R0,
    latent_period_days = moderate_disease_latent_period,
    infectious_duration_days = moderate_disease_measles_like_infectious_duration,
    beta_force = 0.0,
    seasonality = SeasonalityFunction(CosineSeasonality()),
    min_post_burnin_vaccination_coverage = moderate_disease_min_post_burnin_vaccination_coverage,
    max_post_burnin_vaccination_coverage = moderate_disease_max_post_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage = 1.0
)


moderate_disease_measles_like_ensemble_specification = create_ensemble_specifications(
    "moderate-disease-measles-beta",
    time_parameters,
    moderate_disease_population_state_parameters,
    moderate_disease_measles_like_dynamics_parameters,
    common_disease_dynamics_parameters,
    influenza_dynamical_noise_spec
)

#%%
## Medium R0 with COVID beta
moderate_disease_covid_like_infectious_duration = calculate_infectious_duration(
    moderate_disease_R0,
    covid_ensemble_specification.emergent_dynamics_parameter_specification.beta_mean,
    moderate_disease_latent_period,
    calculate_mu(births_per_k_pop)
)

moderate_disease_covid_like_dynamics_parameters = TargetDiseaseDynamicsParameters(;
    R_0 = moderate_disease_R0,
    latent_period_days = moderate_disease_latent_period,
    infectious_duration_days = moderate_disease_covid_like_infectious_duration,
    beta_force = 0.0,
    seasonality = SeasonalityFunction(CosineSeasonality()),
    min_post_burnin_vaccination_coverage = moderate_disease_min_post_burnin_vaccination_coverage,
    max_post_burnin_vaccination_coverage = moderate_disease_max_post_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage = 1.0
)


moderate_disease_covid_like_ensemble_specification = create_ensemble_specifications(
    "moderate-disease-covid-beta",
    time_parameters,
    moderate_disease_population_state_parameters,
    moderate_disease_covid_like_dynamics_parameters,
    common_disease_dynamics_parameters,
    influenza_dynamical_noise_spec
)


#%%
ensemble_specification_vec = [
    measles_ensemble_specification,
    covid_ensemble_specification,
    moderate_disease_measles_like_ensemble_specification,
    moderate_disease_covid_like_ensemble_specification,
]

#%%
noise_level_vec = [1.0, 4.0, 8.0, 12.0]
noise_type_description_vec = [:static, :dynamic]

test_specification_vec = [
    IndividualTestSpecification(val, val, 0) for val in [0.8, 0.9, 0.95:0.01:1.0...]
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
ews_threshold_quantile_vec = collect(0.5:0.01:0.99)
ews_consecutive_thresholds_vec = collect(1:1:20)

#%%
specification_vecs = GridSearchSpecificationVecs(
    ensemble_specification_vec = ensemble_specification_vec,
    noise_level_vec = noise_level_vec,
    noise_type_description_vec = noise_type_description_vec,
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
gridsearch_scenarios = ews_hyperparam_gridsearch_structvector(
    specification_vecs;
    force = true,
    save_results = true,
    save_checkpoints = true,
    disable_time_check = true
)
