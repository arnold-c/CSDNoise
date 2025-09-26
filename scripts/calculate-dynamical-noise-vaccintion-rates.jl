#%%
using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using UnPack: @unpack
using StatsBase: StatsBase

#%%
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

ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

#%%
dynamical_noise_spec = (;
    R0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 0.15,
)

#%%
calculate_mean_dynamical_noise(
    5.0,
    7,
    14,
    "in-phase",
    0.15,
    0.8734,
    0.2,
    ensemble_specification,
)

#%%
emergent_incidence_arr = get_ensemble_file(
    ensemble_specification, OutbreakSpecification(5, 30, 500)
)["emergent_incidence_arr"]

mean_measles = StatsBase.mean(emergent_incidence_arr[:, 1, :])

#%%
calculate_dynamic_vaccination_coverage(
    7,
    mean_measles,
    dynamical_noise_spec,
    ensemble_specification;
    maxiters = 20,
    vaccination_mean_range = [0.0, 1.0],
    max_vaccination_range = 0.2,
    atol = 0.01,
    showprogress = true,
)

calculate_min_max_vaccination_range(
    0.8734, 0.2
)

#%%
dynamical_noise_coverages = map(
    ((scale, coverage_range),) -> calculate_dynamic_vaccination_coverage(
        scale,
        mean_measles,
        dynamical_noise_spec,
        ensemble_specification;
        maxiters = 20,
        vaccination_mean_range = coverage_range,
        max_vaccination_range = 0.2,
        atol = 0.01,
        showprogress = true,
    ),
    zip(
        [1, 7],
        [[0.8, 0.9], [0.05, 0.15]],
    ),
)

#%%
map(
    ((t, c),) -> t[2] - c * mean_measles,
    zip(dynamical_noise_coverages, [1, 7]),
)
