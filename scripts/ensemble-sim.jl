#%%
using DrWatson
@quickactivate "CSDNoise"

using ProgressMeter
using Chain

using CSDNoise

#%%
model_types_vec = [("seasonal-infectivity-import", "tau-leaping")]

#%%
N_vec = [500_000]
nsims_vec = [100]
init_states_prop_dict = [
    Dict(:s_prop => 0.05, :e_prop => 0.00, :i_prop => 0.00, :r_prop => 0.95)
]

ensemble_state_p_vec = create_combinations_vec(
    StateParameters, (N_vec, init_states_prop_dict)
)

#%%
tmin_vec = [0.0]
tstep_vec = [1.0]
nyears_vec = [10]
tmax_vec = nyears_vec .* 365.0

time_p_vec = vec(
    map(
        Iterators.product(tmin_vec, tstep_vec, tmax_vec)
    ) do (tmin, tstep, tmax)
        SimTimeParameters(; tmin = tmin, tmax = tmax, tstep = tstep)
    end,
)

#%%
beta_force_vec = [0.0]
annual_births_per_k_vec = [27]
seed = 1234

#%%
latent_per_days_vec = [LATENT_PER_DAYS]
dur_inf_days_vec = [DUR_INF_DAYS]
R_0_vec = collect(16.0)
sigma_vec = 1 ./ latent_per_days_vec
gamma_vec = 1 ./ dur_inf_days_vec
vaccination_coverage_vec = [0.8]

#%%
ensemble_spec_vec = create_ensemble_spec_combinations(
    beta_force_vec,
    [cos],
    sigma_vec,
    gamma_vec,
    annual_births_per_k_vec,
    R_0_vec,
    vaccination_coverage_vec,
    N_vec,
    init_states_prop_dict,
    model_types_vec,
    time_p_vec,
    nsims_vec,
)

#%%
outbreak_threshold_vec = [5]
min_outbreak_dur_vec = [30]
min_outbreak_size_vec = [500]

outbreak_spec_vec = create_combinations_vec(
    OutbreakSpecification,
    (outbreak_threshold_vec, min_outbreak_dur_vec, min_outbreak_size_vec),
)

outbreak_spec_dict = Vector{Dict}(undef, length(outbreak_spec_vec))
for (i, spec) in pairs(outbreak_spec_vec)
    outbreak_spec_dict[i] = Dict{Symbol,Any}(:outbreak_spec => spec)
end

#%%
# TODO: How to make the mean of the noise a proportion of the mean of incidence
# when incidence depends on the dynamics parameters?
# Could pass variables to ensemble function and calculate each simulations and
# scenario's noise mean, but that would break implementation using NoiseSpecification
# struct currently
poisson_noise_mean_scaling_vec = [8.0]

poisson_noise_spec_vec = create_combinations_vec(
    PoissonNoiseSpecification,
    (["poisson"], poisson_noise_mean_scaling_vec)
)

dynamical_noise_R0 = [5.0]
dynamical_noise_latent_period = [7]
dynamical_noise_duration_infection = [14]
dynamical_noise_correlation = ["in-phase", "out-of-phase", "none"]
dynamical_noise_mean_scaling_vec = [1.0]
dynamical_noise_spec_vec = create_combinations_vec(
    DynamicalNoiseSpecification,
    (
        ["dynamical"],
        dynamical_noise_R0,
        dynamical_noise_latent_period,
        dynamical_noise_duration_infection,
        dynamical_noise_correlation,
        dynamical_noise_mean_scaling_vec,
    ),
)

noise_spec_vec = vcat(poisson_noise_spec_vec, dynamical_noise_spec_vec)

#%%
alertthreshold_vec = collect(10:4:50)
moveavglag_vec = [7]
perc_clinic_vec = [0.6]
perc_clinic_test_vec = [collect(0.1:0.1:0.6)..., 1.0]
alert_method_vec = ["inferred_movingavg"]

outbreak_detection_spec_vec = create_combinations_vec(
    OutbreakDetectionSpecification,
    (
        alertthreshold_vec,
        moveavglag_vec,
        perc_clinic_vec,
        perc_clinic_test_vec,
        alert_method_vec,
    ),
)

#%%
test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 3),
    IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

#%%
ews_method_vec = ["centered", "backward"]
ews_bandwidth_vec = collect(7:1:7)

ews_spec_vec = create_combinations_vec(
    EWSMetricSpecification,
    (ews_method_vec, ews_bandwidth_vec)
)

#%%
base_param_dict = @dict(
    ensemble_spec = ensemble_spec_vec,
    seed = seed,
)

sol_param_dict = dict_list(
    base_param_dict
)

for dict in sol_param_dict
    dict[:quantile_vec] = [95]
    dict[:outbreak_spec_dict] = outbreak_spec_dict
    dict[:noise_spec_vec] = noise_spec_vec
    dict[:outbreak_detection_spec_vec] = outbreak_detection_spec_vec
    dict[:test_spec_vec] = test_spec_vec
    dict[:ews_spec_vec] = ews_spec_vec
end

#%%
run_ensemble_jump_prob(sol_param_dict; force = true)
