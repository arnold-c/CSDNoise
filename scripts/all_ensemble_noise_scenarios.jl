## High R0 (Measles + Rubella) scenario
#
# - Measles
#     - R0 = 16 ()
#     - Duration of infection = 8 days ()
#     - Latent period = 10 days ()
#     - @guerraBasicReproductionNumber2017 @gastanaduyMeasles2019
# - Rubella
#     - R0 = 5 ()
#     - Duration of infection = 14 days ()
#     - Latent period = 7 days ()
#     - @papadopoulosEstimatesBasicReproduction2022 @RubellaCDCYellow.
#
# ## Low R0 (COVID + ILI) scenario
#
# - COVID
#     - R0 = 3.3 (https://pmc.ncbi.nlm.nih.gov/articles/PMC7280807/)
#     - Duration of infection = 10 days (https://pmc.ncbi.nlm.nih.gov/articles/PMC7547320/)
#     - Latent period = 5.5 days (https://academic.oup.com/cid/article/74/9/1678/6359063?login=false)
# - Seasonal influenza
#     - R0 = 1.28 (https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-14-480)
#     - Duration of infection = 4.8 days (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
#     - Latent period = 1 day (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
#
# ## Medium R0 with Measles beta
#
# - Target Disease
#     - R0 = 5
#     - beta = 4.002367123287671e-6
#     - Duration of infection = 1 / calculated_gamma
#     - Latent period = 5 days
# - Seasonal influenza
#     - R0 = 1.28 (https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-14-480)
#     - Duration of infection = 4.8 days (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
#     - Latent period = 1 day (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
#
# ## Medium R0 with COVID beta
#
# - Target Disease
#     - R0 = 5
#     - beta = 6.604882191780822e-7
#     - Duration of infection = 1 / calculated_gamma
#     - Latent period = 5 days
# - Seasonal influenza
#     - R0 = 1.28 (https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-14-480)
#     - Duration of infection = 4.8 days (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
#     - Latent period = 1 day (https://academic.oup.com/aje/article/167/7/775/83777#86199749)
#
#

#%%
using CSDNoise
using Dates
using GLMakie

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
    nsims = 100,
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

#%%
calculate_vaccination_rate_to_achieve_Reff(
    time_parameters,
    population_state_parameters,
    1.0,
    16.0,
    calculate_mu(27)
)

#%%
measles_ensemble_data = generate_single_ensemble(
    measles_ensemble_specification;
    seed = 1234
)

#%%
measles_Reff_thresholds = calculate_all_Reff_thresholds(measles_ensemble_data)

enddates_vec = calculate_all_ews_enddates(
    measles_Reff_thresholds,
    EWSEndDateType(ReffStart())
)

inc = trim_seir_results(
    measles_ensemble_data.emergent_seir_run,
    enddates_vec,
)

#%%
target_noise_scaling = 8.0

optim_res = optimize_dynamic_noise_params_wrapper(
    measles_ensemble_specification,
    inc,
    enddates_vec,
    target_noise_scaling,
    NoiseVaccinationOptimizationParameters(;
        atol = 0.1
    );
    verbose = false
)

#%%
optim_noise_run = CSDNoise.recreate_noise_vecs(
    measles_ensemble_specification,
    enddates_vec,
    optim_res.location[1],
)

fig = Figure()
ax = Axis(fig[1, 1]; title = "Rubella Noise Incidence ($(target_noise_scaling)x)")
for res in optim_noise_run.incidence[1:end]
    lines!(ax, res)
end
ax2 = Axis(fig[2, 1]; title = "Measles Incidence")
for res in inc.incidence[1:end]
    lines!(ax2, res)
end
ax3 = Axis(fig[3, 1]; title = "Measles Reff")
for res in inc.Reff[1:end]
    lines!(ax3, res)
end
ax4 = Axis(fig[4, 1]; title = "Noise Reff")
for res in optim_noise_run.Reff[1:end]
    lines!(ax4, res)
end
linkxaxes!(ax, ax2, ax3, ax4)
linkyaxes!(ax, ax2)
display(fig)

#%%
optim_res.location[1]

target_noise_scaling * mean_measles - optim_noise_run.mean_noise

(target_noise_scaling * mean_measles - optim_noise_run.mean_noise)^2 == optim_res[1]


#%%
enddates_vec = calculate_all_ews_enddates(
    measles_ensemble_data.Reff_thresholds,
    EWSEndDateType(ReffStart())
)

#%%
optim_noise_run = CSDNoise.recreate_noise_vecs(
    measles_dynamical_noise_spec,
    optim_res.location[1],
    measles_spec[1],
    enddates_vec
)

#%%
filtered_seir_results = trim_seir_results(
    measles_ensemble_data.seir_results,
    enddates_vec
)

# Calculate filtered mean incidence up to endpoints
mean_measles = calculate_mean_incidence(filtered_seir_results)

mean_measles * 8.0 - optim_noise_run.mean_noise

#%%
@report_opt target_modules = (CSDNoise,) trim_seir_results(
    measles_ensemble_data.seir_results,
    enddates_vec
)

#%%
measles_noise_spec_8x = NoiseSpecification(
    DynamicalNoise(
        measles_dynamical_noise_spec.R0,
        measles_dynamical_noise_spec.latent_period,
        measles_dynamical_noise_spec.duration_infection,
        measles_dynamical_noise_spec.correlation,
        measles_dynamical_noise_spec.poisson_component,
        calculate_min_max_vaccination_range(
            0.8796, # mean vaccination for 4x
            0.2
        )...
    )
)

measles_noise_8x_ensemble_spec = EnsembleSpecification(
    StateParameters(;
        N = measles_spec[1].state_parameters.init_states.N,
        s_prop = 0.022,
        e_prop = measles_spec[1].state_parameters.init_state_props.e_prop,
        i_prop = measles_spec[1].state_parameters.init_state_props.i_prop,
    ),
    measles_spec[1].dynamics_parameter_specification,
    measles_spec[1].time_parameters,
    measles_spec[1].nsims,
    measles_spec[1].dirpath
)

#%%
measles_noise_8x_vecs = create_noise_vecs(
    measles_noise_spec_8x,
    measles_noise_8x_ensemble_spec,
    enddates_vec
)

inc = trim_seir_results(
    generate_single_ensemble(measles_spec[1]).seir_results,
    enddates_vec,
)

using GLMakie

#%%
fig = Figure()
ax = Axis(fig[1, 1]; title = "Rubella Noise Incidence (8x)")
for res in optim_noise_run.incidence
    lines!(ax, res)
end
ax2 = Axis(fig[2, 1]; title = "Measles Incidence")
for res in inc.incidence
    lines!(ax2, res)
end
ax3 = Axis(fig[3, 1]; title = "Measles Reff")
for res in inc.Reff
    lines!(ax3, res)
end
linkxaxes!(ax, ax2, ax3)
linkyaxes!(ax, ax2)
display(fig)

#%%
for scaling in (1.0, 2.0, 4.0, 6.0)
    measles_res = calculate_dynamic_vaccination_coverage(
        scaling,  # target_scaling
        measles_mean,
        measles_dynamical_noise_spec,
        measles_spec[1],
        enddates
    )
    # covid_res = calculate_dynamic_vaccination_coverage_multistart(
    #     scaling,  # target_scaling
    #     mean_covid,
    #     covid_dynamical_noise_spec,
    #     covid_spec[1]
    # )
    println("\nScaling: $scaling\n")
    println("-"^20)
    println("Measles")
    println(measles_res)
    # println("COVID")
    # println(covid_res)
    println()
    println("="^60)
end


#%%
fig = Figure()
ax = Axis(fig[1, 1])
for k in axes(measles_emergent_incidence_arr, 3)
    lines!(measles_emergent_incidence_arr[:, 1, k]; color = (:blue, 0.01))
    lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
mean_measles_vec = map(
    i -> StatsBase.median(@view(measles_emergent_incidence_arr[i, 1, :])),
    axes(measles_emergent_incidence_arr, 1)
)
mean_noise_vec = map(
    i -> StatsBase.median(@view(measles_noise_arr[i, :])),
    axes(measles_noise_arr, 1)
)
lines!(mean_measles_vec; color = :black)
lines!(mean_noise_vec; color = :black)
ax2 = Axis(fig[2, 1])
for k in axes(emergent_incidence_arr, 3)
    lines!(emergent_incidence_arr[:, 1, k]; color = (:blue, 0.01))
    lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
mean_measles_vec2 = map(
    i -> StatsBase.median(@view(emergent_incidence_arr[i, 1, :])),
    axes(emergent_incidence_arr, 1)
)
lines!(mean_measles_vec2; color = :black)
lines!(mean_noise_vec; color = :black)
display(fig)

#%%
fig = Figure()
ax = Axis(fig[1, 1])
for k in axes(measles_emergent_incidence_arr, 3)
    lines!(measles_emergent_incidence_arr[:, 1, k]; color = (:blue, 0.01))
    # lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
mean_measles_vec = map(
    i -> StatsBase.median(@view(measles_emergent_incidence_arr[i, 1, :])),
    axes(measles_emergent_incidence_arr, 1)
)
# mean_noise_vec = map(
#     i -> StatsBase.median(@view(measles_noise_arr[i, :])),
#     axes(measles_noise_arr, 1)
# )
lines!(mean_measles_vec; color = :black)
# lines!(mean_noise_vec; color = :black)
ax2 = Axis(fig[2, 1])
for k in eachindex(measles_incs)
    lines!(measles_incs[k]; color = (:blue, 0.01))
    # lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
lines!(mean_measles_incs; color = :black)
# lines!(mean_noise_vec; color = :black)
linkaxes!(ax, ax2)
display(fig)

#%%
for i in eachindex(measles_incs)
    local enddate = measles_enddates[i]
    @assert measles_incs[i] == measles_emergent_incidence_arr[1:enddate, 1, i] "failed $i: $(measles_incs[i]) != $(measles_emergent_incidence_arr[1:enddate, 1, i])"
    println("i = $i: $(round(mean(measles_incs[i]); digits = 2)), $(round(mean(measles_emergent_incidence_arr[1:enddate, 1, i]); digits = 2))")
end
