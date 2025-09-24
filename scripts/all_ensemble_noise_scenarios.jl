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
using StatsBase


#%%
## Measles + Rubella
measles_spec = create_ensemble_specs(
    EnsembleSpecsParameters(
        burnin_years = 5,
        nyears = 20,
        annual_births_per_k = ANNUAL_BIRTHS_PER_K,
        ensemble_state_specification = StateParameters(
            500_000,
            Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
        ),
        R_0 = R0,
        gamma = GAMMA,
        sigma = SIGMA,
        target_Reff = 0.9,
        target_years = 10,
        min_vaccination_coverage = 0.6,
        max_vaccination_coverage = 0.8,
        nsims = 100
    )
)

#%%
covid_spec = create_ensemble_specs(
    EnsembleSpecsParameters(
        burnin_years = 5,
        nyears = 20,
        annual_births_per_k = ANNUAL_BIRTHS_PER_K,
        ensemble_state_specification = StateParameters(
            500_000,
            Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
        ),
        R_0 = 3.3,
        gamma = 1 / 10,
        sigma = 1 / 5.5,
        target_Reff = 0.9,
        target_years = 10,
        min_vaccination_coverage = 0.1,
        max_vaccination_coverage = 0.12,
        nsims = 100
    )
)

# Is this the max vaccination rate for the emergent time series??
calculate_vaccination_rate_to_achieve_Reff(
    0.9,
    10,
    500_000 * 0.05,
    500_000,
    3.3,
    MU
)

#%%
measles_spec[1].dynamics_parameter_specification.min_burnin_vaccination_coverage
measles_spec[1].dynamics_parameter_specification.min_vaccination_coverage
measles_spec[1].dynamics_parameter_specification.max_vaccination_coverage

covid_spec[1].dynamics_parameter_specification.min_burnin_vaccination_coverage
covid_spec[1].dynamics_parameter_specification.min_burnin_vaccination_coverage

measles_spec[3] == covid_spec[3]
ensemble_outbreak_specification = measles_spec[3]

#%%
measles_ensemble_data = generate_single_ensemble(measles_spec[1]; seed = 1234)

#%%
measles_inc_arr = getindex.(measles_ensemble_data.ensemble_inc_vecs, 1)
measles_enddates = Vector{Int64}(undef, length(measles_ensemble_data.ensemble_Reff_thresholds_vec));
measles_incs = Vector{Vector{Int64}}(undef, length(measles_enddates))

for i in eachindex(measles_enddates)
    local enddate = Try.@? calculate_ews_enddate(measles_ensemble_data.ensemble_Reff_thresholds_vec[i], EWSEndDateType(Reff_start()))
    measles_enddates[i] = enddate
    measles_incs[i] = measles_inc_arr[1:enddate, i]
end

#%%
mean_measles_incs = Vector{Float64}(undef, maximum(measles_enddates))
for i in eachindex(mean_measles_incs)
    # Only include vectors that have data at position i
    local values_at_i = Int64[
        measles_incs[j][i] for j in eachindex(measles_incs) if i <= length(measles_incs[j])
    ]
    mean_measles_incs[i] = mean(values_at_i)
end

#%%
# Generate incidence arrays for the outbreak specification using the incidence vectors
measles_ensemble_inc_arr = create_inc_infec_arr(
    measles_ensemble_data[:ensemble_inc_vecs], ensemble_outbreak_specification
)[1]


mean_measles = StatsBase.mean(measles_ensemble_inc_arr[:, 1, :])
mean_measles = StatsBase.mean(measles_ensemble_inc_arr[1:3651, 1, :])
mean_measles2 = StatsBase.mean(ensemble_inc_arr[:, 1, :])

#%%
covid_ensemble_data = generate_single_ensemble(covid_spec[1]; seed = 1234)

# Generate incidence arrays for the outbreak specification using the incidence vectors
covid_ensemble_inc_arr = create_inc_infec_arr(
    covid_ensemble_data[:ensemble_inc_vecs], ensemble_outbreak_specification
)[1]


mean_covid = StatsBase.mean(covid_ensemble_inc_arr[:, 1, :])

#%%
measles_dynamical_noise_spec = (;
    R0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 0.15,
)

covid_dynamical_noise_spec = (;
    R0 = 1.28,
    latent_period = 1,
    duration_infection = 4.8,
    correlation = "in-phase",
    poisson_component = 0.15,
)


#%%
using JLD2
ensemble_inc_arr = JLD2.load("/Users/cfa5228/Documents/Repos/CSDNoise/out/ensemble/seasonal-infectivity-import/tau-leaping/N_500000/r_0.95/nsims_100/R0_16.0/latent_period_10.0/infectious_period_8.0/min_burnin_vaccination_coverage_0.9269/max_burnin_vaccination_coverage_1.0/min_vaccination_coverage_0.6/max_vaccination_coverage_0.8/births_per_k_27/beta_force_0.0/burnin_1825.0/tmax_7300.0/tstep_1.0/min_outbreak_dur_30/min_outbreak_size_500/outbreak_threshold_5/ensemble-incidence-array.jld2")["ensemble_inc_arr"]

#%%
size(measles_ensemble_inc_arr) == size(ensemble_inc_arr)

#%%
measles_dynamical_noise_coverages = map(
    ((scale, coverage_range),) -> calculate_dynamic_vaccination_coverage(
        scale,
        mean_measles,
        measles_dynamical_noise_spec,
        measles_spec[1];
        maxiters = 20,
        vaccination_mean_range = coverage_range,
        max_vaccination_range = 0.2,
        atol = 0.02,
        showprogress = false,
    ),
    zip(
        [1, 7],
        [[0.1, 0.9], [0.05, 0.15]],
    ),
)

#%%
measles_dynamical_noise_coverages2 = map(
    ((scale, coverage_range),) -> calculate_dynamic_vaccination_coverage(
        scale,
        mean_measles2,
        measles_dynamical_noise_spec,
        measles_spec[1];
        maxiters = 20,
        vaccination_mean_range = coverage_range,
        max_vaccination_range = 0.2,
        atol = 0.02,
        showprogress = false,
    ),
    zip(
        [1, 7],
        [[0.1, 0.9], [0.05, 0.15]],
    ),
)

# #%%
# map(
#     ((t, c),) -> t[2] - c * mean_measles,
#     zip(measles_dynamical_noise_coverages, [1, 7]),
# )

#%%
for scaling in (1.0, 2.0, 4.0, 6.0)
    measles_res = calculate_dynamic_vaccination_coverage_multistart(
        scaling,  # target_scaling
        mean_measles,
        measles_dynamical_noise_spec,
        measles_spec[1]
    )
    covid_res = calculate_dynamic_vaccination_coverage_multistart(
        scaling,  # target_scaling
        mean_covid,
        covid_dynamical_noise_spec,
        covid_spec[1]
    )
    println("\nScaling: $scaling\n")
    println("-"^20)
    println("Measles")
    println(measles_res)
    println("COVID")
    println(covid_res)
    println()
    println("="^60)
end

#%%
measles_noise_spec_4x = NoiseSpecification(
    DynamicalNoise(
        measles_dynamical_noise_spec.R0,
        measles_dynamical_noise_spec.latent_period,
        measles_dynamical_noise_spec.duration_infection,
        measles_dynamical_noise_spec.correlation,
        measles_dynamical_noise_spec.poisson_component,
        calculate_min_max_vaccination_range(
            0.0154, # mean vaccination for 4x
            0.2
        )...
    )
)

measles_noise_arr = create_noise_arr(
    measles_noise_spec_4x,
    measles_spec[1];
    seed = 1234,
)[1]

#%%
fig = Figure()
ax = Axis(fig[1, 1])
for k in axes(measles_ensemble_inc_arr, 3)
    lines!(measles_ensemble_inc_arr[:, 1, k]; color = (:blue, 0.01))
    lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
mean_measles_vec = map(
    i -> StatsBase.median(@view(measles_ensemble_inc_arr[i, 1, :])),
    axes(measles_ensemble_inc_arr, 1)
)
mean_noise_vec = map(
    i -> StatsBase.median(@view(measles_noise_arr[i, :])),
    axes(measles_noise_arr, 1)
)
lines!(mean_measles_vec; color = :black)
lines!(mean_noise_vec; color = :black)
ax2 = Axis(fig[2, 1])
for k in axes(ensemble_inc_arr, 3)
    lines!(ensemble_inc_arr[:, 1, k]; color = (:blue, 0.01))
    lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
mean_measles_vec2 = map(
    i -> StatsBase.median(@view(ensemble_inc_arr[i, 1, :])),
    axes(ensemble_inc_arr, 1)
)
lines!(mean_measles_vec2; color = :black)
lines!(mean_noise_vec; color = :black)
display(fig)

#%%
fig = Figure()
ax = Axis(fig[1, 1])
for k in axes(measles_ensemble_inc_arr, 3)
    lines!(measles_ensemble_inc_arr[:, 1, k]; color = (:blue, 0.01))
    # lines!(measles_noise_arr[:, k]; color = (:red, 0.01))
end
mean_measles_vec = map(
    i -> StatsBase.median(@view(measles_ensemble_inc_arr[i, 1, :])),
    axes(measles_ensemble_inc_arr, 1)
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
    @assert measles_incs[i] == measles_ensemble_inc_arr[1:enddate, 1, i] "failed $i: $(measles_incs[i]) != $(measles_ensemble_inc_arr[1:enddate, 1, i])"
    println("i = $i: $(round(mean(measles_incs[i]); digits = 2)), $(round(mean(measles_ensemble_inc_arr[1:enddate, 1, i]); digits = 2))")
end
