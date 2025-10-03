#%%
using CSDNoise
using Random
using FLoops
using JLD2
using Dates
using ProgressMeter

#%%
rubella_dynamical_noise_spec = DynamicalNoiseSpecification(
    R0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 0.15,
)

#%%
vaccination_levels = 0.0:0.01:1.0

nsims = 10_000
burnin_years = 0
tmax_years = 100

ensemble_spec = create_ensemble_specs(
    EnsembleSpecsParameters(
        burnin_years = burnin_years,
        tmax_years = tmax_years,
        annual_births_per_k = ANNUAL_BIRTHS_PER_K,
        ensemble_state_specification = StateParameters(
            500_000,
            Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
        ),
        R_0 = rubella_dynamical_noise_spec.R0,
        gamma = 1 / rubella_dynamical_noise_spec.duration_infection,
        sigma = 1 / rubella_dynamical_noise_spec.latent_period,
        target_Reff = 0.9,
        target_years = 10,
        min_vaccination_coverage = 0.0,
        max_vaccination_coverage = 1.0,
        nsims = nsims
    )
)[1]

enddates_vec = fill(ensemble_spec.time_parameters.tlength, nsims)

#%%
dyn = ensemble_spec.dynamics_parameter_specification
@unpack mu, gamma, R_0, beta_mean = dyn
@unpack N = ensemble_spec.state_parameters.init_states

((mu + gamma) * R_0) / (beta_mean * N) - 1

((mu * N * R_0) / beta_mean) * (1 - 0.0 - (1 / R_0))

(mu / beta_mean) * (R_0 * N * (1 - 0.1) - 1)

calculate_endemic_equilibrium_proportions(dyn, 0.0).e_prop * N

mu * (R_0 * (1 - 0.0) - 1) / beta_mean


#%%
function run_vaccination_sweep(
        vaccination_level,
        dynamical_noise_spec,
        ensemble_spec,
        enddates_vec,
        nsims,
        tmax_years,
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")
    )

    Random.seed!(1234)

    noise_run = recreate_noise_vecs(
        rubella_dynamical_noise_spec,
        vaccination_level,
        1 - vaccination_level,
        ensemble_spec,
        enddates_vec
    )

    results = (
        vaccination_level = vaccination_level,
        mean_noise = noise_run.mean_noise,
        mean_poisson_noise = noise_run.mean_poisson_noise,
        mean_dynamic_noise = noise_run.mean_dynamic_noise,
        incidence = noise_run.incidence,
    )

    output_file = joinpath(
        "out",
        "rubella_vaccination_sweep_results_vaccination-$(vaccination_level)_$(timestamp).jld2"
    )

    JLD2.jldsave(
        output_file;
        results = results,
        vaccination_level = vaccination_level,
        noise_spec = rubella_dynamical_noise_spec,
        nsims = nsims,
        tmax_years = tmax_years,
        timestamp = timestamp
    )

    println("Results saved to $output_file")
    return nothing
end

#%%
timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HHMMSS")

@showprogress for  vacc_level in vaccination_levels

    run_vaccination_sweep(
        vacc_level,
        rubella_dynamical_noise_spec,
        ensemble_spec,
        enddates_vec,
        nsims,
        tmax_years,
        timestamp
    )

end

#%%
vac_0_perc_res = JLD2.load(outdir("rubella_vaccination_sweep_results_vaccination-0.0_2025-10-02_230723.jld2"))

#%%
final_inc = StatsBase.mean(
    map(
        inc -> inc[end],
        vac_0_perc_res["results"].incidence
    )
)

final_inc = StatsBase.mean(
    map(
        inc -> inc[end],
        vac_2_perc_res["results"].incidence
    )
)


#%%
vac_2_perc_res = JLD2.load(outdir("rubella_vaccination_sweep_results_vaccination-0.02_2025-10-02_230723.jld2"))

#%%
inc = vac_2_perc_res["results"].incidence

#%%
ensemble_spec

#%%
fig = Figure()
ax = Axis(fig[1, 1])
for i in 1:1000
    lines!(ax, @view(inc[i][(end - 3000):end]); color = (:blue, 0.1))
end
display(fig)
