#%%
using DrWatson
@quickactivate "OutbreakDetection"

using UnPack
using ColorSchemes

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec = IndividualTestSpecification(0.8, 0.8, 0)

#%%
for noise_specification in ensemble_noise_specification_vec
    scenario_specification = ScenarioSpecification(
        ensemble_specification,
        ensemble_outbreak_specification,
        noise_specification,
        OutbreakDetectionSpecification(5, 7, 0.6, 0.3),
        ensemble_single_individual_test_spec,
    )

    ensemble_solution_dict = get_ensemble_file(
        scenario_specification.ensemble_specification
    )

    incarr = get_ensemble_file(
        ensemble_specification, ensemble_outbreak_specification
    )

    ensemble_single_scenario_detection = get_ensemble_file(
        scenario_specification
    )

    noisearr, poisson_noise_prop = create_noise_arr(
        noise_specification,
        incarr["ensemble_inc_arr"];
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )

    noisedir = getdirpath(noise_specification)
    plotpath = joinpath(
        plotsdir(),
        "ensemble",
        "single-scenario",
        noisedir,
    )
    mkpath(plotpath)

    ensemble_noise_fig = visualize_ensemble_noise(
        noisearr,
        poisson_noise_prop,
        ensemble_time_specification,
        noisedir
    )

    save(
        joinpath(
            plotpath,
            "ensemble-sim_single-scenario_noise.png",
        ),
        ensemble_noise_fig; size = (2200, 1600),
    )
end
