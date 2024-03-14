#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes

using CSDNoise

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec_no_lag = IndividualTestSpecification(
    1.0, 1.0, 0
)

ensemble_single_individual_test_spec_lag = IndividualTestSpecification(
    1.0, 1.0, 14
)

ensemble_single_outbreak_detection_spec = OutbreakDetectionSpecification(
    10, 7, 0.6, 0.2, "inferred_movingavg"
)

#%%
ensemble_single_seir_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_seir_arr"]

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

#%%
mkpath(plotsdir("ensemble/single-scenario"))

ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_incarr,
    ensemble_single_seir_arr,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    threshold = ensemble_outbreak_specification.outbreak_threshold,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_incidence_prevalence.png",
    ),
    ensemble_single_scenario_incidence_prevalence_plot,
)

#%%
noise_specification = ensemble_noise_specification_vec[1]

scenario_specification_no_lag = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    noise_specification,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_no_lag,
)

scenario_specification_lag = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    noise_specification,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_lag,
)

ensemble_solution_dict_no_lag = get_ensemble_file(scenario_specification_no_lag)
ensemble_solution_dict_lag = get_ensemble_file(scenario_specification_lag)

OT_chars_no_lag = ensemble_solution_dict_no_lag["OT_chars"]
OT_chars_lag = ensemble_solution_dict_lag["OT_chars"]

noisearr, poisson_noise_prop = create_noise_arr(
    noise_specification,
    ensemble_single_incarr;
    ensemble_specification = ensemble_specification,
    seed = 1234,
)
noisedir = getdirpath(noise_specification)

testarr_no_lag, test_movingvg_arr_no_lag, inferred_positives_arr_no_lag = create_testing_arrs(
    ensemble_single_incarr,
    noisearr,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_no_lag,
)

testarr_lag, test_movingvg_arr_lag, inferred_positives_arr_lag = create_testing_arrs(
    ensemble_single_incarr,
    noisearr,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_lag,
)

inferred_positive_no_lag_plot = lines(ensemble_time_specification.trange, ensemble_single_incarr[:, 1, 1]; color = :orange)
lines!(ensemble_time_specification.trange, inferred_positives_arr_no_lag[:, 1]; color = :black)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_no_lag.png",
    ),
    inferred_positive_no_lag_plot,
)

inferred_positive_lag_plot = lines(ensemble_time_specification.trange, ensemble_single_incarr[:, 1, 1]; color = :orange)
lines!(ensemble_time_specification.trange, inferred_positives_arr_lag[:, 1]; color = :black)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_lag.png",
    ),
    inferred_positive_lag_plot,
)

compare_inferred_positives_plot = lines(ensemble_time_specification.trange, inferred_positives_arr_no_lag[:, 1]; color = :orange)
lines!(ensemble_time_specification.trange, inferred_positives_arr_lag[:, 1]; color = :black)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_compare_inferred_positives.png",
    ),
    compare_inferred_positives_plot,
)

# plot_all_single_scenarios(
#     noisearr,
#     poisson_noise_prop,
#     noisedir,
#     OT_chars,
#     ensemble_single_incarr,
#     testarr,
#     test_movingvg_arr,
#     ensemble_single_individual_test_spec_no_lag,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_time_specification,
# )

GC.gc(true)
@info "Finished plotting the single scenario for $(noisedir)"
println("=================================================================")
