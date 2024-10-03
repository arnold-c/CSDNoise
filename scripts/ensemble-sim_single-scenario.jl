#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes

using CSDNoise
using StructArrays

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec = IndividualTestSpecification(
    1.0, 1.0, 0
)
ensemble_single_outbreak_detection_spec = OutbreakDetectionSpecification(
    10, 7, 1.0, 1.0, "movingavg"
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

ensemble_single_Reff_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_arr"]

ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

#%%
nsims_plot = 5
ensemble_incarr_Reff_plot(
    ensemble_single_incarr[:, :, 1:nsims_plot],
    ensemble_single_Reff_arr[:, 1:nsims_plot],
    ensemble_single_Reff_thresholds_vec[1:nsims_plot],
    ensemble_single_periodsum_vecs[1:nsims_plot], ;
    outbreak_alpha = 0.1,
    Reff_alpha = 1,
)

#%%
ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_incarr,
    ensemble_single_seir_arr,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    threshold = ensemble_outbreak_specification.outbreak_threshold,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_incidence_prevalence.png"
    ),
    ensemble_single_scenario_incidence_prevalence_plot,
)

#%%
ews_method_vec = [Centered]
ews_aggregation_vec = [7]
ews_bandwidth_vec = [35]
ews_lag_vec = [1]

ews_spec_vec = create_combinations_vec(
    EWSMetricSpecification,
    (ews_method_vec, ews_aggregation_vec, ews_bandwidth_vec, ews_lag_vec),
)

# for (noise_specification, ews_metric_specification) in Iterators.product(
#     ensemble_noise_specification_vec,
#     ews_spec_vec,
# )
noise_specification = ensemble_noise_specification_vec[1]
ews_metric_specification = ews_spec_vec[1]

scenario_specification = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    noise_specification,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec,
    ews_metric_specification,
)

noisearr, poisson_noise_prop = create_noise_arr(
    noise_specification,
    ensemble_single_incarr;
    ensemble_specification = ensemble_specification,
    seed = 1234,
)
noisedir = getdirpath(noise_specification)

testarr, ewsvec, test_movingavg_arr, inferred_positives_arr = create_testing_arrs(
    ensemble_single_incarr,
    noisearr,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec,
    ensemble_time_specification,
    ews_metric_specification;
)

plot_all_single_scenarios(
    noisearr,
    poisson_noise_prop,
    noisedir,
    ensemble_single_incarr,
    testarr,
    test_movingavg_arr,
    ensemble_single_individual_test_spec,
    ensemble_single_outbreak_detection_spec,
    ensemble_time_specification,
)

#%%
GC.gc(true)
@info "Finished plotting the single scenario for $(noisedir)"
println("=================================================================")
# end

#%%
Reff_ews_plot(
    ensemble_single_incarr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    StructArray(ewsvec),
    :variance,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    aggregation = 7,
)
