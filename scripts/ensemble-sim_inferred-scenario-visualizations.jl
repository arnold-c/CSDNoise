#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes
using CSV: CSV
using DelimitedFiles: DelimitedFiles
using DataFrames: DataFrames
using RCall
using StructArrays
using SumTypes

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
    10, 7, 1.0, 1.0, "inferred_movingavg"
)

ews_metric_specification_1d = EWSMetricSpecification(Centered, 1, 35, 1)
ews_metric_specification_30d = EWSMetricSpecification(Centered, 30, 35, 1)

#%%
ensemble_single_seir_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_seir_arr"]

ensemble_single_Reff_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_arr"]

ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

emergent_incidence_arr = ensemble_single_scenario_inc_file["emergent_incidence_arr"]
emergent_outbreak_threshold_vec = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

ensemble_single_inc_ews_1d = get_ensemble_file(
    ensemble_specification,
    ensemble_outbreak_specification,
    ews_metric_specification_1d,
)["inc_ewsmetrics"]

ensemble_single_inc_ews_30d = get_ensemble_file(
    ensemble_specification,
    ensemble_outbreak_specification,
    ews_metric_specification_30d,
)["inc_ewsmetrics"]

#%%
mkpath(plotsdir("ensemble/single-scenario"))

ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    emergent_incidence_arr,
    ensemble_single_seir_arr,
    emergent_outbreak_threshold_vec,
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
for sim_num in [10, 20, 55]
    ensemble_single_scenario_Reffective_plot = Reff_plot(
        emergent_incidence_arr,
        ensemble_single_Reff_arr,
        ensemble_single_Reff_thresholds_vec,
        emergent_outbreak_threshold_vec,
        ensemble_time_specification;
        sim = sim_num,
        threshold = ensemble_outbreak_specification.outbreak_threshold,
    )

    save(
        plotsdir(
            "ensemble/single-scenario/ensemble_single_scenario_Reffective_sim_$(sim_num).png"
        ),
        ensemble_single_scenario_Reffective_plot,
    )
end

#%%
noise_specification = ensemble_noise_specification_vec[1]

scenario_specification_no_lag = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    noise_specification,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_no_lag,
    ews_metric_specification_1d,
)

# scenario_specification_lag = ScenarioSpecification(
#     ensemble_specification,
#     ensemble_outbreak_specification,
#     noise_specification,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_single_individual_test_spec_lag,
#     ews_metric_specification,
# )

ensemble_solution_dict_no_lag = get_ensemble_file(scenario_specification_no_lag)
# ensemble_solution_dict_lag = get_ensemble_file(scenario_specification_lag)

OT_chars_no_lag = ensemble_solution_dict_no_lag["OT_chars"]
# OT_chars_lag = ensemble_solution_dict_lag["OT_chars"]

ewsmetrics_no_lag = ensemble_solution_dict_no_lag["test_ewsmetrics"]
# ewsmetrics_lag = ensemble_solution_dict_lag["test_ewsmetrics"]

noisearr, poisson_noise_prop = create_noise_arr(
    noise_specification,
    emergent_incidence_arr;
    ensemble_specification = ensemble_specification,
    seed = 1234,
)
noisedir = getdirpath(noise_specification)

testarr_no_lag, test_movingvg_arr_no_lag, inferred_positives_arr_no_lag = create_testing_arrs(
    emergent_incidence_arr,
    noisearr,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_no_lag,
    ensemble_time_specification,
    ews_metric_specification_1d,
)[[1, 3, 4]]

# testarr_lag, test_movingvg_arr_lag, inferred_positives_arr_lag = create_testing_arrs(
#     emergent_incidence_arr,
#     noisearr,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_single_individual_test_spec_lag,
#     ensemble_time_specification,
#     ews_metric_specification,
# )[[1, 3, 4]]

#%%
inferred_positive_no_lag_plot = lines(
    ensemble_time_specification.trange,
    emergent_incidence_arr[:, 1, 1];
    color = :orange,
)
lines!(
    ensemble_time_specification.trange,
    inferred_positives_arr_no_lag[:, 1];
    color = :black,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_no_lag.png"
    ),
    inferred_positive_no_lag_plot,
)

#%%
# inferred_positive_lag_plot = lines(
#     ensemble_time_specification.trange,
#     emergent_incidence_arr[:, 1, 1];
#     color = :orange,
# )
# lines!(
#     ensemble_time_specification.trange,
#     inferred_positives_arr_lag[:, 1];
#     color = :black,
# )
#
# save(
#     plotsdir(
#         "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_lag.png",
#     ),
#     inferred_positive_lag_plot,
# )

#%%
# compare_inferred_positives_plot = lines(
#     ensemble_time_specification.trange,
#     inferred_positives_arr_no_lag[:, 1];
#     color = :orange,
# )
# lines!(
#     ensemble_time_specification.trange,
#     inferred_positives_arr_lag[:, 1];
#     color = :black,
# )
#
# save(
#     plotsdir(
#         "ensemble/single-scenario/ensemble_single_scenario_compare_inferred_positives.png",
#     ),
#     compare_inferred_positives_plot,
# )

#%%
sim_num = 1

#%%
# Can use RCall to @rput, but write to csv so can put R code in own file
# that can be run independently (as well as provide linting etc)
DelimitedFiles.writedlm(
    outdir(
        "tycho",
        "tests",
        "incidence-array_1d_sim-$(sim_num).csv",
    ),
    emergent_incidence_arr[:, 1, sim_num],
    ',',
)

DelimitedFiles.writedlm(
    outdir(
        "tycho",
        "tests",
        "incidence-array_30d_sim-$(sim_num).csv",
    ),
    CSDNoise.aggregate_timeseries(emergent_incidence_arr[:, 1, sim_num], 30),
    ',',
)

#%%
R"""
source(here::here("scripts", "spaero-ews.R"))
"""

@rget spaero_ews_backward_1d spaero_ews_backward_30d spaero_ews_centered_1d spaero_ews_centered_30d

#%%
backward_ews_1d = StructArray(
    EWSMetrics[
        EWSMetrics(
                EWSMetricSpecification(Backward, 1, 35, 1),
                emergent_incidence_arr[:, 1, sim],
            ) for sim in axes(emergent_incidence_arr, 3)
    ],
)

backward_ews_30d = StructArray(
    EWSMetrics[
        EWSMetrics(
                EWSMetricSpecification(Backward, 30, 35, 1),
                emergent_incidence_arr[:, 1, sim],
            ) for sim in axes(emergent_incidence_arr, 3)
    ],
)

centered_ews_1d = StructArray(
    EWSMetrics[
        EWSMetrics(
                EWSMetricSpecification(Centered, 1, 35, 1),
                emergent_incidence_arr[:, 1, sim],
            ) for sim in axes(emergent_incidence_arr, 3)
    ],
)

centered_ews_30d = StructArray(
    EWSMetrics[
        EWSMetrics(
                EWSMetricSpecification(Centered, 30, 35, 1),
                emergent_incidence_arr[:, 1, sim],
            ) for sim in axes(emergent_incidence_arr, 3)
    ],
)

#%%
ewsmetrics = [
    :autocorrelation,
    :autocovariance,
    :coefficient_of_variation,
    :index_of_dispersion,
    :kurtosis,
    :mean,
    :skewness,
    :variance,
]

#%%
compare_against_spaero(
    spaero_ews_centered_1d,
    ensemble_single_inc_ews_1d[sim_num];
    tolerance = 1.0e-10,
    showwarnings = false,
    ews = [:mean],
)

compare_against_spaero(
    spaero_ews_backward_1d,
    backward_ews_1d[sim_num];
    tolerance = 1.0e-10,
    showdiffs = true,
    ews = ewsmetrics,
)

compare_against_spaero(
    spaero_ews_backward_30d,
    backward_ews_30d[sim_num];
    tolerance = 1.0e-10,
    showdiffs = true,
    ews = ewsmetrics,
)

#%%
Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    backward_ews_1d,
    spaero_ews_backward_1d,
    :autocorrelation,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification
)

#%%
autocov_backward_1d = Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    backward_ews_1d,
    spaero_ews_backward_1d,
    :autocovariance,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    plottitle = "Backward 1d Autocovariance",
    ylims_metric = (0, 0.5),
)

save(
    plotsdir(
        "ensemble/single-scenario/ewsmetrics/ensemble_single_scenario_autocovar_backward_1d.png"
    ),
    autocov_backward_1d,
)
#%%

autocov_centered_1d = Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    centered_ews_1d,
    spaero_ews_centered_1d,
    :autocovariance,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    plottitle = "centered 1d Autocovariance",
    ylims_metric = (0, 0.5),
)

save(
    plotsdir(
        "ensemble/single-scenario/ewsmetrics/ensemble_single_scenario_autocovar_centered_1d.png"
    ),
    autocov_centered_1d,
)

#%%
Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    backward_ews_30d,
    spaero_ews_backward_30d,
    :autocorrelation,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    aggregation = 30,
)

#%%
Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    backward_ews_30d,
    spaero_ews_backward_30d,
    :autocovariance,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    aggregation = 30,
)

#%%
var_backward_30d = Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    backward_ews_30d,
    spaero_ews_backward_30d,
    :autocovariance,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    aggregation = 30,
    plottitle = "Backward 30d Autocovariance",
    ylims_metric = (0, 100),
)

save(
    plotsdir(
        "ensemble/single-scenario/ewsmetrics/ensemble_single_scenario_autocovar_backward_30d.png"
    ),
    var_backward_30d,
)

#%%
var_centered_30d = Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    centered_ews_30d,
    spaero_ews_centered_30d,
    :autocovariance,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    aggregation = 30,
    plottitle = "Centered 30d Autocovariance",
)

save(
    plotsdir(
        "ensemble/single-scenario/ewsmetrics/ensemble_single_scenario_autocovar_centered_30d.png"
    ),
    var_centered_30d,
)

#%%
Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    centered_ews_30d,
    spaero_ews_centered_30d,
    :autocovariance,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    aggregation = 30,
)

#%%
Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    centered_ews_30d,
    spaero_ews_centered_30d,
    :coefficient_of_variation,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification;
    aggregation = 30,
)

#%%
Reff_ews_plot(
    emergent_incidence_arr,
    ensemble_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    backward_ews_1d,
    spaero_ews_backward_1d,
    :autocorrelation,
    emergent_outbreak_threshold_vec,
    ensemble_time_specification
)

#%%
for ewsmetric in ewsmetrics
    for (ews_metric_specification, inc_ews) in zip(
            [ews_metric_specification_1d, ews_metric_specification_30d],
            [ensemble_single_inc_ews_1d, ensemble_single_inc_ews_30d],
        )
        basedir = plotsdir(
            "ensemble/single-scenario/ewsmetrics/$(ews_metric_specification.dirpath)"
        )
        mkpath(basedir)
        testdir = joinpath(basedir, "test")
        mkpath(testdir)
        incdir = joinpath(basedir, "incidence")
        spaero_comparison_dir = joinpath(incdir, "spaero-comparison")
        mkpath(spaero_comparison_dir)

        inc_ews_metric_plot = Reff_ews_plot(
            emergent_incidence_arr,
            ensemble_single_Reff_arr,
            ensemble_single_Reff_thresholds_vec,
            inc_ews,
            ewsmetric,
            emergent_outbreak_threshold_vec,
            ensemble_time_specification;
            sim = sim_num,
            threshold = ensemble_outbreak_specification.outbreak_threshold,
            plottitle = "Incidence $(ewsmetric): $(ews_metric_specification.method)",
        )

        save(
            joinpath(
                incdir,
                "ensemble_single_scenario_incidence_metric_$(String(ewsmetric)).png",
            ),
            inc_ews_metric_plot,
        )
    end

    for (ewsmetric_sa, lag_label) in
        zip(
            (
                # ewsmetrics_lag,
                [ewsmetrics_no_lag]
            ), (
                # "lag",
                ["no_lag"]
            )
        )
        test_ews_metric_plot = Reff_ews_plot(
            emergent_incidence_arr,
            ensemble_single_Reff_arr,
            ensemble_single_Reff_thresholds_vec,
            ewsmetric_sa,
            ewsmetric,
            emergent_outbreak_threshold_vec,
            ensemble_time_specification;
            sim = sim_num,
            threshold = ensemble_outbreak_specification.outbreak_threshold,
            plottitle = "Test Positive $(ewsmetric)\t$(lag_label): $(ews_metric_specification.method)",
        )

        save(
            joinpath(
                testdir,
                "ensemble_single_scenario_testpositive_metric_$(String(ewsmetric))_$lag_label.png",
            ),
            test_ews_metric_plot,
        )
    end
end

# plot_all_single_scenarios(
#     noisearr,
#     poisson_noise_prop,
#     noisedir,
#     OT_chars,
#     emergent_incidence_arr,
#     testarr,
#     test_movingvg_arr,
#     ensemble_single_individual_test_spec_no_lag,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_time_specification,
# )

# GC.gc(true)
# @info "Finished plotting the single scenario for $(noisedir)"
# println("=================================================================")
