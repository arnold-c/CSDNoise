#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes
using CSV: CSV
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
    10, 7, 0.6, 0.2, "inferred_movingavg"
)

ews_metric_specification = EWSMetricSpecification("backward", 35, 1)

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

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

ensemble_single_inc_ews = get_ensemble_file(
    ensemble_specification,
    ensemble_outbreak_specification,
    ews_metric_specification,
)["inc_ewsmetrics"]

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
for sim_num in [10, 20, 55]
    ensemble_single_scenario_Reffective_plot = Reff_plot(
        ensemble_single_incarr,
        ensemble_single_Reff_arr,
        ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs,
        ensemble_time_specification;
        sim = sim_num,
        threshold = ensemble_outbreak_specification.outbreak_threshold,
    )

    save(
        plotsdir(
            "ensemble/single-scenario/ensemble_single_scenario_Reffective_sim_$(sim_num).png",
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
    ews_metric_specification,
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
    ensemble_time_specification,
    ews_metric_specification,
)[[1, 3, 4]]

# testarr_lag, test_movingvg_arr_lag, inferred_positives_arr_lag = create_testing_arrs(
#     ensemble_single_incarr,
#     noisearr,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_single_individual_test_spec_lag,
#     ensemble_time_specification,
#     ews_metric_specification,
# )[[1, 3, 4]]

#%%
inferred_positive_no_lag_plot = lines(
    ensemble_time_specification.trange,
    ensemble_single_incarr[:, 1, 1];
    color = :orange,
)
lines!(
    ensemble_time_specification.trange,
    inferred_positives_arr_no_lag[:, 1];
    color = :black,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_no_lag.png",
    ),
    inferred_positive_no_lag_plot,
)

#%%
# inferred_positive_lag_plot = lines(
#     ensemble_time_specification.trange,
#     ensemble_single_incarr[:, 1, 1];
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
CSV.write(
    outdir("incidence-array_sim-$(sim_num).csv"),
    DataFrames.DataFrame(
        ensemble_single_incarr[:, :, sim_num],
        [:incidence, :above_threshold, :outbreak],
    ),
)

#%%
R"""
source(here::here("scripts","spaero-ews.R"))
"""

@rget spaero_ews_backward spaero_ews_centered

#%%
testinc = ensemble_single_incarr[:, 1, sim_num]
emptyvec = zeros(Float64, length(testinc))

#%%
centered_ews = SpaeroEWSMetrics(
    SpaeroEWSMetricSpecification(Centered, 30, 1),
    testinc,
    1.0
)
backward_ews = SpaeroEWSMetrics(
    SpaeroEWSMetricSpecification(Backward, 30, 1),
    testinc,
    1.0
)

#%%
compare_against_spaero(
    spaero_ews_centered,
    centered_ews;
    tolerance = 1e-10,
    showwarnings = false
)
compare_against_spaero(
    spaero_ews_backward,
    backward_ews;
    tolerance = 1e-10,
    showwarnings = false
)

#%%

# iod_df_centered = compare_against_spaero(
#     spaero_ews_centered.index_of_dispersion,
#     spaero_iod(Centered, testinc, 30)
# )
# filter_spaero_comparison(iod_df_centered)
#
# #%%
# skew_df_centered = compare_against_spaero(
#     spaero_ews_centered.skewness,
#     spaero_skew(Centered, testinc, 30)
# )
# filter_spaero_comparison(skew_df_centered)
#
# #%%
# kurtosis_df_centered = compare_against_spaero(
#     spaero_ews_centered.kurtosis,
#     spaero_kurtosis(Centered, testinc, 30)
# )
# filter_spaero_comparison(kurtosis_df_centered; tolerance = 1e-10)
#
# #%%
# autocov_df_centered = compare_against_spaero(
#     spaero_ews_centered.autocovariance,
#     spaero_autocov(Centered, testinc, 30)
# )
# filter_spaero_comparison(autocov_df_centered; tolerance = 1e-10)
#
# #%%
# autocor_df_centered = compare_against_spaero(
#     spaero_ews_centered.autocorrelation,
#     spaero_autocor(Centered, testinc, 30)
# )
# filter_spaero_comparison(autocor_df_centered; tolerance = 1e-10)
#
# #%%
# spaero_ensemble_single_inc_ews = StructArray([
#     EWSMetrics(
#         ensemble_time_specification.tstep,
#         ews_metric_specification,
#         spaero_ensemble_single_inc_ews_arr.mean,
#         spaero_ensemble_single_inc_ews_arr.variance,
#         spaero_ensemble_single_inc_ews_arr.coefficient_of_variation,
#         spaero_ensemble_single_inc_ews_arr.index_of_dispersion,
#         spaero_ensemble_single_inc_ews_arr.skewness,
#         spaero_ensemble_single_inc_ews_arr.kurtosis,
#         convert(
#             Vector{Float64},
#             replace!(
#                 spaero_ensemble_single_inc_ews_arr.autocovariance,
#                 missing => NaN,
#             ),
#         ),
#         convert(
#             Vector{Float64},
#             replace!(
#                 spaero_ensemble_single_inc_ews_arr.autocorrelation,
#                 missing => NaN,
#             ),
#         ),
#     ),
# ])
#
# #%%
# mean_df = DataFrames.DataFrame(
#     [
#         ensemble_single_inc_ews[1].mean,
#         spaero_ensemble_single_inc_ews[1].mean,
#     ],
#     [:mine, :spaero],
# )
#
# mean_df.diff = mean_df.mine .- mean_df.spaero;
# mean_df
#
# lines(mean_df.diff)
#
# #%%
# var_df = DataFrames.DataFrame(
#     [
#         ensemble_single_inc_ews[1].variance,
#         spaero_ensemble_single_inc_ews[1].variance,
#     ],
#     [:mine, :spaero],
# )
#
# var_df.diff = var_df.mine .- var_df.spaero;
# var_df
#
# lines(var_df.diff)
#
# #%%
# ewsmetrics = [
#     :autocorrelation,
#     :autocovariance,
#     :coefficient_of_variation,
#     :index_of_dispersion,
#     :kurtosis,
#     :mean,
#     :skewness,
#     :variance,
# ]
#
# basedir = plotsdir(
#     "ensemble/single-scenario/ewsmetrics/$(ews_metric_specification.dirpath)"
# )
# mkpath(basedir)
# testdir = joinpath(basedir, "test")
# mkpath(testdir)
# incdir = joinpath(basedir, "incidence")
# spaero_comparison_dir = joinpath(incdir, "spaero-comparison")
# mkpath(spaero_comparison_dir)
#
# for ewsmetric in ewsmetrics
#     inc_ews_metric_plot = Reff_ews_plot(
#         ensemble_single_incarr,
#         ensemble_single_Reff_arr,
#         ensemble_single_Reff_thresholds_vec,
#         ensemble_single_inc_ews,
#         ewsmetric,
#         ensemble_single_periodsum_vecs,
#         ensemble_time_specification;
#         sim = sim_num,
#         threshold = ensemble_outbreak_specification.outbreak_threshold,
#         plottitle = "Incidence $(ewsmetric): $(ews_metric_specification.method)",
#     )
#
#     save(
#         joinpath(
#             incdir,
#             "ensemble_single_scenario_incidence_metric_$(String(ewsmetric)).png",
#         ),
#         inc_ews_metric_plot,
#     )
#
#     spaero_comparison_inc_ews_metric_plot = Reff_ews_plot(
#         ensemble_single_incarr,
#         ensemble_single_Reff_arr,
#         ensemble_single_Reff_thresholds_vec,
#         ensemble_single_inc_ews,
#         spaero_ensemble_single_inc_ews,
#         ewsmetric,
#         ensemble_single_periodsum_vecs,
#         ensemble_time_specification;
#         sim = sim_num,
#         threshold = ensemble_outbreak_specification.outbreak_threshold,
#         plottitle = "Incidence $(ewsmetric): $(ews_metric_specification.method)",
#     )
#
#     save(
#         joinpath(
#             spaero_comparison_dir,
#             "ensemble_single_scenario_incidence_metric_$(String(ewsmetric)).png",
#         ),
#         spaero_comparison_inc_ews_metric_plot,
#     )
#
#     for (ewsmetric_sa, lag_label) in
#         zip((
#         # ewsmetrics_lag,
#         [ewsmetrics_no_lag]
#     ), (
#         # "lag",
#         ["no_lag"]
#     ))
#         test_ews_metric_plot = Reff_ews_plot(
#             ensemble_single_incarr,
#             ensemble_single_Reff_arr,
#             ensemble_single_Reff_thresholds_vec,
#             ewsmetric_sa,
#             ewsmetric,
#             ensemble_single_periodsum_vecs,
#             ensemble_time_specification;
#             sim = sim_num,
#             threshold = ensemble_outbreak_specification.outbreak_threshold,
#             plottitle = "Test Positive $(ewsmetric)\t$(lag_label): $(ews_metric_specification.method)",
#         )
#
#         save(
#             joinpath(
#                 testdir,
#                 "ensemble_single_scenario_testpositive_metric_$(String(ewsmetric))_$lag_label.png",
#             ),
#             test_ews_metric_plot,
#         )
#     end
# end
#
# # plot_all_single_scenarios(
# #     noisearr,
# #     poisson_noise_prop,
# #     noisedir,
# #     OT_chars,
# #     ensemble_single_incarr,
# #     testarr,
# #     test_movingvg_arr,
# #     ensemble_single_individual_test_spec_no_lag,
# #     ensemble_single_outbreak_detection_spec,
# #     ensemble_time_specification,
# # )
#
# # GC.gc(true)
# # @info "Finished plotting the single scenario for $(noisedir)"
# # println("=================================================================")
