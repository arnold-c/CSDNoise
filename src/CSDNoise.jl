"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module CSDNoise

# using Reexport
include("DrWatson-helpers.jl")
export outdir

include("transmission-functions.jl")
export calculate_beta, calculate_beta_amp,
    calculateReffective, calculateReffective_t!, calculateR0,
    calculate_import_rate, calculate_mu
# @reexport using .TransmissionFunctions

include("dynamics-constants.jl")
export POPULATION_N, LATENT_PER_DAYS, DUR_INF_DAYS, R0, SIGMA, GAMMA,
    LIFE_EXPECTANCY_YEARS, ANNUAL_BIRTHS_PER_K, VACCINATION_COVERAGE,
    MU, BETA_MEAN, BETA_FORCE, EPSILON

# @reexport using .EWSMetrics

include("structs.jl")
export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
    StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
    OutbreakSpecification, IndividualTestSpecification, get_test_description,
    PoissonNoiseSpecification, DynamicalNoiseSpecification, NoiseSpecification,
    get_noise_description, get_noise_magnitude, get_noise_magnitude_description,
    getdirpath,
    EWSMethod, Backward, Centered,
    EWSMetricSpecification,
    ScenarioSpecification, TestPositivity, OptimalThresholdCharacteristics,
    EWSMetrics
# @reexport using .ODStructs

include("ews-metrics.jl")
export aggregate_timeseries,
    spaero_mean, spaero_mean!,
    spaero_var,
    spaero_cov,
    spaero_iod,
    spaero_skew,
    spaero_kurtosis,
    spaero_autocov,
    spaero_autocor,
    compare_against_spaero, filter_spaero_comparison,
    EWSMetricSpecification, EWSMetrics

include("ews-functions.jl")
export EWSThresholdWindow,
    Rolling, Expanding, expanding_ews_thresholds,
    tycho_testing_plots,
    calculate_ews_lead_time, calculate_ews_trigger_index

include("tycho-cleaning.jl")
export cdc_week_to_date,
    calculate_aggregation_cases, fill_aggregation_values,
    calculate_ews_enddate

include("test-constants.jl")
export CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC, CLINICAL_TEST_SPECS

include("SEIR-model.jl")
export seir_mod, seir_mod!, seir_mod_loop!,
    convert_svec_to_matrix, convert_svec_to_matrix!, convert_svec_to_array
# @reexport using .SEIRModel

include("cleaning-functions.jl")
export create_sir_df
# @reexport using .CleaningFunctions

include("detection-thresholds.jl")
export create_inc_infec_arr,
    create_inc_infec_arr!,
    calculate_outbreak_thresholds,
    Reff_ge_than_one
# @reexport using .DetectionThresholds

include("diag-testing-functions.jl")
export create_testing_arrs, create_testing_arrs!, calculate_tested!,
    calculate_positives!,
    calculate_true_positives!, calculate_noise_positives!,
    infer_true_positives, calculate_test_positivity_rate,
    calculate_movingavg, calculate_movingavg!,
    detectoutbreak, detectoutbreak!, calculate_daily_detection_characteristics,
    calculate_noutbreaks, calculate_OutbreakThresholdChars,
    calculate_test_positivity, calculate_outbreak_detection_characteristics,
    filter_first_matched_bounds, calculate_first_matched_bounds_index,
    calculate_cases_before_after_alert!, calculate_cases_before_after_alert
# @reexport using .DiagTestingFunctions

include("ensemble-functions.jl")
export create_combinations_vec, create_ensemble_spec_combinations,
    run_ensemble_jump_prob, run_jump_prob,
    run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation,
    get_ensemble_file
# @reexport using .EnsembleFunctions

include("noise-functions.jl")
export create_noise_arr, add_poisson_noise_arr!
# @reexport using .NoiseFunctions

include(
    "plotting-functions.jl"
)
export seircolors,
    seir_state_labels, create_sir_plot, draw_sir_plot,
    ACCURACY_COLOR, DAILY_SENSITIVITY_COLOR, DAILY_SPECIFICITY_COLOR,
    DAILY_PPV_COLOR, DAILY_NPV_COLOR,
    DETECTION_DELAY_COLOR, N_ALERTS_PER_OUTBREAK_COLOR,
    N_FALSE_ALERTS_COLOR, N_ALERTS_COLOR, N_OUTBREAKS_COLOR,
    N_MISSED_OUTBREAKS_COLOR, PERC_OUTBREAKS_DETECTED_COLOR,
    PERC_OUTBREAKS_MISSED_COLOR, PERC_ALERTS_CORRECT_COLOR,
    PERC_ALERTS_FALSE_COLOR,
    sir_quantiles_array_base_plot, create_sir_quantiles_plot,
    incidence_prevalence_plot, Reff_plot,
    visualize_ensemble_noise, incidence_testing_plot,
    testing_plot, ensemble_outbreak_distribution_plot, ensemble_OTChars_plot,
    singlescenario_test_positivity_plot, test_positivity_distribution_plot,
    ensemble_outbreak_detect_diff_plot, save_compare_ensemble_OTchars_plot,
    compare_ensemble_OTchars_plots,
    compare_optimal_thresholds_chars_plot, create_optimal_thresholds_chars_plot,
    compare_optimal_thresholds_test_chars_plot,
    create_optimal_thresholds_test_chars_plot,
    Reff_ews_plot, tycho_epicurve, tycho_noise_components_epicurve,
    tycho_test_positive_components_epicurve,
    tycho_tau_distribution
# @reexport using .PlottingFunctions

include("ensemble-sim_single-scenario_plots.jl")
export plot_all_single_scenarios

include("threshold_comparison_plots.jl")
export collect_threshold_char_vec, plot_all_threshold_comparisons

include("optimal-threshold-functions.jl")
export calculate_optimal_threshold, calculate_OptimalThresholdCharacteristics,
    calculate_optimal_threshold_summaries,
    create_optimal_thresholds_df, create_wide_optimal_thresholds_df,
    create_and_save_xlsx_optimal_threshold_summaries,
    create_optimal_threshold_summary_df,
    create_wide_optimal_threshold_summary_df,
    create_all_wide_optimal_threshold_summary_dfs,
    save_xlsx_optimal_threshold_summaries,
    create_and_save_xlsx_optimal_threshold_summaries

@static if false
    include("scripts/ensemble-sim_inferred-scenario-visualizations.jl")
    include("scripts/tycho-visualizations.jl")
end

end
