"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module CSDNoise
using DispatchDoctor
# @stable default_mode = "warn" default_union_limit = 2 begin

include("helpers.jl")

# Constants
include("./constants/dynamics-constants.jl")

# Types
include("./types/time-parameters.jl")
include("./types/dynamics-parameters.jl")
include("./types/state-parameters.jl")
include("./types/ensemble-specifications.jl")
include("./types/test-specifications.jl")
include("./types/outbreak-specifications.jl")
include("./types/noise-specifications.jl")
include("./types/ews-specifications.jl")
include("./types/scenario-specifications.jl")
include("./types/simulation-results.jl")
include("./types/optimization-specifications.jl")

# Type descriptions
include("./vaccination-rate-calculations.jl")
include("./test-description.jl")
include("./noise-description.jl")

include("./transmission-functions.jl")

# Constants that rely on transmisison function definitions
include("./constants/dynamics-constants_calculated.jl")


# EWS Metric functions
include("./ews-metrics/ews-alerts.jl")
include("./ews-metrics/ews-bandwidths.jl")
include("./ews-metrics/ews-enddates.jl")
include("./ews-metrics/ews-lead-time.jl")
include("./ews-metrics/ews-metrics.jl")
include("./ews-metrics/ews-summaries.jl")
include("./ews-metrics/ews-thresholds.jl")
include("./ews-metrics/ews-timeseries-aggregation.jl")
include("./ews-metrics/ews-to-df.jl")


include("optimization-setup-functions.jl")
include("optimization-checkpointing-functions.jl")
include("optimization-functions.jl")

include("ews-hyperparam-optimization.jl")
export
    load_most_recent_hyperparam_file,
    get_most_recent_hyperparam_filepath,
    optimal_ews_heatmap_df

include("ews-multistart-optimization.jl")
export ews_multistart_optimization,
    optimize_single_scenario,
    ews_objective_function_with_tracking,
    map_continuous_to_ews_parameters,
    create_optimization_scenarios,
    OptimizationTracker,
    OptimizedValues,
    confirm_optimization_run_structvector,
    optimize_scenarios_in_batches_structvector,
    load_previous_multistart_results_structvector,
    create_scenarios_structvector

include("ews-hyperparam-gridsearch-structvector.jl")
export ews_hyperparam_gridsearch_structvector,
    create_gridsearch_scenarios_structvector,
    evaluate_single_gridsearch_scenario,
    evaluate_gridsearch_scenarios,
    evaluate_gridsearch_scenarios2,
    evaluate_gridsearch_scenarios_optimized,
    evaluate_gridsearch_scenarios_bumper,
    evaluate_gridsearch_scenarios_multistart,
    optimize_quantile_threshold_multistart,
    quantile_only_objective_function_with_tracking,
    find_missing_scenarios,
    load_previous_gridsearch_results_structvector,
    confirm_gridsearch_run_structvector

include("ews-survival.jl")
export simulate_and_plot_ews_survival,
    simulate_ews_survival_data,
    create_ews_survival_data

include("test-constants.jl")
export CLINICAL_CASE_TEST_SPEC,
    EPI_LINKED_CASE_TEST_SPEC,
    CLINICAL_TEST_SPECS

include("SEIR-model.jl")
export seir_mod,
    seir_mod!,
    seir_mod_loop!,
    convert_svec_to_matrix,
    convert_svec_to_matrix!,
    convert_svec_to_array

include("detection-thresholds.jl")
export calculate_outbreak_thresholds,
    calculate_outbreak_thresholds!,
    calculate_above_threshold_bounds,
    Reff_ge_than_one

include("diag-testing-functions.jl")
export create_testing_arrs,
    create_testing_arrs!,
    calculate_tested!,
    calculate_positives!,
    calculate_true_positives!,
    calculate_noise_positives!,
    infer_true_positives,
    calculate_test_positivity_rate,
    calculate_movingavg,
    calculate_movingavg!,
    calculate_test_positivity

include("diag-testing-structvector-functions.jl")
export create_testing_vecs,
    calculate_tested_vec!,
    calculate_positives_vec!

include("ensemble-functions.jl")
export create_combinations_vec,
    create_ensemble_spec_combinations,
    run_ensemble_jump_prob,
    run_jump_prob,
    get_ensemble_file

include("noise-functions.jl")
export create_noise_vecs

include("filter-seir-results.jl")

include("calculate-noise-vacc-level-functions.jl")
export calculate_dynamic_vaccination_coverage,
    calculate_mean_dynamical_noise,
    calculate_dynamic_vaccination_coverage,
    calculate_ews_endpoints,
    calculate_mean_incidence,
    calculate_dynamic_vaccination_coverage_multistart_with_endpoints,
    calculate_mean_dynamical_noise_variable_length

# Plotting Functions
include("./plotting-functions/helpers_plots.jl")
include("./plotting-functions/tau-heatmap_plots.jl")
include("./plotting-functions/lead-time_plots.jl")


include("plotting-functions/single-simulation_plots.jl")
export incidence_prevalence_plot,
    Reff_plot,
    incidence_testing_plot,
    ensemble_incarr_Reff_plot

include("plotting-functions/simulation-ews_plots.jl")
export Reff_ews_plot,
    simulation_tau_distribution,
    Reff_plot,
    Reff_inc_plot,
    Reff_noise_plot,
    Reff_noise_inc_plot,
    Reff_noise_inc_test_plot

include("ensemble-sim_single-scenario_plots.jl")
export plot_all_single_scenarios

include("plotting-functions/simulation-optimal-ews_plots.jl")
export create_optimal_ews_plots

include("plotting-functions/hyperparam-debugging_plots.jl")
export hyperparam_debugging_Reff_plot

include("plotting-functions/accuracy-lines_plots.jl")
export prepare_line_plot_df!,
    line_plot

include("plotting-functions/ews-heatmap_plots.jl")

include("plotting-functions/survival_plots.jl")
export ews_survival_plot,
    ews_reff_histogram_plot

include("benchmark-functions.jl")

include("logging-utilities.jl")

@static if false
    include("../scripts/ensemble-sim.jl")
    include("../scripts/ensemble-sim_inferred-scenario-visualizations.jl")
    include("../scripts/tycho-visualization.jl")
    include("../scripts/ensemble-sim_ews-optimization.jl")
    include("../scripts/ensemble-sim_ews-multistart-optimization.jl")
    include("../scripts/benchmark_optimization_speed.jl")
    include("../manuscript/scripts/optimal-thresholds.jl")
end

# end
end
