"""
Note that each file handles exporting its local function for the API.
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
include("./types/optimization-types.jl")

# Type descriptions
include("./vaccination-rate-calculations.jl")
include("./test-description.jl")
include("./noise-description.jl")

# Simulation core
include("./simulation/transmission-functions.jl")

# Constants that rely on transmission function definitions
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

# Optimization functions
## Generic
include("./optimization-functions/results_retrieval.jl")
include("./optimization-functions/results_classification.jl")
## Multistart
include("./optimization-functions/multistart/checkpoint_cleanup.jl")
include("./optimization-functions/multistart/checkpoint_save.jl")
include("./optimization-functions/multistart/checkpoint_load.jl")
include("./optimization-functions/multistart/helpers_ews-calculation.jl")
include("./optimization-functions/multistart/objective-function.jl")
include("./optimization-functions/multistart/optimization_cached-data.jl")
include("./optimization-functions/multistart/optimization_batches.jl")
include("./optimization-functions/multistart/optimization_wrapper.jl")
include("./optimization-functions/multistart/optimization_single-scenario.jl")
include("./optimization-functions/multistart/results_create-empty-df.jl")
include("./optimization-functions/multistart/results_summary.jl")
include("./optimization-functions/multistart/results_validation.jl")
include("./optimization-functions/multistart/results_retrieval.jl")
include("./optimization-functions/multistart/results_merge-dfs.jl")
include("./optimization-functions/multistart/results_df-conversion.jl")
include("./optimization-functions/multistart/results_save.jl")
include("./optimization-functions/multistart/scenarios_creation.jl")
include("./optimization-functions/multistart/scenarios_find-missing.jl")
include("./optimization-functions/multistart/scenarios_confirmation.jl")
include("./optimization-functions/multistart/scenarios_df-row-conversion.jl")
include("./optimization-functions/multistart/scenarios_equality-check.jl")
## Structvector
### Grid search
include("./optimization-functions/gridsearch/checkpoint_save.jl")
include("./optimization-functions/gridsearch/gridsearch_wrapper.jl")
include("./optimization-functions/gridsearch/helpers_ews-calculation.jl")
include("./optimization-functions/gridsearch/results_retrieval.jl")
include("./optimization-functions/gridsearch/results_save.jl")
include("./optimization-functions/gridsearch/scenario_confirmation.jl")
include("./optimization-functions/gridsearch/scenario_creation.jl")
include("./optimization-functions/gridsearch/scenario_evaluation.jl")
include("./optimization-functions/gridsearch/scenario_find-missing.jl")


# Survival
## Data
include("./survival/create-survival-data.jl")
include("./survival/simulate-survival-data.jl")
include("./survival/survival-detection-indices.jl")
## Plots
include("./plotting-functions/survival/ews-survival_facet-parameter-preparation.jl")
include("./plotting-functions/survival/ews-survival_facet.jl")
include("./plotting-functions/survival/ews-survival_plot-wrapper.jl")
include("./plotting-functions/survival/ews-survival_Reff-histogram.jl")
include("./plotting-functions/survival/ews-survival_simulate-and-plot.jl")
include("./plotting-functions/survival/ews-survival_survival-lines.jl")
include("./plotting-functions/survival/Reff_histogram_plot.jl")


include("./constants/test-constants.jl")

include("./simulation/seir-model.jl")

# Detection
include("./detection/outbreak-thresholds.jl")
include("./detection/Reff_thresholds.jl")
include("./detection/shared_threshold-functions.jl")


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

include("./simulation/ensemble-functions.jl")
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
