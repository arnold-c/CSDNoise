"""
Note that each file handles exporting its local function for the API.
"""
module CSDNoise
using DispatchDoctor
# @stable default_mode = "warn" default_union_limit = 2 begin

# Utilities
include("./utilities/helpers.jl")
include("./utilities/test-description.jl")
include("./utilities/benchmark-functions.jl")
include("./utilities/logging-utilities.jl")

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


# Testing
include("./testing/calculate-moving-average.jl")
include("./testing/calculate-positive.jl")
include("./testing/calculate-test-positivity.jl")
include("./testing/calculate-tested.jl")
include("./testing/infer-test-positive.jl")
include("./testing/testing-arrays.jl")
include("./testing/testing-vectors.jl")

include("./simulation/ensemble-functions.jl")
export create_combinations_vec,
    create_ensemble_spec_combinations,
    run_ensemble_jump_prob,
    run_jump_prob,
    get_ensemble_file

# Shared utilities (used across multiple modules)
include("filter-seir-results.jl")

# Noise
include("./noise/noise-description.jl")
include("./noise/noise-generation.jl")
include("./noise/noise-mean-incidence.jl")
include("./noise/noise-parameters-optimization.jl")
include("./noise/noise-recreation.jl")

# Vaccination
include("./vaccination/vaccination-distribution-sample.jl")
include("./vaccination/vaccination-emergent-level.jl")
include("./vaccination/vaccination-range.jl")

# Plotting Functions
include("./plotting-functions/helpers_plots.jl")
include("./plotting-functions/tau-heatmap_plots.jl")
include("./plotting-functions/lead-time_plots.jl")
include("./plotting-functions/single-simulation_plots.jl")
include("./plotting-functions/simulation-ews_plots.jl")
include("./plotting-functions/single-scenario_plots.jl")
include("./plotting-functions/simulation-optimal-ews_plots.jl")
include("./plotting-functions/hyperparam-debugging_plots.jl")
include("./plotting-functions/accuracy-lines_data-preparation.jl")
include("./plotting-functions/accuracy-lines_plots.jl")
include("./plotting-functions/ews-heatmap_plots.jl")

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
