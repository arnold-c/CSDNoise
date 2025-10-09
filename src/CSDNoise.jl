"""
Note that each file handles exporting its local function for the API.
"""
module CSDNoise
using DispatchDoctor: @unstable
using Base: Base.@kwdef
using DrWatson: DrWatson
using Random: Random
using Distributions: Distributions
using LinearAlgebra: LinearAlgebra
using StatsBase: StatsBase
using StaticArrays: StaticArrays
using StructArrays: StructVector
using LabelledArrays: LabelledArrays
import DataFrames as DF
using Dates: Dates
using LightSumTypes: LightSumTypes
using UnPack: @unpack
using Bumper: @no_escape, @alloc
using Printf: @sprintf
using Try: Try
using TryExperimental: TryExperimental
using MultistartOptimization: MultistartOptimization
using NLopt: NLopt
using FLoops: FLoops
using ProgressMeter: ProgressMeter
using BangBang: BangBang
using JLD2: JLD2
using StyledStrings: StyledStrings
import REPL.TerminalMenus as RTM
using Logging: Logging
using LoggingExtras: LoggingExtras
using GLMakie

# Constants
include("./constants/dynamics-constants.jl")

# Types
include("./types/test-specifications.jl")
include("./constants/test-constants.jl") # depends on test specifications, and some types depend on them
include("./types/time-parameters.jl")
include("./types/dynamics-parameters.jl")
include("./types/state-parameters.jl")
include("./types/noise-specifications.jl")
include("./types/ensemble-specifications.jl")
include("./types/outbreak-specifications.jl")
include("./types/ews-specifications.jl")
include("./types/scenario-specifications.jl")
include("./types/simulation-results.jl")
include("./types/optimization-types.jl")

# Simulation core
include("./simulation/transmission-functions.jl")
include("./constants/dynamics-constants_calculated.jl") # depends on transmission-functions
include("./simulation/endemic-equilibrium.jl")
include("./simulation/seir-model.jl")
include("./simulation/seir-mean-incidence.jl")
include("./simulation/trim-seir-results.jl")

# Shared utilities (used across multiple modules)
include("./utilities/benchmark-functions.jl")
include("./utilities/create-combinations.jl")
include("./utilities/helpers.jl")
include("./utilities/create_emergent_null_specifications.jl")
include("./utilities/logging-utilities.jl")
include("./utilities/test-description.jl")
include("./utilities/test_multistart_vs_gridsearch.jl")

# Noise
include("./noise/noise-description.jl")
include("./noise/noise-generation.jl")
include("./noise/noise-mean-incidence.jl")
include("./noise/noise-parameters-optimization.jl")
include("./noise/noise-recreation.jl")
include("./noise/noise_dynamics_parameters-creation.jl")

# Vaccination
include("./vaccination/vaccination-distribution-sample.jl")
include("./vaccination/vaccination-emergent-level.jl")
include("./vaccination/vaccination-range.jl")

# Testing
include("./diagnostic-testing/calculate-num-positive.jl")
include("./diagnostic-testing/calculate-num-tested.jl")
include("./diagnostic-testing/create-test-positive-vectors.jl")

# Detection
include("./detection/outbreak-thresholds.jl")
include("./detection/Reff_thresholds.jl")
include("./detection/shared_threshold-functions.jl")

# EWS Metrics
include("./ews-metrics/ews-alerts.jl")
include("./ews-metrics/ews-bandwidths.jl")
include("./ews-metrics/ews-enddates.jl")
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
## Grid search
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
include("./survival/create-survival-data.jl")
include("./survival/simulate-survival-data.jl")
include("./survival/survival-detection-indices.jl")

# Plotting Functions
include("./plotting-functions/plotting-setup_makie.jl")
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
include("./plotting-functions/survival/ews-survival_facet-parameter-preparation.jl")
include("./plotting-functions/survival/ews-survival_facet.jl")
include("./plotting-functions/survival/ews-survival_plot-wrapper.jl")
include("./plotting-functions/survival/ews-survival_Reff-histogram.jl")
include("./plotting-functions/survival/ews-survival_simulate-and-plot.jl")
include("./plotting-functions/survival/ews-survival_survival-lines.jl")
include("./plotting-functions/survival/Reff_histogram_plot.jl")

@static if false
    include("../scripts/gridsearch-optimization.jl")
    include("../scripts/ensemble-sim_inferred-scenario-visualizations.jl")
    include("../scripts/ensemble-sim_ews-optimization.jl")
    include("../scripts/ensemble-sim_ews-multistart-optimization.jl")
    include("../scripts/benchmark_optimization_speed.jl")
    include("../manuscript/scripts/optimal-thresholds.jl")
end

end
