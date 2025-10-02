export CachedSimulationData,
    OptimizationScenario,
    GridSearchScenario,
    EWSClassificationResults,
    OptimizationTracker,
    OptimizedValues,
    OptimizationResult

# TODO: See if this should be removed or updated
"""
    CachedSimulationData

Pre-computed simulation data that can be reused across parameter evaluations.
This avoids expensive recomputation of noise arrays and test arrays.
"""
Base.@kwdef struct CachedSimulationData
    testarr::Array{Int64, 3}
    null_testarr::Array{Int64, 3}
    thresholds::Vector{Matrix{Int64}}
    ews_metrics::Vector{EWSMetrics}
    null_ews_metrics::Vector{EWSMetrics}
end

"""
    OptimizationScenario

Struct representing a single optimization scenario with all necessary parameters
for EWS hyperparameter optimization.
"""
Base.@kwdef struct OptimizationScenario
    ensemble_specification::EnsembleSpecification
    null_specification::EnsembleSpecification
    noise_specification::NoiseSpecification
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_threshold_burnin::Dates.Day
    ews_metric::String
end

function OptimizationScenario(
        ensemble_specification::EnsembleSpecification,
        null_specification::EnsembleSpecification,
        noise_specification::NoiseSpecification,
        test_specification::IndividualTestSpecification,
        percent_tested::Float64,
        ews_metric_specification::EWSMetricSpecification,
        ews_enddate_type::EWSEndDateType,
        ews_threshold_window::EWSThresholdWindowType,
        threshold_burnin::P,
        ews_metric::String,
    ) where {P <: Dates.Period}

    ews_threshold_burnin = Dates.Day(round(Dates.days(threshold_burnin)))

    return OptimizationScenario(
        ensemble_specification,
        null_specification,
        noise_specification,
        test_specification,
        percent_tested,
        ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window,
        ews_threshold_burnin,
        ews_metric,
    )
end

"""
    GridSearchScenario

Scenario for grid search including both base scenario and grid parameters.
"""
Base.@kwdef struct GridSearchScenario
    ensemble_specification::EnsembleSpecification
    null_specification::EnsembleSpecification
    noise_specification::NoiseSpecification
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_threshold_burnin::Dates.Day
    ews_metric::String
    threshold_quantile::Float64
    consecutive_thresholds::Int64
end

function GridSearchScenario(
        ensemble_specification::EnsembleSpecification,
        null_specification::EnsembleSpecification,
        noise_specification::NoiseSpecification,
        test_specification::IndividualTestSpecification,
        percent_tested::Float64,
        ews_metric_specification::EWSMetricSpecification,
        ews_enddate_type::EWSEndDateType,
        ews_threshold_window::EWSThresholdWindowType,
        threshold_burnin::P,
        ews_metric::String,
        threshold_quantile::Float64,
        consecutive_thresholds::Int64
    ) where {P <: Dates.Period}

    ews_threshold_burnin = Dates.Day(round(Dates.days(threshold_burnin)))

    return GridSearchScenario(
        ensemble_specification,
        null_specification,
        noise_specification,
        test_specification,
        percent_tested,
        ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window,
        ews_threshold_burnin,
        ews_metric,
        threshold_quantile,
        consecutive_thresholds
    )
end


"""
    EWSClassificationResults

Stores binary classification results from Early Warning Signal (EWS) detection analysis.

Contains the confusion matrix components and total counts for evaluating EWS performance
across emergent (positive class) and null (negative class) simulations.

# Fields
- `true_positives::Float64`: Number of emergent simulations correctly identified by EWS metric
- `true_negatives::Float64`: Number of null simulations correctly identified as by EWS metric
- `false_positives::Float64`: Number of null simulations incorrectly identified by EWS metric
- `false_negatives::Float64`: Number of emergent simulations incorrectly identified by EWS metric
- `n_emergent_sims::Int64`: Total number of emergent (positive class) simulations
- `n_null_sims::Int64`: Total number of null (negative class) simulations

# Notes
The classification counts are stored as Float64 to support weighted or fractional classifications,
while the total simulation counts remain as integers. This struct serves as an intermediate
representation for calculating performance metrics like sensitivity, specificity, and accuracy.
"""
Base.@kwdef struct EWSClassificationResults
    true_positives::Float64
    true_negatives::Float64
    false_positives::Float64
    false_negatives::Float64
    n_emergent_sims::Int64
    n_null_sims::Int64
end

"""
    OptimizationTracker

Mutable struct to track the best solution and its metrics during optimization.
"""
Base.@kwdef mutable struct OptimizationTracker
    best_loss::Float64 = Inf
    best_accuracy::Float64 = 0.0
    best_sensitivity::Float64 = 0.0
    best_specificity::Float64 = 0.0
    best_params::Vector{Float64} = Float64[]
end

Base.@kwdef struct OptimizedValues
    threshold_quantile::Float64
    consecutive_thresholds::Int64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
end

Base.@kwdef struct OptimizationResult
    ensemble_specification::EnsembleSpecification
    null_specification::EnsembleSpecification
    noise_specification::NoiseSpecification
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_threshold_burnin::Dates.Day
    ews_metric::String
    threshold_quantile::Float64
    consecutive_thresholds::Int64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
end
