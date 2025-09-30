using Dates: Dates

export CachedSimulationData,
    OptimizationScenario,
    GridSearchScenario,
    OptimizationResult

# TODO: See if this should be removed or updated
"""
    CachedSimulationData

Pre-computed simulation data that can be reused across parameter evaluations.
This avoids expensive recomputation of noise arrays and test arrays.
"""
struct CachedSimulationData
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
struct OptimizationScenario
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
struct GridSearchScenario
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

struct OptimizationResult
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
