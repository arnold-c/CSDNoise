export CachedSimulationData,
    OptimizationScenario,
    GridSearchSpecificationVecs,
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
    noise_level::Float64
    noise_type_description::Symbol
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_metric::String
    function OptimizationScenario(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            ews_metric_specification,
            ews_enddate_type,
            ews_threshold_window,
            ews_metric,
        )
        @assert noise_type_description in [:static, :dynamic] "The noise type must be either :static or :dynamic. Received $noise_type_description"
        return new(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            ews_metric_specification,
            ews_enddate_type,
            ews_threshold_window,
            ews_metric,
        )
    end
end

"""
    GridSearchSpecificationVecs

Vectors of parameter values for comprehensive grid search optimization of EWS hyperparameters.

This struct contains vectors of all possible parameter values that will be used to generate
the complete grid search space for Early Warning Signal (EWS) optimization. Each vector
represents the range of values to explore for a specific parameter, and the Cartesian
product of all vectors defines the full search space.

The grid search will evaluate EWS performance across all combinations of these parameter
values to identify optimal hyperparameter settings for outbreak detection.

# Fields
- `ensemble_specification_vec::Vector{EnsembleSpecification}`: Ensemble simulation configurations to test
- `noise_level_vec::Vector{Float64}`: Noise levels to apply to simulation data
- `test_specification_vec::Vector{IndividualTestSpecification}`: Diagnostic testing configurations
- `percent_tested_vec::Vector{Float64}`: Percentages of population tested (0.0 to 1.0)
- `ews_metric_specification_vec::Vector{EWSMetricSpecification}`: EWS metric calculation specifications
- `ews_enddate_type_vec::Vector{EWSEndDateType}`: Methods for determining EWS calculation end dates
- `ews_threshold_window_vec::Vector{EWSThresholdWindowType}`: Window types for threshold calculations
- `ews_metric_vec::Vector{String}`: Names of EWS metrics to evaluate
- `ews_threshold_quantile_vec::Vector{Float64}`: Quantile values for threshold determination
- `ews_consecutive_thresholds_vec::Vector{Int64}`: Numbers of consecutive threshold exceedances required

# Example
```julia
# Create grid search specification vectors
grid_specs = GridSearchSpecificationVecs(
    ensemble_specification_vec = [ensemble_spec1, ensemble_spec2],
    noise_level_vec = [1.0, 2.0, 3.0],
    test_specification_vec = [test_spec],
    percent_tested_vec = [0.05, 0.10, 0.15],
    ews_metric_specification_vec = [ews_spec],
    ews_enddate_type_vec = [EWSEndDateType("fixed")],
    ews_threshold_window_vec = [EWSThresholdWindowType("rolling")],
    ews_metric_vec = ["variance", "autocorrelation"],
    ews_threshold_quantile_vec = [0.90, 0.95, 0.99],
    ews_consecutive_thresholds_vec = [1, 2, 3]
)

# Total combinations: 2 × 3 × 1 × 3 × 1 × 1 × 1 × 2 × 3 × 3 = 324 scenarios
```

# See Also
- [`GridSearchScenario`](@ref): Individual scenario generated from these specification vectors
- [`OptimizationScenario`](@ref): Base optimization scenario structure
"""
Base.@kwdef struct GridSearchSpecificationVecs
    ensemble_specification_vec::Vector{EnsembleSpecification}
    noise_level_vec::Vector{Float64}
    noise_type_description_vec::Vector{Symbol}
    test_specification_vec::Vector{IndividualTestSpecification}
    percent_tested_vec::Vector{Float64}
    ews_metric_specification_vec::Vector{EWSMetricSpecification}
    ews_enddate_type_vec::Vector{EWSEndDateType}
    ews_threshold_window_vec::Vector{EWSThresholdWindowType}
    ews_metric_vec::Vector{String}
    ews_threshold_quantile_vec::Vector{Float64}
    ews_consecutive_thresholds_vec::Vector{Int64}
end

"""
    GridSearchScenario

Scenario for grid search including both base scenario and grid parameters.
"""
Base.@kwdef struct GridSearchScenario
    ensemble_specification::EnsembleSpecification
    noise_level::Float64
    noise_type_description::Symbol
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_metric::String
    threshold_quantile::Float64
    consecutive_thresholds::Int64
    function GridSearchScenario(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            ews_metric_specification,
            ews_enddate_type,
            ews_threshold_window,
            ews_metric,
            threshold_quantile,
            consecutive_thresholds,
        )
        @assert noise_type_description in [:static, :dynamic] "The noise type must be either :static or :dynamic. Received $noise_type_description"
        return new(
            ensemble_specification,
            noise_level,
            noise_type_description,
            test_specification,
            percent_tested,
            ews_metric_specification,
            ews_enddate_type,
            ews_threshold_window,
            ews_metric,
            threshold_quantile,
            consecutive_thresholds,

        )
    end
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

"""
    OptimizedValues

Optimal hyperparameter values and performance metrics from EWS optimization.

This struct stores the results of hyperparameter optimization for Early Warning Signal (EWS)
detection, containing both the optimal parameter values found during optimization and the
corresponding performance metrics achieved with those parameters.

The optimization process searches for the combination of threshold quantile and consecutive
threshold requirements that maximizes classification performance for outbreak detection.

# Fields
- `threshold_quantile::Float64`: Optimal quantile value for EWS threshold determination (0.0 to 1.0)
- `consecutive_thresholds::Int64`: Optimal number of consecutive threshold exceedances required for alert
- `accuracy::Float64`: Overall classification accuracy achieved with optimal parameters (0.0 to 1.0)
- `sensitivity::Float64`: True positive rate (outbreak detection rate) with optimal parameters (0.0 to 1.0)
- `specificity::Float64`: True negative rate (false alarm avoidance) with optimal parameters (0.0 to 1.0)

# Example
```julia
# Create optimized values from optimization results
optimized = OptimizedValues(
    threshold_quantile = 0.95,
    consecutive_thresholds = 2,
    accuracy = 0.87,
    sensitivity = 0.82,
    specificity = 0.91
)

println("Optimal threshold: \$(optimized.threshold_quantile) quantile")
println("Consecutive alerts needed: \$(optimized.consecutive_thresholds)")
println("Performance - Accuracy: \$(optimized.accuracy), Sensitivity: \$(optimized.sensitivity)")
```

# See Also
- [`OptimizationResult`](@ref): Complete optimization results including scenario details
- [`OptimizationTracker`](@ref): Mutable tracker used during optimization process
- [`EWSClassificationResults`](@ref): Raw classification results used to calculate performance metrics
"""
Base.@kwdef struct OptimizedValues
    threshold_quantile::Float64
    consecutive_thresholds::Int64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
end

"""
    OptimizationResult

Complete results from EWS hyperparameter optimization including scenario details and performance.

This struct provides a comprehensive record of an optimization run, storing both the complete
scenario configuration used during optimization and the optimal hyperparameter values and
performance metrics discovered. It serves as the primary output format for optimization
functions and enables reproducible analysis of optimization results.

The result includes all parameters needed to recreate the optimization scenario, making it
suitable for result storage, comparison across different scenarios, and further analysis.
This is also used when new scenarios are prepared to be run, checking if the results already
exist, and if they do, load the existing results to avoid redundant computations.

# Fields

## Scenario Configuration
- `ensemble_specification::EnsembleSpecification`: Emergent outbreak ensemble configuration
- `null_specification::EnsembleSpecification`: Null (no outbreak) ensemble configuration
- `noise_specification::NoiseSpecification`: Noise model applied to simulation data
- `test_specification::IndividualTestSpecification`: Diagnostic testing configuration
- `percent_tested::Float64`: Percentage of population tested (0.0 to 1.0)

## EWS Configuration
- `ews_metric_specification::EWSMetricSpecification`: EWS metric calculation specification
- `ews_enddate_type::EWSEndDateType`: Method for determining EWS calculation end dates
- `ews_threshold_window::EWSThresholdWindowType`: Window type for threshold calculations
- `ews_threshold_burnin::Dates.Day`: Burn-in period before threshold calculations begin
- `ews_metric::String`: Name of the EWS metric used

## Optimization Results
- `threshold_quantile::Float64`: Optimal quantile value for threshold determination
- `consecutive_thresholds::Int64`: Optimal number of consecutive threshold exceedances required
- `accuracy::Float64`: Overall classification accuracy with optimal parameters
- `sensitivity::Float64`: True positive rate (outbreak detection rate) with optimal parameters
- `specificity::Float64`: True negative rate (false alarm avoidance) with optimal parameters

# Example
```julia
# Create optimization result
result = OptimizationResult(
    ensemble_specification = emergent_ensemble,
    null_specification = null_ensemble,
    noise_specification = noise_spec,
    test_specification = test_spec,
    percent_tested = 0.10,
    ews_metric_specification = ews_spec,
    ews_enddate_type = EWSEndDateType("fixed"),
    ews_threshold_window = EWSThresholdWindowType("rolling"),
    ews_threshold_burnin = Dates.Day(30),
    ews_metric = "variance",
    threshold_quantile = 0.95,
    consecutive_thresholds = 2,
    accuracy = 0.87,
    sensitivity = 0.82,
    specificity = 0.91
)

println("Optimization achieved \$(result.accuracy) accuracy")
println("Optimal parameters: \$(result.threshold_quantile) quantile, \$(result.consecutive_thresholds) consecutive alerts")
```

# See Also
- [`OptimizedValues`](@ref): Subset containing only the optimal parameters and performance metrics
- [`OptimizationScenario`](@ref): Input scenario specification for optimization
- [`GridSearchScenario`](@ref): Extended scenario specification for grid search optimization
"""
Base.@kwdef struct OptimizationResult
    ensemble_specification::EnsembleSpecification
    noise_level::Float64
    noise_type_description::Symbol
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_metric::String
    threshold_quantile::Float64
    consecutive_thresholds::Int64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
end
