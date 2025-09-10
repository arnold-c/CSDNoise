using CSDNoise
using InteractiveUtils

# Create a vector of OptimizationScenario instances
function create_optimization_scenarios()
    scenarios = OptimizationScenario[]

    # Add different scenarios with different concrete types
    push!(
        scenarios, OptimizationScenario(
            noise_specification = NoiseSpecification(PoissonNoise(1.0)),
            test_specification = IndividualTestSpecification(0.9, 0.9, 0),
            percent_tested = 0.1,
            ews_metric_specification = EWSMetricSpecification(Backward, Dates.Day(7), Dates.Day(28), 1),
            ews_enddate_type = Reff_start,
            ews_threshold_window = ExpandingThresholdWindow,
            ews_threshold_burnin = Dates.Day(30),
            ews_metric = "variance"
        )
    )

    push!(
        scenarios, OptimizationScenario(
            noise_specification = NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.1, 0.1, 0.3)),
            test_specification = IndividualTestSpecification(0.95, 0.95, 1),
            percent_tested = 0.2,
            ews_metric_specification = EWSMetricSpecification(Centered, Dates.Day(14), Dates.Day(35), 2),
            ews_enddate_type = Outbreak_start,
            ews_threshold_window = RollingThresholdWindow,
            ews_threshold_burnin = Dates.Year(1),
            ews_metric = "autocovariance"
        )
    )

    return scenarios
end

# Test type stability
@code_warntype create_optimization_scenarios()
