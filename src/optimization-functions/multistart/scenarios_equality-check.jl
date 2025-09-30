export scenarios_equal

"""
    scenarios_equal(scenario1, scenario2)

Compare two scenarios for equality, handling floating point comparisons properly.
"""
function scenarios_equal(scenario1::OptimizationScenario, scenario2::OptimizationScenario)
    scenario1 == scenario2 && return true

    # If the previous equality fails, check all struct properties
    # Doing in a loop results in runtime dispatch (type groundedness as type in loop not determined just
    # by the method argument types)
    scenario1.ensemble_specification != scenario2.ensemble_specification && return false
    scenario1.null_specification != scenario2.null_specification && return false
    scenario1.noise_specification != scenario2.noise_specification && return false
    scenario1.test_specification != scenario2.test_specification && return false
    # Check the error isn't due to minor floating point errors
    !isapprox(scenario1.percent_tested, scenario2.percent_tested) && return false
    scenario1.ews_metric_specification != scenario2.ews_metric_specification && return false
    scenario1.ews_enddate_type != scenario2.ews_enddate_type && return false
    scenario1.ews_threshold_window != scenario2.ews_threshold_window && return false
    scenario1.ews_threshold_burnin != scenario2.ews_threshold_burnin && return false
    scenario1.ews_metric != scenario2.ews_metric && return false
    return true
end
