using Test
using CSDNoise
using StructArrays
using Dates

@testset "find_missing_scenarios" begin
    time_parameters = SimTimeParameters(
        burnin = 365.0 * 5,
        tmin = 0.0,
        tmax = 365.0 * 20.0,
        tstep = 1.0
    )

    births_per_k_pop = 27.0
    burnin_target_Reff = 0.9

    common_disease_dynamics_parameters = CommonDiseaseDynamicsParameters(
        births_per_k_pop = births_per_k_pop,
        nsims = 100,
        burnin_target_Reff = burnin_target_Reff
    )

    initial_Reff = 0.8
    test_R0 = 5.0
    test_initial_s_prop = round(initial_Reff / test_R0; digits = 2)

    state_parameters_1 = StateParameters(
        100_000,
        Dict(
            :s_prop => test_initial_s_prop,
            :e_prop => 0.0,
            :i_prop => 0.0,
            :r_prop => 1.0 - test_initial_s_prop
        )
    )

    state_parameters_2 = StateParameters(
        200_000,
        Dict(
            :s_prop => test_initial_s_prop,
            :e_prop => 0.0,
            :i_prop => 0.0,
            :r_prop => 1.0 - test_initial_s_prop
        )
    )

    dynamics_parameters = TargetDiseaseDynamicsParameters(
        R_0 = test_R0,
        latent_period_days = 5.0,
        infectious_duration_days = 10.0,
        beta_force = 0.0,
        seasonality = SeasonalityFunction(CosineSeasonality()),
        min_post_burnin_vaccination_coverage = 0.4,
        max_post_burnin_vaccination_coverage = 0.6,
        max_burnin_vaccination_coverage = 1.0
    )

    dynamical_noise_spec = DynamicalNoiseSpecification(
        R_0 = 1.28,
        latent_period = 1.0,
        duration_infection = 4.8,
        correlation = "in-phase",
        poisson_component = 0.15
    )

    ensemble_spec_1 = create_ensemble_specifications(
        time_parameters,
        state_parameters_1,
        dynamics_parameters,
        common_disease_dynamics_parameters,
        dynamical_noise_spec
    )

    ensemble_spec_2 = create_ensemble_specifications(
        time_parameters,
        state_parameters_2,
        dynamics_parameters,
        common_disease_dynamics_parameters,
        dynamical_noise_spec
    )

    test_spec = IndividualTestSpecification(
        sensitivity = 0.9,
        specificity = 0.95,
        test_result_lag = 1
    )

    ews_spec = EWSMetricSpecification(
        EWSMethod(Backward()),
        Day(7),
        Day(28),
        0
    )

    @testset "GridSearchScenario - empty completed results" begin
        scenarios = StructArray(
            [
                GridSearchScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1
                ),
                GridSearchScenario(
                    ensemble_specification = ensemble_spec_2,
                    noise_level = 0.2,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1
                ),
            ]
        )

        completed = StructArray(OptimizationResult[])

        missing = find_missing_scenarios(scenarios, completed)

        @test length(missing) == 2
        @test missing == scenarios
    end

    @testset "GridSearchScenario - all scenarios completed" begin
        scenarios = StructArray(
            [
                GridSearchScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1
                ),
            ]
        )

        completed = StructArray(
            [
                OptimizationResult(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1,
                    accuracy = 0.85,
                    sensitivity = 0.8,
                    specificity = 0.9
                ),
            ]
        )

        missing = find_missing_scenarios(scenarios, completed)

        @test length(missing) == 0
    end

    @testset "GridSearchScenario - partial completion" begin
        scenarios = StructArray(
            [
                GridSearchScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1
                ),
                GridSearchScenario(
                    ensemble_specification = ensemble_spec_2,
                    noise_level = 0.2,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1
                ),
                GridSearchScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.3,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "variance",
                    threshold_quantile = 0.9,
                    consecutive_thresholds = 2
                ),
            ]
        )

        completed = StructArray(
            [
                OptimizationResult(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1,
                    accuracy = 0.85,
                    sensitivity = 0.8,
                    specificity = 0.9
                ),
            ]
        )

        missing = find_missing_scenarios(scenarios, completed)

        @test length(missing) == 2
        @test missing[1].ensemble_specification == ensemble_spec_2
        @test missing[1].noise_level == 0.2
        @test missing[2].ensemble_specification == ensemble_spec_1
        @test missing[2].noise_level == 0.3
    end

    @testset "OptimizationScenario - empty completed results" begin
        scenarios = StructArray(
            [
                OptimizationScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean"
                ),
            ]
        )

        completed = StructArray(OptimizationResult[])

        missing = find_missing_scenarios(scenarios, completed)

        @test length(missing) == 1
        @test missing == scenarios
    end

    @testset "OptimizationScenario - all scenarios completed" begin
        scenarios = StructArray(
            [
                OptimizationScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean"
                ),
            ]
        )

        completed = StructArray(
            [
                OptimizationResult(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1,
                    accuracy = 0.85,
                    sensitivity = 0.8,
                    specificity = 0.9
                ),
            ]
        )

        missing = find_missing_scenarios(scenarios, completed)

        @test length(missing) == 0
    end

    @testset "OptimizationScenario - partial completion" begin
        scenarios = StructArray(
            [
                OptimizationScenario(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean"
                ),
                OptimizationScenario(
                    ensemble_specification = ensemble_spec_2,
                    noise_level = 0.2,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "variance"
                ),
            ]
        )

        completed = StructArray(
            [
                OptimizationResult(
                    ensemble_specification = ensemble_spec_1,
                    noise_level = 0.1,
                    noise_type_description = :static,
                    test_specification = test_spec,
                    percent_tested = 0.5,
                    ews_metric_specification = ews_spec,
                    ews_enddate_type = EWSEndDateType(ReffStart()),
                    ews_threshold_window = EWSThresholdWindowType(ExpandingThresholdWindow()),
                    ews_metric = "mean",
                    threshold_quantile = 0.95,
                    consecutive_thresholds = 1,
                    accuracy = 0.85,
                    sensitivity = 0.8,
                    specificity = 0.9
                ),
            ]
        )

        missing = find_missing_scenarios(scenarios, completed)

        @test length(missing) == 1
        @test missing[1].ensemble_specification == ensemble_spec_2
        @test missing[1].ews_metric == "variance"
    end
end
