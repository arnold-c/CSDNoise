#!/usr/bin/env julia

"""
Enhanced benchmark script to compare grid search vs multistart optimization
with comprehensive accuracy verification across all scenarios.

Usage:
    julia benchmark_optimization_speed.jl
"""

using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using DataFrames
using Dates
using Random
using Printf
using StyledStrings
using UnPack
using StructArrays
using SumTypes
using Try
using StaticArrays
using Distributions: Distributions
using ProgressMeter
using FLoops
using StatsBase
using CSV

function main()
    println(styled"{bold blue:Optimization Speed & Accuracy Benchmark}")
    println("="^60)

    # Set random seed for reproducibility
    Random.seed!(42)

    # Setup ensemble configuration (using realistic parameters from ensemble-sim_ews-optimization.jl)
    ensemble_model_type = ("seasonal-infectivity-import", "tau-leaping")

    burnin_years = 5
    nyears = 20
    burnin_time = 365.0 * burnin_years
    ensemble_time_specification = SimTimeParameters(;
        burnin = 365.0 * burnin_years, tmin = 0.0, tmax = 365.0 * nyears,
        tstep = 1.0,
    )

    ensemble_state_specification = StateParameters(
        500_000,
        Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
    )

    mu = calculate_mu(27)
    beta_mean = calculate_beta(
        R0, GAMMA, mu, 1, ensemble_state_specification.init_states.N
    )
    epsilon = calculate_import_rate(
        mu, R0, ensemble_state_specification.init_states.N
    )

    min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
        0.9,
        burnin_years * 2,
        ensemble_state_specification.init_states.S,
        ensemble_state_specification.init_states.N,
        R0,
        mu,
    )

    max_burnin_vaccination_coverage = 1.0

    ensemble_dynamics_specification = DynamicsParameterSpecification(
        beta_mean,
        0.0,
        cos,
        SIGMA,
        GAMMA,
        mu,
        27,
        epsilon,
        R0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        0.6,
        0.8,
    )

    null_dynamics_specification = DynamicsParameterSpecification(
        map(
            pn -> getproperty(ensemble_dynamics_specification, pn),
            filter(
                name ->
                name != :min_vaccination_coverage &&
                    name != :max_vaccination_coverage,
                propertynames(ensemble_dynamics_specification),
            ),
        )...,
        nothing,
        nothing,
    )

    # Reduced ensemble size for benchmarking speed
    ensemble_nsims = 30

    ensemble_specification = EnsembleSpecification(
        ensemble_model_type,
        ensemble_state_specification,
        ensemble_dynamics_specification,
        ensemble_time_specification,
        ensemble_nsims,
    )

    null_specification = EnsembleSpecification(
        ensemble_model_type,
        ensemble_state_specification,
        null_dynamics_specification,
        ensemble_time_specification,
        ensemble_nsims,
    )

    ensemble_outbreak_specification = OutbreakSpecification(
        5, 30, 500
    )

    # Define comprehensive test parameter space (reduced for speed)
    noise_specification_vec = [
        NoiseSpecification(PoissonNoise(1.0)),
        NoiseSpecification(PoissonNoise(7.0)),
    ]

    test_specification_vec = [
        IndividualTestSpecification(0.9, 0.9, 0),
        IndividualTestSpecification(0.95, 0.95, 0),
    ]

    percent_tested_vec = [1.0]

    ews_method_vec = [Backward]
    ews_aggregation_vec = [Day(28)]
    ews_bandwidth_vec = [Week(52)]
    ews_lag_days_vec = [1]

    ews_metric_specification_vec = create_combinations_vec(
        EWSMetricSpecification,
        (
            ews_method_vec,
            ews_aggregation_vec,
            ews_bandwidth_vec,
            ews_lag_days_vec,
        ),
    )

    ews_metric_vec = [
        "autocorrelation",
        "autocovariance",
        "coefficient_of_variation",
        "index_of_dispersion",
        "kurtosis",
        "mean",
        "skewness",
        "variance",
    ]

    ews_enddate_type_vec = [Reff_start]
    ews_threshold_window_vec = [ExpandingThresholdWindow]
    ews_threshold_percentile_vec = collect(0.5:0.02:0.99)  # Reduced resolution for speed
    ews_consecutive_thresholds_vec = collect(2:2:30)       # Reduced resolution for speed
    ews_threshold_burnin_vec = [Year(5)]

    specification_vecs = (;
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_threshold_percentile_vec,
        ews_consecutive_thresholds_vec,
        ews_metric_vec,
    )

    println(styled"Generating ensemble data...")
    data_arrs = generate_ensemble_data(ensemble_specification, null_specification, ensemble_outbreak_specification)

    # Calculate search space size
    n_scenarios = calculate_scenarios(specification_vecs)
    n_grid_points = calculate_grid_points(specification_vecs)
    total_grid_evaluations = n_scenarios * n_grid_points

    println(styled"Search space:")
    println(styled"  Scenarios: {cyan:$n_scenarios}")
    println(styled"  Grid points per scenario: {cyan:$n_grid_points}")
    println(styled"  Total grid evaluations: {yellow:$total_grid_evaluations}")

    # Compilation run - single scenario to pre-compile functions
    println(styled"\n{bold yellow:COMPILATION RUN}")
    println("-"^40)
    println(styled"Running single scenario to compile functions...")

    compilation_spec_vecs = (;
        noise_specification_vec = [specification_vecs.noise_specification_vec[1]],
        test_specification_vec = [specification_vecs.test_specification_vec[1]],
        percent_tested_vec = [specification_vecs.percent_tested_vec[1]],
        ews_metric_specification_vec = [specification_vecs.ews_metric_specification_vec[1]],
        ews_enddate_type_vec = [specification_vecs.ews_enddate_type_vec[1]],
        ews_threshold_window_vec = [specification_vecs.ews_threshold_window_vec[1]],
        ews_threshold_burnin_vec = [specification_vecs.ews_threshold_burnin_vec[1]],
        ews_threshold_percentile_vec = [specification_vecs.ews_threshold_percentile_vec[1]],
        ews_consecutive_thresholds_vec = [specification_vecs.ews_consecutive_thresholds_vec[1]],
        ews_metric_vec = [specification_vecs.ews_metric_vec[1]],
    )

    compilation_time = @elapsed _ = ews_hyperparam_optimization(
        compilation_spec_vecs,
        data_arrs;
        disable_time_check = true,
        force = true,
        return_df = true,
        save_results = false,
        verbose = false,
    )

    println(styled"Compilation completed in {yellow:$(round(compilation_time, digits=2))} seconds")

    # Benchmark grid search
    println(styled"\n{bold green:GRID SEARCH BENCHMARK}")
    println("-"^40)

    grid_time = @elapsed grid_results = ews_hyperparam_optimization(
        specification_vecs,
        data_arrs;
        disable_time_check = true,
        force = true,
        return_df = true,
        save_results = false,
        verbose = false,
    )

    println("Grid search complete")

    # Compilation run for multistart optimization
    println(styled"\n{bold yellow:MULTISTART COMPILATION RUN}")
    println("-"^40)
    println(styled"Running single multistart scenario to compile functions...")

    multistart_compilation_time = @elapsed _ = ews_multistart_optimization(
        compilation_spec_vecs,
        data_arrs;
        n_sobol_points = 1,
        maxeval = 10,
        executor = FLoops.ThreadedEx(),
        force = true,
        return_df = true,
        save_results = false,
        verbose = false,
        disable_time_check = true,
    )

    println(styled"Multistart compilation completed in {yellow:$(round(multistart_compilation_time, digits=2))} seconds")

    # Benchmark multistart with different configurations
    multistart_configs = [
        (n_sobol = 20, name = "Fast"),
        (n_sobol = 50, name = "Balanced"),
        (n_sobol = 100, name = "Thorough"),
        (n_sobol = 150, name = "Very Thorough"),
        # (n_sobol = 250, name = "Extremely Thorough"),
    ]

    println(styled"\n{bold green:MULTISTART OPTIMIZATION BENCHMARK}")
    println("-"^40)

    multistart_results = Dict()

    for config in multistart_configs
        n_sobol = config.n_sobol
        name = config.name

        println(styled"\n{bold:$name} ({cyan:$n_sobol} Sobol points):")

        ms_time = @elapsed results = ews_multistart_optimization(
            specification_vecs,
            data_arrs;
            percentile_bounds = (0.5, 0.99),  # Custom bounds
            consecutive_bounds = (2.0, 30.0),  # Custom bounds
            n_sobol_points = n_sobol,
            maxeval = 1000,
            executor = FLoops.ThreadedEx(),          # Executor for parallel processing
            force = true,
            return_df = true,
            save_results = false,
            verbose = false,
            disable_time_check = true
        )

        multistart_results[n_sobol] = (; time = ms_time, res = results)

        println("Multistart $name complete")
    end

    # Summary comparison - fix to specific scenario for fair comparison
    println(styled"\n{bold blue:SUMMARY - FIXED SCENARIO COMPARISON}")
    println("="^60)

    # Define reference scenario for comparison (first noise, first test spec, autocovariance metric)
    reference_scenario = (
        noise_specification = noise_specification_vec[1],  # NoiseSpecification(PoissonNoise(1.0))
        test_specification = test_specification_vec[1],    # IndividualTestSpecification(0.9, 0.9, 0)
        percent_tested = percent_tested_vec[1],            # 1.0
        ews_metric = "autocovariance",
        ews_metric_specification = ews_metric_specification_vec[1],
        ews_enddate_type = ews_enddate_type_vec[1],
        ews_threshold_window = ews_threshold_window_vec[1],
        ews_threshold_burnin = ews_threshold_burnin_vec[1],
    )

    println(styled"Reference scenario:")
    println(styled"  Noise: {cyan:$(get_noise_magnitude_description(reference_scenario.noise_specification))}")
    println(styled"  Test spec: {cyan:$(get_test_description(reference_scenario.test_specification))}")
    println(styled"  EWS metric: {cyan:$(reference_scenario.ews_metric)}")
    println()

    # Filter grid results to reference scenario
    grid_scenario_results = filter(
        row ->
        row.noise_specification == reference_scenario.noise_specification &&
            row.test_specification == reference_scenario.test_specification &&
            row.percent_tested == reference_scenario.percent_tested &&
            row.ews_metric == reference_scenario.ews_metric &&
            row.ews_metric_specification == reference_scenario.ews_metric_specification &&
            row.ews_enddate_type == reference_scenario.ews_enddate_type &&
            row.ews_threshold_window == reference_scenario.ews_threshold_window &&
            row.ews_threshold_burnin == reference_scenario.ews_threshold_burnin,
        grid_results
    )

    best_grid_acc = maximum(grid_scenario_results.accuracy)
    best_grid_idx = argmax(grid_scenario_results.accuracy)
    best_grid_row = grid_scenario_results[best_grid_idx, :]

    println(styled"Grid Search (reference scenario):")
    println(styled"  Time: {yellow:$(round(grid_time, digits=2))}s")
    println(styled"  Best accuracy: {green:$(round(best_grid_acc, digits=4))}")
    println(styled"  Best parameters (optimized):")
    println(styled"    Threshold percentile: {cyan:$(round(best_grid_row.ews_threshold_percentile, digits=3))}")
    println(styled"    Consecutive thresholds: {cyan:$(best_grid_row.ews_consecutive_thresholds)}")
    println(styled"    Sensitivity: {green:$(round(best_grid_row.sensitivity, digits=4))}")
    println(styled"    Specificity: {green:$(round(best_grid_row.specificity, digits=4))}")

    for config in multistart_configs
        n_sobol = config.n_sobol
        name = config.name
        results = multistart_results[n_sobol][:res]
        ms_time = multistart_results[n_sobol][:time]

        # Filter multistart results to same reference scenario
        ms_scenario_results = filter(
            row ->
            row.noise_specification == reference_scenario.noise_specification &&
                row.test_specification == reference_scenario.test_specification &&
                row.percent_tested == reference_scenario.percent_tested &&
                row.ews_metric == reference_scenario.ews_metric &&
                row.ews_metric_specification == reference_scenario.ews_metric_specification &&
                row.ews_enddate_type == reference_scenario.ews_enddate_type &&
                row.ews_threshold_window == reference_scenario.ews_threshold_window &&
                row.ews_threshold_burnin == reference_scenario.ews_threshold_burnin,
            results
        )

        if nrow(ms_scenario_results) == 0
            println(styled"\nMultistart ($name - $n_sobol Sobol points): {red:No results for reference scenario}")
            continue
        end

        best_ms_acc = maximum(ms_scenario_results.accuracy)
        acc_diff = best_ms_acc - best_grid_acc

        speedup = round(grid_time / ms_time, digits = 1)
        reduction_pct = round(100 * (1 - ms_time / grid_time), digits = 1)

        # Find best accuracy row in filtered results
        best_idx = argmax(ms_scenario_results.accuracy)
        best_row = ms_scenario_results[best_idx, :]

        println(styled"\nMultistart ($name - $n_sobol Sobol points, reference scenario):")
        println(styled"  Time: {yellow:$(round(ms_time, digits=2))}s")
        println(styled"  Speedup: {green:$(speedup)}x")
        println(styled"  Time reduction: {green:$(reduction_pct)}%")
        println(styled"  Best accuracy: {green:$(round(best_ms_acc, digits=4))}")
        println(styled"  Δ Accuracy vs Grid: {yellow:$(round(acc_diff, digits=4))}")

        # Print optimized parameters (only the ones that vary)
        println(styled"  Best parameters (optimized):")
        println(styled"    Threshold percentile: {cyan:$(round(best_row.ews_threshold_percentile, digits=3))}")
        println(styled"    Consecutive thresholds: {cyan:$(best_row.ews_consecutive_thresholds)}")
        println(styled"    Sensitivity: {green:$(round(best_row.sensitivity, digits=4))}")
        println(styled"    Specificity: {green:$(round(best_row.specificity, digits=4))}")
    end

    # NEW: Comprehensive accuracy comparison across all scenarios
    println(styled"\n{bold blue:COMPREHENSIVE ACCURACY VERIFICATION}")
    println("="^60)

    # Perform detailed comparison
    comparison_results = compare_all_scenarios(
        grid_results,
        multistart_results,
        specification_vecs
    )

    # Display summary statistics
    display_accuracy_comparison_summary(comparison_results)

    # Generate detailed report
    generate_accuracy_verification_report(
        comparison_results,
        grid_time,
        multistart_results
    )

    return grid_results, multistart_results, comparison_results
end

# Run benchmark if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
