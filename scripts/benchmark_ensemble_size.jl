#!/usr/bin/env julia

"""
Benchmark multistart optimization with 100 Sobol points comparing:
- 100 simulations vs 1000 simulations
- Single scenario focus for clear comparison

Usage:
    julia --startup-file=no --history-file=no enchmark_ensemble_size.jl
"""

using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using DataFrames
using Dates
using Random
using Printf
using StyledStrings
using Statistics
using CSV
using FLoops
using UnPack

function main()
    println(styled"{bold blue:Ensemble Size Benchmark (100 vs 1000 vs 10000 simulations)}")
    println("="^70)

    # Set seed for reproducibility
    Random.seed!(42)

    # Define single test scenario (autocovariance with moderate noise)
    test_scenario = create_single_test_scenario()

    # Define ensemble sizes to compare
    ensemble_sizes = [100, 1000, 10000]
    n_sobol_points = 100

    println(styled"Configuration:")
    println(styled"  Scenario: {cyan:Autocovariance, Poisson(1.0x), Test(0.9,0.9)}")
    println(styled"  Sobol points: {cyan:$n_sobol_points}")
    println(styled"  Ensemble sizes: {cyan:$(join(ensemble_sizes, \", \"))}")
    println()

    results = Dict()

    # Benchmark each ensemble size
    for nsims in ensemble_sizes
        println(styled"\n{bold green:Testing with $nsims simulations}")
        println("-"^50)

        # Create ensemble specifications
        ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(nsims)

        # Generate data
        println(styled"Generating ensemble data...")
        data_time = @elapsed data_arrs = generate_ensemble_data(
            ensemble_spec, null_spec, outbreak_spec
        )
        println(styled"  Data generation: {yellow:$(round(data_time, digits=2))}s")

        # Run multistart optimization
        println(styled"Running multistart optimization...")
        opt_time = @elapsed opt_results = ews_multistart_optimization(
            test_scenario,
            data_arrs;
            n_sobol_points = n_sobol_points,
            maxeval = 1000,
            executor = FLoops.ThreadedEx(),
            percentile_bounds = (0.5, 0.99),
            consecutive_bounds = (2.0, 30.0),
            force = true,
            return_df = true,
            save_results = false,
            verbose = false,
            disable_time_check = true
        )
        println(styled"  Optimization: {yellow:$(round(opt_time, digits=2))}s")

        # Store results
        best_accuracy = maximum(opt_results.accuracy)
        best_idx = argmax(opt_results.accuracy)
        best_params = opt_results[best_idx, :]

        results[nsims] = (
            data_time = data_time,
            opt_time = opt_time,
            total_time = data_time + opt_time,
            best_accuracy = best_accuracy,
            best_params = best_params,
            results_df = opt_results,
            nsims = nsims,
        )

        println(styled"  Best accuracy: {green:$(round(best_accuracy, digits=4))}")
        println(styled"  Total time: {yellow:$(round(data_time + opt_time, digits=2))}s")
    end

    # Compare results
    println(styled"\n{bold blue:COMPARISON ANALYSIS}")
    println("="^70)

    display_accuracy_comparison(results)
    display_performance_analysis(results)
    display_parameter_comparison(results)

    # Save results
    save_benchmark_comparison_results(results, "ensemble_size_comparison_$(n_sobol_points)sobol")

    return results
end

function create_single_test_scenario()
    """Create a single focused test scenario"""

    # Use moderate noise and good test characteristics
    noise_specification_vec = [PoissonNoiseSpecification(1.0)]
    test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)]
    percent_tested_vec = [1.0]

    # Focus on autocovariance metric
    ews_method_vec = [Backward]
    ews_aggregation_vec = [Day(28)]
    ews_bandwidth_vec = [Week(52)]
    ews_lag_days_vec = [1]

    ews_metric_specification_vec = create_combinations_vec(
        EWSMetricSpecification,
        (ews_method_vec, ews_aggregation_vec, ews_bandwidth_vec, ews_lag_days_vec)
    )

    ews_metric_vec = ["autocovariance"]
    ews_enddate_type_vec = [Reff_start]
    ews_threshold_window_vec = [ExpandingThresholdWindow]
    ews_threshold_percentile_vec = collect(0.5:0.02:0.99)
    ews_consecutive_thresholds_vec = collect(2:2:30)
    ews_threshold_burnin_vec = [Year(5)]

    return (;
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
end

# create_ensemble_specs function moved to src/benchmark-functions.jl

function display_accuracy_comparison(results)
    println(styled"\n{bold:Accuracy Comparison}")
    println("-"^40)

    acc_100 = results[100].best_accuracy
    acc_1000 = results[1000].best_accuracy
    acc_10000 = results[10000].best_accuracy

    println(styled"100 simulations:   {green:$(round(acc_100, digits=4))}")
    println(styled"1000 simulations:  {green:$(round(acc_1000, digits=4))}")
    println(styled"10000 simulations: {green:$(round(acc_10000, digits=4))}")

    # Compare improvements
    diff_1000_vs_100 = acc_1000 - acc_100
    diff_10000_vs_1000 = acc_10000 - acc_1000
    diff_10000_vs_100 = acc_10000 - acc_100

    rel_imp_1000 = (diff_1000_vs_100 / acc_100) * 100
    rel_imp_10000_vs_1000 = (diff_10000_vs_1000 / acc_1000) * 100
    rel_imp_10000_vs_100 = (diff_10000_vs_100 / acc_100) * 100

    println()
    println(styled"Improvements:")
    println(styled"  1000 vs 100:   {yellow:$(round(diff_1000_vs_100, digits=4))} ({cyan:$(round(rel_imp_1000, digits=2))%})")
    println(styled"  10000 vs 1000: {yellow:$(round(diff_10000_vs_1000, digits=4))} ({cyan:$(round(rel_imp_10000_vs_1000, digits=2))%})")
    println(styled"  10000 vs 100:  {yellow:$(round(diff_10000_vs_100, digits=4))} ({cyan:$(round(rel_imp_10000_vs_100, digits=2))%})")

    # Determine best performer
    best_acc = max(acc_100, acc_1000, acc_10000)
    if best_acc == acc_10000
        println(styled"\n{green:→ 10000 simulations provides the best accuracy}")
    elseif best_acc == acc_1000
        println(styled"\n{yellow:→ 1000 simulations provides the best accuracy}")
    else
        println(styled"\n{red:→ 100 simulations provides the best accuracy (possibly due to noise)}")
    end

    # Check for diminishing returns
    return if abs(diff_10000_vs_1000) < abs(diff_1000_vs_100) / 2
        println(styled"{blue:→ Diminishing returns: 1000→10000 improvement is less than half of 100→1000}")
    end
end

function display_performance_analysis(results)
    println(styled"\n{bold:Performance Analysis}")
    println("-"^40)

    time_100 = results[100].total_time
    time_1000 = results[1000].total_time
    time_10000 = results[10000].total_time

    time_per_sim_100 = time_100 / 100
    time_per_sim_1000 = time_1000 / 1000
    time_per_sim_10000 = time_10000 / 10000

    println(styled"Total Times:")
    println(styled"  100 simulations:   {yellow:$(round(time_100, digits=2))}s")
    println(styled"  1000 simulations:  {yellow:$(round(time_1000, digits=2))}s")
    println(styled"  10000 simulations: {yellow:$(round(time_10000, digits=2))}s")

    println(styled"\nTime per simulation:")
    println(styled"  100 simulations:   {cyan:$(round(time_per_sim_100*1000, digits=2))}ms")
    println(styled"  1000 simulations:  {cyan:$(round(time_per_sim_1000*1000, digits=2))}ms")
    println(styled"  10000 simulations: {cyan:$(round(time_per_sim_10000*1000, digits=2))}ms")

    println(styled"\nScaling factors:")
    scaling_1000_vs_100 = time_1000 / time_100
    scaling_10000_vs_1000 = time_10000 / time_1000
    scaling_10000_vs_100 = time_10000 / time_100

    println(styled"  1000 vs 100:   {yellow:$(round(scaling_1000_vs_100, digits=2))}x")
    println(styled"  10000 vs 1000: {yellow:$(round(scaling_10000_vs_1000, digits=2))}x")
    println(styled"  10000 vs 100:  {yellow:$(round(scaling_10000_vs_100, digits=2))}x")

    # Efficiency analysis
    efficiency_100 = results[100].best_accuracy / time_100
    efficiency_1000 = results[1000].best_accuracy / time_1000
    efficiency_10000 = results[10000].best_accuracy / time_10000

    println(styled"\nEfficiency (accuracy/second):")
    println(styled"  100 simulations:   {green:$(round(efficiency_100, digits=4))}")
    println(styled"  1000 simulations:  {green:$(round(efficiency_1000, digits=4))}")
    println(styled"  10000 simulations: {green:$(round(efficiency_10000, digits=4))}")

    # Find most efficient
    efficiencies = [efficiency_100, efficiency_1000, efficiency_10000]
    sizes = [100, 1000, 10000]
    best_eff_idx = argmax(efficiencies)
    best_size = sizes[best_eff_idx]

    println(styled"\n{blue:→ Most efficient: $best_size simulations}")

    # Check for linear scaling
    expected_linear_10000 = time_100 * 100  # Perfect linear scaling
    actual_vs_linear = time_10000 / expected_linear_10000

    return if actual_vs_linear < 1.2
        println(styled"{green:→ Excellent scaling: $(round(actual_vs_linear, digits=2))x of linear}")
    elseif actual_vs_linear < 2.0
        println(styled"{yellow:→ Good scaling: $(round(actual_vs_linear, digits=2))x of linear}")
    else
        println(styled"{red:→ Poor scaling: $(round(actual_vs_linear, digits=2))x of linear}")
    end
end

function display_parameter_comparison(results)
    println(styled"\n{bold:Parameter Comparison}")
    println("-"^40)

    params_100 = results[100].best_params
    params_1000 = results[1000].best_params
    params_10000 = results[10000].best_params

    println(styled"Threshold Percentile:")
    println(styled"  100 sims:   {cyan:$(round(params_100.ews_threshold_percentile, digits=3))}")
    println(styled"  1000 sims:  {cyan:$(round(params_1000.ews_threshold_percentile, digits=3))}")
    println(styled"  10000 sims: {cyan:$(round(params_10000.ews_threshold_percentile, digits=3))}")

    perc_diff_1000_100 = abs(params_100.ews_threshold_percentile - params_1000.ews_threshold_percentile)
    perc_diff_10000_1000 = abs(params_1000.ews_threshold_percentile - params_10000.ews_threshold_percentile)
    perc_diff_10000_100 = abs(params_100.ews_threshold_percentile - params_10000.ews_threshold_percentile)

    println(styled"  Differences:")
    println(styled"    1000 vs 100:   {yellow:$(round(perc_diff_1000_100, digits=3))}")
    println(styled"    10000 vs 1000: {yellow:$(round(perc_diff_10000_1000, digits=3))}")
    println(styled"    10000 vs 100:  {yellow:$(round(perc_diff_10000_100, digits=3))}")

    println(styled"\nConsecutive Thresholds:")
    println(styled"  100 sims:   {cyan:$(params_100.ews_consecutive_thresholds)}")
    println(styled"  1000 sims:  {cyan:$(params_1000.ews_consecutive_thresholds)}")
    println(styled"  10000 sims: {cyan:$(params_10000.ews_consecutive_thresholds)}")

    cons_diff_1000_100 = abs(params_100.ews_consecutive_thresholds - params_1000.ews_consecutive_thresholds)
    cons_diff_10000_1000 = abs(params_1000.ews_consecutive_thresholds - params_10000.ews_consecutive_thresholds)
    cons_diff_10000_100 = abs(params_100.ews_consecutive_thresholds - params_10000.ews_consecutive_thresholds)

    println(styled"  Differences:")
    println(styled"    1000 vs 100:   {yellow:$(cons_diff_1000_100)}")
    println(styled"    10000 vs 1000: {yellow:$(cons_diff_10000_1000)}")
    println(styled"    10000 vs 100:  {yellow:$(cons_diff_10000_100)}")

    # Overall convergence assessment
    max_perc_diff = max(perc_diff_1000_100, perc_diff_10000_1000, perc_diff_10000_100)
    max_cons_diff = max(cons_diff_1000_100, cons_diff_10000_1000, cons_diff_10000_100)

    println()
    if max_perc_diff < 0.01 && max_cons_diff <= 1
        println(styled"{green:→ Parameter estimates are very consistent across ensemble sizes}")
    elseif max_perc_diff < 0.05 && max_cons_diff <= 3
        println(styled"{yellow:→ Parameter estimates are reasonably consistent}")
    else
        println(styled"{red:→ Parameter estimates show significant variation with ensemble size}")
    end

    # Check for convergence pattern
    return if perc_diff_10000_1000 < perc_diff_1000_100 && cons_diff_10000_1000 < cons_diff_1000_100
        println(styled"{blue:→ Parameters are converging with larger ensemble sizes}")
    elseif perc_diff_10000_1000 > perc_diff_1000_100 || cons_diff_10000_1000 > cons_diff_1000_100
        println(styled"{orange:→ Parameter estimates may not be fully converged}")
    end
end

# save_benchmark_results function moved to src/benchmark-functions.jl as save_benchmark_comparison_results

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
