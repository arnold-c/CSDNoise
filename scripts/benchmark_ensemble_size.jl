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
using StatsBase
using CSV
using FLoops
using UnPack
using Logging
using LoggingExtras

# Import logging utilities
using CSDNoise: styled_to_markdown, format_markdown_table, setup_dual_logging, cleanup_logging, log_both

function main()
    # Set seed for reproducibility
    Random.seed!(42)

    # Define ensemble sizes to compare
    ensemble_sizes = [100] #, 1000, 10000]
    n_sobol_points = 100

    # Create consistent filename base for both CSV and markdown files
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    filename_base = "ensemble-size-comparison_$(n_sobol_points)_sobol_$(timestamp)"

    # Setup dual logging in the same directory as CSV results
    benchmark_dir = outdir("benchmark")
    markdown_filename = setup_dual_logging("benchmark-ensemble-size"; title = "Ensemble Size Benchmark Report", output_dir = benchmark_dir, filename_base = filename_base)

    try
        log_both(styled"{bold blue:Ensemble Size Benchmark (100 vs 1000 vs 10000 simulations)}")
        log_both("="^70)

        # Define single test scenario (autocovariance with moderate noise)
        test_scenario = create_single_test_scenario()

        log_both(styled"Configuration:")
        log_both(styled"  Scenario: {cyan:Autocovariance, Poisson(1.0x), Test(0.9,0.9)}")
        log_both(styled"  Sobol points: {cyan:$n_sobol_points}")
        ensemble_sizes_str = join(ensemble_sizes, ", ")
        log_both(styled"  Ensemble sizes: {cyan:$ensemble_sizes_str}")
        log_both("")

        results = Dict()

        # Benchmark each ensemble size
        for nsims in ensemble_sizes
            log_both(styled"\n{bold green:Testing with $nsims simulations}")
            log_both("-"^50)

            # Create ensemble specifications
            ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(nsims)

            # Generate data
            log_both(styled"Generating ensemble data...")
            data_time = @elapsed data_arrs = generate_ensemble_data(
                ensemble_spec, null_spec, outbreak_spec
            )
            log_both(styled"  Data generation: {yellow:$(round(data_time, digits=2))}s")

            # Run multistart optimization
            log_both(styled"Running multistart optimization...")
            opt_time = @elapsed opt_results = ews_multistart_optimization(
                test_scenario,
                data_arrs;
                n_sobol_points = n_sobol_points,
                maxeval = 1000,
                executor = FLoops.SequentialEx(),
                percentile_bounds = (0.5, 0.99),
                consecutive_bounds = (2.0, 30.0),
                force = true,
                return_df = true,
                save_results = false,
                verbose = false,
                disable_time_check = true
            )
            log_both(styled"  Optimization: {yellow:$(round(opt_time, digits=2))}s")

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

            log_both(styled"  Best accuracy: {green:$(round(best_accuracy, digits=4))}")
            log_both(styled"  Total time: {yellow:$(round(data_time + opt_time, digits=2))}s")
        end

        # Compare results
        log_both(styled"\n{bold blue:COMPARISON ANALYSIS}")
        log_both("="^70)

        # display_accuracy_comparison(results)
        # display_performance_analysis(results)
        # display_parameter_comparison(results)

        # Save results using the same filename base
        save_benchmark_comparison_results(results, "ensemble-size-comparison_$(n_sobol_points)_sobol"; full_filename = filename_base)

        @info ""
        @info "---"
        @info ""
        @info "**Report saved to:** `$markdown_filename`"

        return results
    finally
        cleanup_logging()
    end
end

function create_single_test_scenario()
    """Create a single focused test scenario"""

    # Use moderate noise and good test characteristics
    noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))]
    test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)]
    percent_tested_vec = [1.0]

    # Focus on autocovariance metric
    ews_method_vec = [EWSMethod(Backward())]
    ews_aggregation_vec = [Day(28)]
    ews_bandwidth_vec = [Week(52)]
    ews_lag_days_vec = [1]

    ews_metric_specification_vec = create_combinations_vec(
        EWSMetricSpecification,
        (ews_method_vec, ews_aggregation_vec, ews_bandwidth_vec, ews_lag_days_vec)
    )

    ews_metric_vec = ["autocovariance"]
    ews_enddate_type_vec = [EWSEndDateType(Reff_start())]
    ews_threshold_window_vec = [EWSThresholdWindowType(ExpandingThresholdWindow())]
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
        ews_metric_vec,
    )
end

# create_ensemble_specs function moved to src/benchmark-functions.jl

function display_accuracy_comparison(results)
    log_both(
        styled"\n{bold:Accuracy Comparison}",
        "### Accuracy Comparison"
    )
    log_both("-"^40, "")

    acc_100 = results[100].best_accuracy
    acc_1000 = results[1000].best_accuracy
    acc_10000 = results[10000].best_accuracy

    # Create accuracy table for markdown
    accuracy_headers = ["Ensemble Size", "Best Accuracy"]
    accuracy_rows = [
        ["100 simulations", "$(round(acc_100, digits = 4))"],
        ["1000 simulations", "$(round(acc_1000, digits = 4))"],
        ["10000 simulations", "$(round(acc_10000, digits = 4))"],
    ]
    accuracy_table = format_markdown_table(accuracy_headers, accuracy_rows)

    log_both(styled"100 simulations:   {green:$(round(acc_100, digits=4))}", "")
    log_both(styled"1000 simulations:  {green:$(round(acc_1000, digits=4))}", "")
    log_both(styled"10000 simulations: {green:$(round(acc_10000, digits=4))}", "")
    @info accuracy_table

    # Compare improvements
    diff_1000_vs_100 = acc_1000 - acc_100
    diff_10000_vs_1000 = acc_10000 - acc_1000
    diff_10000_vs_100 = acc_10000 - acc_100

    rel_imp_1000 = (diff_1000_vs_100 / acc_100) * 100
    rel_imp_10000_vs_1000 = (diff_10000_vs_1000 / acc_1000) * 100
    rel_imp_10000_vs_100 = (diff_10000_vs_100 / acc_100) * 100

    log_both("", "")
    log_both(styled"Improvements:", "#### Improvements")

    # Create improvements table for markdown
    improvements_headers = ["Comparison", "Absolute Difference", "Relative Improvement (%)"]
    improvements_rows = [
        ["1000 vs 100", "$(round(diff_1000_vs_100, digits = 4))", "$(round(rel_imp_1000, digits = 2))%"],
        ["10000 vs 1000", "$(round(diff_10000_vs_1000, digits = 4))", "$(round(rel_imp_10000_vs_1000, digits = 2))%"],
        ["10000 vs 100", "$(round(diff_10000_vs_100, digits = 4))", "$(round(rel_imp_10000_vs_100, digits = 2))%"],
    ]
    improvements_table = format_markdown_table(improvements_headers, improvements_rows)

    log_both(styled"  1000 vs 100:   {yellow:$(round(diff_1000_vs_100, digits=4))} ({cyan:$(round(rel_imp_1000, digits=2))%})", "")
    log_both(styled"  10000 vs 1000: {yellow:$(round(diff_10000_vs_1000, digits=4))} ({cyan:$(round(rel_imp_10000_vs_1000, digits=2))%})", "")
    log_both(styled"  10000 vs 100:  {yellow:$(round(diff_10000_vs_100, digits=4))} ({cyan:$(round(rel_imp_10000_vs_100, digits=2))%})", "")
    @info improvements_table

    # Determine best performer
    best_acc = max(acc_100, acc_1000, acc_10000)
    if best_acc == acc_10000
        log_both(
            styled"\n{green:→ 10000 simulations provides the best accuracy}",
            "**→ 10000 simulations provides the best accuracy**"
        )
    elseif best_acc == acc_1000
        log_both(
            styled"\n{yellow:→ 1000 simulations provides the best accuracy}",
            "**→ 1000 simulations provides the best accuracy**"
        )
    else
        log_both(
            styled"\n{red:→ 100 simulations provides the best accuracy (possibly due to noise)}",
            "**→ 100 simulations provides the best accuracy (possibly due to noise)**"
        )
    end

    # Check for diminishing returns
    return if abs(diff_10000_vs_1000) < abs(diff_1000_vs_100) / 2
        log_both(
            styled"{blue:→ Diminishing returns: 1000→10000 improvement is less than half of 100→1000}",
            "**→ Diminishing returns:** 1000→10000 improvement is less than half of 100→1000"
        )
    end
end

function display_performance_analysis(results)
    log_both(
        styled"\n{bold:Performance Analysis}",
        "### Performance Analysis"
    )
    log_both("-"^40, "")

    time_100 = results[100].total_time
    time_1000 = results[1000].total_time
    time_10000 = results[10000].total_time

    time_per_sim_100 = time_100 / 100
    time_per_sim_1000 = time_1000 / 1000
    time_per_sim_10000 = time_10000 / 10000

    # Create performance table for markdown
    perf_headers = ["Ensemble Size", "Total Time (s)", "Time per Simulation (ms)"]
    perf_rows = [
        ["100 simulations", "$(round(time_100, digits = 2))", "$(round(time_per_sim_100 * 1000, digits = 2))"],
        ["1000 simulations", "$(round(time_1000, digits = 2))", "$(round(time_per_sim_1000 * 1000, digits = 2))"],
        ["10000 simulations", "$(round(time_10000, digits = 2))", "$(round(time_per_sim_10000 * 1000, digits = 2))"],
    ]
    perf_table = format_markdown_table(perf_headers, perf_rows)

    log_both(styled"Total Times:", "#### Total Times")
    log_both(styled"  100 simulations:   {yellow:$(round(time_100, digits=2))}s", "")
    log_both(styled"  1000 simulations:  {yellow:$(round(time_1000, digits=2))}s", "")
    log_both(styled"  10000 simulations: {yellow:$(round(time_10000, digits=2))}s", "")

    log_both(styled"\nTime per simulation:", "#### Time per Simulation")
    log_both(styled"  100 simulations:   {cyan:$(round(time_per_sim_100*1000, digits=2))}ms", "")
    log_both(styled"  1000 simulations:  {cyan:$(round(time_per_sim_1000*1000, digits=2))}ms", "")
    log_both(styled"  10000 simulations: {cyan:$(round(time_per_sim_10000*1000, digits=2))}ms", "")
    @info perf_table

    scaling_1000_vs_100 = time_1000 / time_100
    scaling_10000_vs_1000 = time_10000 / time_1000
    scaling_10000_vs_100 = time_10000 / time_100

    # Create scaling table for markdown
    scaling_headers = ["Comparison", "Scaling Factor"]
    scaling_rows = [
        ["1000 vs 100", "$(round(scaling_1000_vs_100, digits = 2))x"],
        ["10000 vs 1000", "$(round(scaling_10000_vs_1000, digits = 2))x"],
        ["10000 vs 100", "$(round(scaling_10000_vs_100, digits = 2))x"],
    ]
    scaling_table = format_markdown_table(scaling_headers, scaling_rows)

    log_both(styled"\nScaling factors:", "#### Scaling Factors")
    log_both(styled"  1000 vs 100:   {yellow:$(round(scaling_1000_vs_100, digits=2))}x", "")
    log_both(styled"  10000 vs 1000: {yellow:$(round(scaling_10000_vs_1000, digits=2))}x", "")
    log_both(styled"  10000 vs 100:  {yellow:$(round(scaling_10000_vs_100, digits=2))}x", "")
    @info scaling_table

    # Efficiency analysis
    efficiency_100 = results[100].best_accuracy / time_100
    efficiency_1000 = results[1000].best_accuracy / time_1000
    efficiency_10000 = results[10000].best_accuracy / time_10000

    # Create efficiency table for markdown
    eff_headers = ["Ensemble Size", "Efficiency (accuracy/second)"]
    eff_rows = [
        ["100 simulations", "$(round(efficiency_100, digits = 4))"],
        ["1000 simulations", "$(round(efficiency_1000, digits = 4))"],
        ["10000 simulations", "$(round(efficiency_10000, digits = 4))"],
    ]
    eff_table = format_markdown_table(eff_headers, eff_rows)

    log_both(styled"\nEfficiency (accuracy/second):", "#### Efficiency (accuracy/second)")
    log_both(styled"  100 simulations:   {green:$(round(efficiency_100, digits=4))}", "")
    log_both(styled"  1000 simulations:  {green:$(round(efficiency_1000, digits=4))}", "")
    log_both(styled"  10000 simulations: {green:$(round(efficiency_10000, digits=4))}", "")
    @info eff_table

    # Find most efficient
    efficiencies = [efficiency_100, efficiency_1000, efficiency_10000]
    sizes = [100, 1000, 10000]
    best_eff_idx = argmax(efficiencies)
    best_size = sizes[best_eff_idx]

    log_both(
        styled"\n{blue:→ Most efficient: $best_size simulations}",
        "**→ Most efficient:** $best_size simulations"
    )

    # Check for linear scaling
    expected_linear_10000 = time_100 * 100  # Perfect linear scaling
    actual_vs_linear = time_10000 / expected_linear_10000

    return if actual_vs_linear < 1.2
        log_both(
            styled"{green:→ Excellent scaling: $(round(actual_vs_linear, digits=2))x of linear}",
            "**→ Excellent scaling:** $(round(actual_vs_linear, digits = 2))x of linear"
        )
    elseif actual_vs_linear < 2.0
        log_both(
            styled"{yellow:→ Good scaling: $(round(actual_vs_linear, digits=2))x of linear}",
            "**→ Good scaling:** $(round(actual_vs_linear, digits = 2))x of linear"
        )
    else
        log_both(
            styled"{red:→ Poor scaling: $(round(actual_vs_linear, digits=2))x of linear}",
            "**→ Poor scaling:** $(round(actual_vs_linear, digits = 2))x of linear"
        )
    end
end

function display_parameter_comparison(results)
    log_both(
        styled"\n{bold:Parameter Comparison}",
        "### Parameter Comparison"
    )
    log_both("-"^40, "")

    params_100 = results[100].best_params
    params_1000 = results[1000].best_params
    params_10000 = results[10000].best_params

    # Create parameter comparison table for markdown
    param_headers = ["Ensemble Size", "Threshold Percentile", "Consecutive Thresholds"]
    param_rows = [
        ["100 sims", "$(round(params_100.ews_threshold_percentile, digits = 3))", "$(params_100.ews_consecutive_thresholds)"],
        ["1000 sims", "$(round(params_1000.ews_threshold_percentile, digits = 3))", "$(params_1000.ews_consecutive_thresholds)"],
        ["10000 sims", "$(round(params_10000.ews_threshold_percentile, digits = 3))", "$(params_10000.ews_consecutive_thresholds)"],
    ]
    param_table = format_markdown_table(param_headers, param_rows)

    log_both(styled"Threshold Percentile:", "#### Threshold Percentile")
    log_both(styled"  100 sims:   {cyan:$(round(params_100.ews_threshold_percentile, digits=3))}", "")
    log_both(styled"  1000 sims:  {cyan:$(round(params_1000.ews_threshold_percentile, digits=3))}", "")
    log_both(styled"  10000 sims: {cyan:$(round(params_10000.ews_threshold_percentile, digits=3))}", "")

    perc_diff_1000_100 = abs(params_100.ews_threshold_percentile - params_1000.ews_threshold_percentile)
    perc_diff_10000_1000 = abs(params_1000.ews_threshold_percentile - params_10000.ews_threshold_percentile)
    perc_diff_10000_100 = abs(params_100.ews_threshold_percentile - params_10000.ews_threshold_percentile)

    log_both(styled"  Differences:", "**Percentile Differences:**")
    log_both(styled"    1000 vs 100:   {yellow:$(round(perc_diff_1000_100, digits=3))}", "- 1000 vs 100: *$(round(perc_diff_1000_100, digits = 3))*")
    log_both(styled"    10000 vs 1000: {yellow:$(round(perc_diff_10000_1000, digits=3))}", "- 10000 vs 1000: *$(round(perc_diff_10000_1000, digits = 3))*")
    log_both(styled"    10000 vs 100:  {yellow:$(round(perc_diff_10000_100, digits=3))}", "- 10000 vs 100: *$(round(perc_diff_10000_100, digits = 3))*")

    log_both(styled"\nConsecutive Thresholds:", "#### Consecutive Thresholds")
    log_both(styled"  100 sims:   {cyan:$(params_100.ews_consecutive_thresholds)}", "")
    log_both(styled"  1000 sims:  {cyan:$(params_1000.ews_consecutive_thresholds)}", "")
    log_both(styled"  10000 sims: {cyan:$(params_10000.ews_consecutive_thresholds)}", "")

    cons_diff_1000_100 = abs(params_100.ews_consecutive_thresholds - params_1000.ews_consecutive_thresholds)
    cons_diff_10000_1000 = abs(params_1000.ews_consecutive_thresholds - params_10000.ews_consecutive_thresholds)
    cons_diff_10000_100 = abs(params_100.ews_consecutive_thresholds - params_10000.ews_consecutive_thresholds)

    log_both(styled"  Differences:", "**Consecutive Threshold Differences:**")
    log_both(styled"    1000 vs 100:   {yellow:$(cons_diff_1000_100)}", "- 1000 vs 100: *$(cons_diff_1000_100)*")
    log_both(styled"    10000 vs 1000: {yellow:$(cons_diff_10000_1000)}", "- 10000 vs 1000: *$(cons_diff_10000_1000)*")
    log_both(styled"    10000 vs 100:  {yellow:$(cons_diff_10000_100)}", "- 10000 vs 100: *$(cons_diff_10000_100)*")

    @info param_table

    # Overall convergence assessment
    max_perc_diff = max(perc_diff_1000_100, perc_diff_10000_1000, perc_diff_10000_100)
    max_cons_diff = max(cons_diff_1000_100, cons_diff_10000_1000, cons_diff_10000_100)

    log_both("", "")
    if max_perc_diff < 0.01 && max_cons_diff <= 1
        log_both(
            styled"{green:→ Parameter estimates are very consistent across ensemble sizes}",
            "**→ Parameter estimates are very consistent across ensemble sizes**"
        )
    elseif max_perc_diff < 0.05 && max_cons_diff <= 3
        log_both(
            styled"{yellow:→ Parameter estimates are reasonably consistent}",
            "**→ Parameter estimates are reasonably consistent**"
        )
    else
        log_both(
            styled"{red:→ Parameter estimates show significant variation with ensemble size}",
            "**→ Parameter estimates show significant variation with ensemble size**"
        )
    end

    # Check for convergence pattern
    return if perc_diff_10000_1000 < perc_diff_1000_100 && cons_diff_10000_1000 < cons_diff_1000_100
        log_both(
            styled"{blue:→ Parameters are converging with larger ensemble sizes}",
            "**→ Parameters are converging with larger ensemble sizes**"
        )
    elseif perc_diff_10000_1000 > perc_diff_1000_100 || cons_diff_10000_1000 > cons_diff_1000_100
        log_both(
            styled"{orange:→ Parameter estimates may not be fully converged}",
            "**→ Parameter estimates may not be fully converged**"
        )
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
