#!/usr/bin/env julia
"""
Profiling Script for EWS Multistart Optimization

This script profiles the multistart optimization function to identify performance 
bottlenecks and optimization opportunities. It uses multiple profiling approaches:
- BenchmarkTools for precise timing measurements
- PProf for CPU and memory profiling
- Built-in @time and @allocated for quick analysis
- Custom allocation tracking for specific functions

Usage:
    julia --project=. scripts/profile_multistart_optimization.jl
"""

using Pkg
Pkg.activate(".")

# Core dependencies
using CSDNoise
using BenchmarkTools
using PProf
using ProfileCanvas
using ProfileView
using DataFrames
using Statistics
using Printf
using Dates
using NLopt

# Profiling utilities
using InteractiveUtils: @code_warntype
using Base.GC: gc

println("🔍 EWS Multistart Optimization Profiler")
println("="^50)

# =============================================================================
# SETUP TEST DATA AND SCENARIOS
# =============================================================================

"""
Create minimal test data for profiling purposes by generating simulations.
"""
function create_test_data_for_profiling()
    println("📊 Creating test data for profiling...")

    # Create minimal but realistic test scenarios (similar to benchmark_ensemble_size.jl)
    specification_vecs = create_single_test_scenario()

    # Generate small ensemble data for profiling (use small number of simulations)
    nsims = 50  # Small for profiling
    println("  Generating ensemble specifications for $nsims simulations...")

    ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(nsims)

    println("  Generating ensemble data...")
    data_arrs = generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)

    println("✅ Generated test data for profiling")

    return specification_vecs, data_arrs
end

"""
Create a single focused test scenario for profiling (adapted from benchmark_ensemble_size.jl).
"""
function create_single_test_scenario()
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

# =============================================================================
# INDIVIDUAL FUNCTION BENCHMARKS
# =============================================================================

"""
Benchmark individual functions that are likely bottlenecks.
"""
function benchmark_individual_functions(specification_vecs, data_arrs)
    println("\n🎯 Benchmarking Individual Functions")
    println("-"^40)

    # Create a single scenario for testing
    scenarios_vec = create_optimization_scenarios(specification_vecs)
    scenario = scenarios_vec[1]

    println("Function: create_cached_simulation_data")
    bench_cache = @benchmark create_cached_simulation_data($scenario, $data_arrs)
    println("  Time: $(BenchmarkTools.prettytime(median(bench_cache.times)))")
    println("  Memory: $(BenchmarkTools.prettymemory(median(bench_cache.memory)))")
    println("  Allocs: $(median(bench_cache.allocs))")

    # Cache the data for subsequent benchmarks
    cached_data = create_cached_simulation_data(scenario, data_arrs)

    # Benchmark the objective function
    bounds = (lowers = [0.5, 2.0], uppers = [0.99, 30.0])
    config = (
        n_sobol_points = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
    )

    tracker = OptimizationTracker()
    test_params = [0.9, 5.0]  # Reasonable parameters

    println("\nFunction: ews_objective_function_with_tracking")
    bench_objective = @benchmark ews_objective_function_with_tracking(
        $test_params, $scenario, $cached_data, $tracker
    )
    println("  Time: $(BenchmarkTools.prettytime(median(bench_objective.times)))")
    println("  Memory: $(BenchmarkTools.prettymemory(median(bench_objective.memory)))")
    println("  Allocs: $(median(bench_objective.allocs))")

    println("\nFunction: optimize_single_scenario")
    bench_single = @benchmark optimize_single_scenario($scenario, $data_arrs, $bounds, $config)
    println("  Time: $(BenchmarkTools.prettytime(median(bench_single.times)))")
    println("  Memory: $(BenchmarkTools.prettymemory(median(bench_single.memory)))")
    println("  Allocs: $(median(bench_single.allocs))")

    return (
        cache_bench = bench_cache,
        objective_bench = bench_objective,
        single_scenario_bench = bench_single,
    )
end

# =============================================================================
# CPU PROFILING WITH PPROF
# =============================================================================

"""
Profile CPU usage using PProf for detailed analysis.
"""
function profile_with_pprof(specification_vecs, data_arrs)
    println("\n🔥 CPU Profiling with PProf")
    println("-"^40)

    # Reduce problem size for profiling
    small_spec_vecs = (
        noise_specification_vec = specification_vecs.noise_specification_vec[1:1],
        test_specification_vec = specification_vecs.test_specification_vec[1:1],
        percent_tested_vec = specification_vecs.percent_tested_vec[1:1],
        ews_metric_specification_vec = specification_vecs.ews_metric_specification_vec[1:1],
        ews_enddate_type_vec = specification_vecs.ews_enddate_type_vec[1:1],
        ews_threshold_window_vec = specification_vecs.ews_threshold_window_vec[1:1],
        ews_threshold_burnin_vec = specification_vecs.ews_threshold_burnin_vec[1:1],
        ews_metric_vec = specification_vecs.ews_metric_vec[1:1],
    )

    println("Starting CPU profiling...")

    # Profile the main optimization function
    PProf.@pprof begin
        result = ews_multistart_optimization(
            small_spec_vecs,
            data_arrs;
            n_sobol_points = 100,
            maxeval = 1000,
            batch_size = 10,
            save_results = false,    # Don't save during profiling
            verbose = false,         # Reduce output noise
            force = true,             # Skip existing results check
            disable_time_check = true
        )
    end

    println("✅ CPU profiling complete. Check the generated flamegraph.")
    return println("   Open the HTML file that was created to view detailed CPU usage.")
end

# =============================================================================
# TYPE STABILITY ANALYSIS
# =============================================================================

"""
Check type stability of key functions.
"""
function analyze_type_stability(specification_vecs, data_arrs)
    println("\n🔍 Type Stability Analysis")
    println("-"^40)

    scenarios_vec = create_optimization_scenarios(specification_vecs)
    scenario = scenarios_vec[1]
    cached_data = create_cached_simulation_data(scenario, data_arrs)
    tracker = OptimizationTracker()
    test_params = [0.9, 5.0]

    println("Type analysis for ews_objective_function_with_tracking:")
    @code_warntype ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker)

    println("\nType analysis for create_cached_simulation_data:")
    return @code_warntype create_cached_simulation_data(scenario, data_arrs)
end

# =============================================================================
# MAIN PROFILING EXECUTION
# =============================================================================

function main()
    println("Starting comprehensive profiling analysis...")
    println("Timestamp: $(Dates.now())")

    # Setup
    specification_vecs, data_arrs = create_test_data_for_profiling()

    if isnothing(specification_vecs) || isnothing(data_arrs)
        println("❌ Could not create test data. Please check data loading.")
        return
    end

    return try
        # Run individual benchmarks
        individual_benchmarks = benchmark_individual_functions(specification_vecs, data_arrs)

        # CPU profiling (comment out if you don't want the interactive flamegraph)
        # profile_with_pprof(specification_vecs, data_arrs)

        # Type stability analysis
        analyze_type_stability(specification_vecs, data_arrs)

        println("\n✅ Profiling analysis complete!")
        println("\n📋 Summary of Recommendations:")
        println("1. Check the individual function benchmarks for the biggest time consumers")
        println("2. Review memory allocation patterns for optimization opportunities")
        println("3. Consider the performance/accuracy tradeoffs in the configuration comparison")
        println("4. Address any type instabilities found in the analysis")
        println("5. Run PProf profiling (uncomment the line) for detailed CPU analysis")

    catch e
        println("❌ Error during profiling: $e")
        println("Stack trace:")
        showerror(stdout, e, catch_backtrace())
    end
end

# Run the profiling if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
