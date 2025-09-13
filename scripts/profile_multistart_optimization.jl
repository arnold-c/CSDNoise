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

using DrWatson
@quickactivate("CSDNoise")

# Core dependencies
using CSDNoise
using JET
using Dates

println("🔍 EWS Multistart Optimization Profiler")
println("="^50)

# =============================================================================
# MAIN PROFILING EXECUTION
# =============================================================================

function main()
    println("Starting comprehensive profiling analysis...")
    println("Timestamp: $(Dates.now())")

    # Setup
    specification_vecs, data_arrs = create_test_data_for_profiling()

    scenarios_vec = create_optimization_scenarios(specification_vecs)

    scenario = scenarios_vec[1]
    cached_data = create_cached_simulation_data(scenario, data_arrs)

    tracker = OptimizationTracker()
    test_params = [0.9, 5.0]

    # Example EWSMetrics instantiation
    example_ews_spec = EWSMetricSpecification(EWSMethod(Backward()), Dates.Day(1), Dates.Day(30), 1)
    example_timeseries = [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
        21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
        31.0, 32.0, 33.0, 34.0, 35.0,
    ]

    if isnothing(specification_vecs) || isnothing(data_arrs)
        println("❌ Could not create test data. Please check data loading.")
        return
    end

    return try
        # Run individual benchmarks
        individual_benchmarks = benchmark_individual_functions(specification_vecs, data_arrs)

        # CPU profiling (comment out if you don't want the interactive flamegraph)
        # profile_with_pprof(specification_vecs, data_arrs)

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
