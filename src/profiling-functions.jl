using BenchmarkTools
using PProf
using ProfileCanvas
using ProfileView
using JET
using StatsBase
using Dates
using NLopt
using InteractiveUtils: @code_warntype

export create_test_data_for_profiling, create_single_test_scenario, benchmark_individual_functions, profile_with_pprof, analyze_type_stability

# =============================================================================
# SETUP TEST DATA AND SCENARIOS
# =============================================================================

"""
Create minimal test data for profiling purposes by generating simulations.
"""
@unstable function create_test_data_for_profiling()
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
    noise_specification_vec = NoiseSpecification[NoiseSpecification(PoissonNoise(1.0))]
    test_specification_vec = IndividualTestSpecification[IndividualTestSpecification(0.9, 0.9, 0)]
    percent_tested_vec = Float64[1.0]

    # Focus on autocovariance metric
    ews_method_vec = [EWSMethod(Backward())]
    ews_aggregation_vec = Dates.Period[Day(28)]
    ews_bandwidth_vec = Dates.Period[Week(52)]
    ews_lag_days_vec = Int64[1]

    ews_metric_specification_vec = create_combinations_vec(
        EWSMetricSpecification,
        (ews_method_vec, ews_aggregation_vec, ews_bandwidth_vec, ews_lag_days_vec)
    )

    ews_metric_vec = ["autocovariance"]
    ews_enddate_type_vec = [EWSEndDateType(Reff_start())]
    ews_threshold_window_vec = [EWSThresholdWindowType(ExpandingThresholdWindow())]
    ews_threshold_burnin_vec = Dates.Period[Year(5)]

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

    println("\nType analysis for create_cached_simulation_data:")
    # @code_warntype create_cached_simulation_data(scenario, data_arrs)
    println("\n\tDynamics dispatch check:")
    @report_opt create_cached_simulation_data(scenario, data_arrs)
    println("\n\tType error check:")
    @report_call target_modules = (CSDNoise,) create_cached_simulation_data(scenario, data_arrs)
    cached_data = create_cached_simulation_data(scenario, data_arrs)

    tracker = OptimizationTracker()
    test_params = [0.9, 5.0]

    println("Type analysis for ews_objective_function_with_tracking:")
    @code_warntype ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker)
    println("\n\tDynamics dispatch check:")
    @report_opt target_modules = (CSDNoise,) ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker)
    println("\n\tType error check:")
    @report_call target_modules = (CSDNoise,) ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker)
    return nothing
end
