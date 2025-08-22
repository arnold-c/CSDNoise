using CSDNoise
using DataFrames
using Dates
using Random
using Printf
using StyledStrings

"""
Test script to compare multistart optimization vs grid search
"""
function benchmark_optimization_methods()
    println("="^70)
    println("BENCHMARKING: Multistart Optimization vs Grid Search")
    println("="^70)

    # Define small test space for fair comparison
    specification_vecs = (
        noise_specification_vec = [
            PoissonNoiseSpecification(1.0),
            PoissonNoiseSpecification(3.0),
        ],
        test_specification_vec = [
            IndividualTestSpecification(0.8, 0.9, 0),
            IndividualTestSpecification(1.0, 1.0, 0),
        ],
        percent_tested_vec = [0.1, 0.3],
        ews_metric_specification_vec = [
            EWSMetricSpecification(Gaussian, Day(7), Day(28), 1),
        ],
        ews_enddate_type_vec = [Reff_start],
        ews_threshold_window_vec = [ExpandingThresholdWindow],
        ews_metric_vec = ["mean", "variance"],
        # Grid search specific parameters
        ews_threshold_percentile_vec = [0.5, 0.7, 0.9, 0.95],
        ews_consecutive_thresholds_vec = [1, 3, 5, 7, 10],
        ews_threshold_burnin_vec = [Day(7), Day(30), Day(90), Day(180)],
    )

    # Generate test data
    println("\nGenerating test data...")
    Random.seed!(42)

    # Create test ensemble (small for benchmarking)
    ensemble_specification = EnsembleSpecification(
        StateParameters(
            SVector(990, 10, 0, 0),  # S, E, I, R
            0.01,  # p_infection
        ),
        DynamicsParameters(
            0.4,   # β₀
            0.1,   # σ
            0.1,   # γ
            0.0,   # μ
        ),
        TimeParameters(
            Day(1),    # tstep
            500,       # tlength
            Date(2020, 1, 1):Day(1):Date(2021, 4, 15),  # trange
        ),
        50,  # nsims - small for testing
    )

    # Generate test data arrays
    println("Creating ensemble simulations...")
    data_arrs = create_test_data_arrays(ensemble_specification)

    # Calculate expected grid search evaluations
    n_scenarios = length(specification_vecs.noise_specification_vec) *
        length(specification_vecs.test_specification_vec) *
        length(specification_vecs.percent_tested_vec) *
        length(specification_vecs.ews_metric_specification_vec) *
        length(specification_vecs.ews_enddate_type_vec) *
        length(specification_vecs.ews_threshold_window_vec) *
        length(specification_vecs.ews_metric_vec)

    n_grid_points = length(specification_vecs.ews_threshold_percentile_vec) *
        length(specification_vecs.ews_consecutive_thresholds_vec) *
        length(specification_vecs.ews_threshold_burnin_vec)

    total_grid_evaluations = n_scenarios * n_grid_points

    println("\nScenarios: $n_scenarios")
    println("Grid points per scenario: $n_grid_points")
    println("Total grid search evaluations: $total_grid_evaluations")

    # Benchmark grid search
    println("\n" * "=" * 70)
    println("GRID SEARCH")
    println("=" * 70)

    grid_time = @elapsed grid_results = ews_hyperparam_optimization(
        specification_vecs,
        data_arrs;
        disable_time_check = true,
        force = true,
        return_df = true,
    )

    println("Time: $(round(grid_time, digits = 2)) seconds")
    println("Results: $(nrow(grid_results)) parameter combinations")
    println("Best accuracy: $(round(maximum(grid_results.accuracy), digits = 4))")

    # Benchmark multistart with different Sobol point counts
    sobol_configs = [20, 50, 100]

    println("\n" * "=" * 70)
    println("MULTISTART OPTIMIZATION")
    println("=" * 70)

    multistart_results = Dict()

    for n_sobol in sobol_configs
        println("\nSobol points: $n_sobol")
        total_multistart_evaluations = n_scenarios * n_sobol
        println("Total evaluations: $total_multistart_evaluations")
        println("Reduction vs grid: $(round(100 * (1 - total_multistart_evaluations / total_grid_evaluations), digits = 1))%")

        multistart_time = @elapsed results = ews_multistart_optimization(
            specification_vecs,
            data_arrs;
            n_sobol_points = n_sobol,
            maxeval = n_sobol * 10,
            use_threads = false,
            force = true,
            return_df = true,
            verbose = false,
        )

        multistart_results[n_sobol] = results

        println("Time: $(round(multistart_time, digits = 2)) seconds")
        println("Speedup: $(round(grid_time / multistart_time, digits = 1))x")
        println("Best accuracy: $(round(maximum(results.accuracy), digits = 4))")
    end

    # Compare best parameters found
    println("\n" * "=" * 70)
    println("PARAMETER COMPARISON")
    println("=" * 70)

    # Find best configuration from grid search
    best_grid = grid_results[argmax(grid_results.accuracy), :]

    println("\nGrid Search Best:")
    println("  Accuracy: $(round(best_grid.accuracy, digits = 4))")
    println("  Percentile: $(best_grid.ews_threshold_percentile)")
    println("  Consecutive: $(best_grid.ews_consecutive_thresholds)")
    println("  Burnin: $(best_grid.ews_threshold_burnin)")

    # Compare with multistart results
    for n_sobol in sobol_configs
        results = multistart_results[n_sobol]
        best_ms = results[argmax(results.accuracy), :]

        println("\nMultistart ($n_sobol points) Best:")
        println("  Accuracy: $(round(best_ms.accuracy, digits = 4))")
        println("  Percentile: $(round(best_ms.ews_threshold_percentile, digits = 3))")
        println("  Consecutive: $(best_ms.ews_consecutive_thresholds)")
        println("  Burnin: $(best_ms.ews_threshold_burnin)")

        accuracy_diff = best_ms.accuracy - best_grid.accuracy
        println("  Δ Accuracy: $(round(accuracy_diff, digits = 4))")
    end

    return grid_results, multistart_results
end

"""
Create test data arrays for benchmarking
"""
function create_test_data_arrays(ensemble_specification)
    # Generate ensemble simulations
    @unpack state_parameters, dynamics_parameters, time_parameters, nsims = ensemble_specification
    @unpack tstep, tlength, trange = time_parameters

    # Initialize arrays
    ensemble_single_incarr = zeros(Float64, tlength, 5, nsims)  # [time, variables, sims]
    null_single_incarr = zeros(Float64, tlength, 5, nsims)

    # Generate simple test data (placeholder - replace with actual simulation)
    Random.seed!(42)
    for sim in 1:nsims
        # Generate synthetic incidence data
        base_incidence = 5.0 + 2.0 * randn()
        trend = 0.01 * randn()

        for t in 1:tlength
            # Simple synthetic data with trend and noise
            incidence = max(0.0, base_incidence + trend * t + 0.5 * randn())
            ensemble_single_incarr[t, 5, sim] = incidence  # Column 5 is incidence

            # Null data (no trend)
            null_incidence = max(0.0, base_incidence + 0.3 * randn())
            null_single_incarr[t, 5, sim] = null_incidence
        end
    end

    # Generate threshold vectors (simplified)
    ensemble_single_Reff_thresholds_vec = [rand(200:300) for _ in 1:nsims]
    ensemble_single_periodsum_vecs = [rand(150:250) for _ in 1:nsims]

    return (
        ensemble_specification = ensemble_specification,
        ensemble_single_incarr = ensemble_single_incarr,
        null_single_incarr = null_single_incarr,
        ensemble_single_Reff_thresholds_vec = ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs = ensemble_single_periodsum_vecs,
    )
end

"""
Performance comparison test with different configurations
"""
function benchmark_multistart_scaling()
    println("\nBenchmarking Multistart Optimization Scaling")
    println("="^60)

    # Test scaling with different numbers of Sobol points
    sobol_points = [10, 25, 50, 100]
    times = Float64[]
    accuracies = Float64[]

    # Create minimal test configuration
    test_spec_vecs = (
        noise_specification_vec = [PoissonNoiseSpecification(1.0)],
        test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
        percent_tested_vec = [0.2],
        ews_metric_specification_vec = [EWSMetricSpecification(Gaussian, Day(7), Day(28), 1)],
        ews_enddate_type_vec = [Reff_start],
        ews_threshold_window_vec = [ExpandingThresholdWindow],
        ews_metric_vec = ["mean"],
    )

    # Generate minimal test data
    ensemble_spec = EnsembleSpecification(
        StateParameters(SVector(990, 10, 0, 0), 0.01),
        DynamicsParameters(0.4, 0.1, 0.1, 0.0),
        TimeParameters(Day(1), 300, Date(2020, 1, 1):Day(1):Date(2020, 10, 27)),
        25,  # Very small for scaling test
    )

    test_data = create_test_data_arrays(ensemble_spec)

    for n_points in sobol_points
        println("\nTesting with $n_points Sobol points...")

        t = @elapsed results = ews_multistart_optimization(
            test_spec_vecs,
            test_data;
            n_sobol_points = n_points,
            maxeval = n_points * 5,
            use_threads = false,
            force = true,
            verbose = false,
            save_results = false,
        )

        push!(times, t)
        push!(accuracies, maximum(results.accuracy))

        println("  Time: $(round(t, digits = 2))s")
        println("  Best accuracy: $(round(accuracies[end], digits = 3))")
    end

    println("\n" * "=" * 60)
    println("SCALING SUMMARY")
    println("=" * 60)
    println("Sobol Points | Time (s) | Best Accuracy")
    println("-" * 40)
    for i in eachindex(sobol_points)
        println("$(lpad(sobol_points[i], 12)) | $(lpad(round(times[i], digits = 1), 8)) | $(round(accuracies[i], digits = 3))")
    end

    return sobol_points, times, accuracies
end

"""
Quick test of multistart optimization functionality
"""
function test_multistart_basic()
    println("Testing Basic Multistart Optimization Functionality")
    println("="^60)

    # Minimal test configuration
    spec_vecs = (
        noise_specification_vec = [PoissonNoiseSpecification(1.0)],
        test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
        percent_tested_vec = [0.2],
        ews_metric_specification_vec = [EWSMetricSpecification(Gaussian, Day(7), Day(28), 1)],
        ews_enddate_type_vec = [Reff_start],
        ews_threshold_window_vec = [ExpandingThresholdWindow],
        ews_metric_vec = ["mean"],
    )

    # Create minimal ensemble
    ensemble_spec = EnsembleSpecification(
        StateParameters(SVector(990, 10, 0, 0), 0.01),
        DynamicsParameters(0.4, 0.1, 0.1, 0.0),
        TimeParameters(Day(1), 200, Date(2020, 1, 1):Day(1):Date(2020, 7, 19)),
        10,  # Very small ensemble
    )

    test_data = create_test_data_arrays(ensemble_spec)

    println("\nRunning multistart optimization...")

    results = ews_multistart_optimization(
        spec_vecs,
        test_data;
        n_sobol_points = 10,
        maxeval = 50,
        use_threads = false,
        force = true,
        verbose = true,
        save_results = false,
    )

    println("\nResults:")
    println("  Scenarios optimized: $(nrow(results))")
    println("  Best accuracy: $(round(maximum(results.accuracy), digits = 4))")
    println("  Best parameters:")

    best_row = results[argmax(results.accuracy), :]
    println("    Percentile: $(round(best_row.ews_threshold_percentile, digits = 3))")
    println("    Consecutive: $(best_row.ews_consecutive_thresholds)")
    println("    Burnin: $(best_row.ews_threshold_burnin)")
    println("    Evaluations: $(best_row.n_evaluations)")
    println("    Convergence: $(best_row.convergence_status)")

    return results
end

# Export functions for easy testing
export benchmark_optimization_methods, benchmark_multistart_scaling, test_multistart_basic

