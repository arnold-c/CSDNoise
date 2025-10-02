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
using Dates
using Random
using StyledStrings
using Try
using FLoops
using DataFrames
using Logging
using LoggingExtras
using StructArrays: StructVector
using StatsBase
using DataFrames

function main()
    # Create consistent filename base for both CSV and markdown files
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    filename_base = "optimization-speed-comparison_$(timestamp)"

    # Setup dual logging in benchmark/optimization-speed subdirectory
    benchmark_dir = outdir("benchmark", "optimization-speed")
    markdown_filename = setup_dual_logging("benchmark-optimization-speed"; title = "Optimization Speed & Accuracy Benchmark Report", output_dir = benchmark_dir, filename_base = filename_base)

    try
        log_both(styled"{bold blue:Optimization Speed & Accuracy Benchmark}")
        log_both("="^60)

        # Set random seed for reproducibility
        Random.seed!(42)

        # Setup ensemble configuration (using realistic parameters from ensemble-sim_ews-optimization.jl)
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
            SeasonalityFunction(CosineSeasonality()),
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
            ensemble_state_specification,
            ensemble_dynamics_specification,
            ensemble_time_specification,
            ensemble_nsims,
        )

        ensemble_specification_vec = [ensemble_specification]

        null_specification = EnsembleSpecification(
            ensemble_state_specification,
            null_dynamics_specification,
            ensemble_time_specification,
            ensemble_nsims,
        )

        null_specification_vec = [null_specification]

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

        ews_method_vec = [EWSMethod(Backward())]
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

        ews_enddate_type_vec = [EWSEndDateType(ReffStart())]
        ews_threshold_window_vec = [EWSThresholdWindowType(ExpandingThresholdWindow())]
        ews_threshold_quantile_vec = collect(0.5:0.01:0.99)
        ews_consecutive_thresholds_vec = collect(2:1:20)
        ews_threshold_burnin_vec = [Year(5)]

        specification_vecs = (;
            ensemble_specification_vec,
            null_specification_vec,
            noise_specification_vec,
            test_specification_vec,
            percent_tested_vec,
            ews_metric_specification_vec,
            ews_enddate_type_vec,
            ews_threshold_window_vec,
            ews_threshold_burnin_vec,
            ews_threshold_quantile_vec,
            ews_consecutive_thresholds_vec,
            ews_metric_vec,
        )
        # Remove unused vectors for multivariate optimization
        optimization_spec_vecs = Base.structdiff(specification_vecs, (; ews_threshold_quantile_vec, ews_consecutive_thresholds_vec))

        log_both(styled"Generating ensemble data...")
        data_arrs = generate_ensemble_data(ensemble_specification, null_specification, ensemble_outbreak_specification)

        # Calculate search space size
        n_scenarios = calculate_scenarios(specification_vecs)
        n_grid_points = calculate_grid_points(specification_vecs)
        total_grid_evaluations = n_scenarios * n_grid_points

        log_both(styled"Search space:")
        log_both(styled"  Scenarios: {cyan:$n_scenarios}")
        log_both(styled"  Grid points per scenario: {cyan:$n_grid_points}")
        log_both(styled"  Total grid evaluations: {yellow:$total_grid_evaluations}")

        # Compilation run - single scenario to pre-compile functions
        log_both(styled"\n{bold yellow:COMPILATION RUN}")
        log_both("-"^40)
        log_both(styled"Running single scenario to compile functions...")

        compilation_spec_vecs = map(spec -> [spec[1]], specification_vecs)
        compilation_optim_spec_vecs = map(spec -> [spec[1]], optimization_spec_vecs)

        log_both(styled"\n{bold yellow:GRID SEARCH (DATAFRAME) COMPILATION RUN}")
        compilation_time = @elapsed _ = ews_hyperparam_optimization(
            compilation_spec_vecs,
            data_arrs;
            disable_time_check = true,
            force = true,
            return_df = true,
            save_results = false,
            verbose = false,
        )
        log_both(styled"Compilation completed in {yellow:$(round(compilation_time, digits=2))} seconds")

        log_both(styled"\n{bold yellow:GRID SEARCH (STRUCTVECTOR) COMPILATION RUN}")
        structvector_compilation_time = @elapsed structvector_grid_results = ews_hyperparam_gridsearch_structvector(
            compilation_spec_vecs,
            data_arrs;
            executor = FLoops.ThreadedEx(),
            disable_time_check = true,
            force = true,
            return_results = true,
            save_results = false,
            save_checkpoints = false,
            verbose = false,
        )
        log_both(styled"Compilation completed in {yellow:$(round(structvector_compilation_time, digits=2))} seconds")

        log_both(styled"\n{bold yellow:MULTIVARIATE OPTIMIZATION COMPILATION RUN}")
        # Compilation run for multistart optimization
        multistart_compilation_time = @elapsed _ = ews_multistart_optimization(
            compilation_optim_spec_vecs,
            data_arrs;
            n_sobol_points = 10,
            maxeval = 10,
            executor = FLoops.ThreadedEx(),
            force = true,
            return_df = true,
            save_results = false,
            verbose = false,
            disable_time_check = true,
        )
        log_both(styled"Multivariate optimization compilation completed in {yellow:$(round(multistart_compilation_time, digits=2))} seconds")

        # Benchmark grid search
        log_both(styled"\n{bold green:GRID SEARCH BENCHMARK}")
        log_both("-"^40)

        grid_time = @elapsed grid_results = ews_hyperparam_optimization(
            specification_vecs,
            data_arrs,
            disable_time_check = true,
            force = true,
            return_df = true,
            save_results = false,
            verbose = false,
        )

        log_both("Grid search complete ($grid_time)")

        # Compilation run - single scenario to pre-compile functions
        # Benchmark StructVector grid search
        log_both(styled"\n{bold green:STRUCTVECTOR GRID SEARCH BENCHMARK}")
        log_both("-"^40)

        structvector_grid_time = @elapsed structvector_grid_results = ews_hyperparam_gridsearch_structvector(
            specification_vecs,
            data_arrs;
            # executor = FLoops.ThreadedEx(),
            executor = FLoops.SequentialEx(),
            disable_time_check = true,
            force = true,
            return_results = true,
            save_results = false,
            save_checkpoints = false,
            verbose = false,
        )

        log_both("StructVector grid search complete ($structvector_grid_time)")


        # Benchmark multistart with different configurations
        multistart_configs = [
            (n_sobol = 20, name = "Fast"),
            (n_sobol = 50, name = "Balanced"),
            (n_sobol = 100, name = "Thorough"),
            (n_sobol = 150, name = "Very Thorough"),
            # (n_sobol = 250, name = "Extremely Thorough"),
        ]

        log_both(styled"\n{bold green:MULTISTART OPTIMIZATION BENCHMARK}")
        log_both("-"^40)

        multistart_results = Dict()

        for config in multistart_configs
            n_sobol = config.n_sobol
            name = config.name

            log_both(styled"\n{bold:$name} ({cyan:$n_sobol} Sobol points):")

            ms_time = @elapsed results = ews_multistart_optimization(
                optimization_spec_vecs,
                data_arrs;
                quantile_bounds = (0.5, 0.99),  # Custom bounds
                consecutive_bounds = (2.0, 30.0),  # Custom bounds
                n_sobol_points = n_sobol,
                maxeval = 1000,
                # executor = FLoops.SequentialEx(),          # Executor for sequential processing
                executor = FLoops.ThreadedEx(),          # Executor for parallel processing
                force = true,
                return_df = true,
                save_results = false,
                verbose = false,
                disable_time_check = true
            )

            multistart_results[n_sobol] = (; time = ms_time, res = results)

            log_both("Multistart $name complete ($ms_time)")
        end

        # Summary comparison - fix to specific scenario for fair comparison
        log_both(styled"\n{bold blue:SUMMARY - FIXED SCENARIO COMPARISON}")
        log_both("="^60)

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

        log_both(styled"Reference scenario:")
        log_both(styled"  Noise: {cyan:$(get_noise_magnitude_description(reference_scenario.noise_specification))}")
        log_both(styled"  Test spec: {cyan:$(get_test_description(reference_scenario.test_specification))}")
        log_both(styled"  EWS metric: {cyan:$(reference_scenario.ews_metric)}")
        log_both("")

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

        log_both(styled"Grid Search (DataFrame, reference scenario):")
        log_both(styled"  Time: {yellow:$(round(grid_time, digits=2))}s")
        log_both(styled"  Best accuracy: {green:$(round(best_grid_acc, digits=4))}")
        log_both(styled"  Best parameters (optimized):")
        log_both(styled"    Threshold quantile: {cyan:$(round(best_grid_row.ews_threshold_quantile, digits=3))}")
        log_both(styled"    Consecutive thresholds: {cyan:$(best_grid_row.ews_consecutive_thresholds)}")
        log_both(styled"    Sensitivity: {green:$(round(best_grid_row.sensitivity, digits=4))}")
        log_both(styled"    Specificity: {green:$(round(best_grid_row.specificity, digits=4))}")

        # Filter StructVector grid results to reference scenario
        structvector_grid_scenario_results = filter(structvector_grid_results) do row
            reference_scenario_burnin_days = Dates.Day(round(Dates.days(reference_scenario.ews_threshold_burnin)))

            row.noise_specification == reference_scenario.noise_specification &&
                row.test_specification == reference_scenario.test_specification &&
                row.percent_tested == reference_scenario.percent_tested &&
                row.ews_metric == reference_scenario.ews_metric &&
                row.ews_metric_specification == reference_scenario.ews_metric_specification &&
                row.ews_enddate_type == reference_scenario.ews_enddate_type &&
                row.ews_threshold_window == reference_scenario.ews_threshold_window &&
                row.ews_threshold_burnin == reference_scenario_burnin_days
        end

        best_structvector_grid_acc = maximum(structvector_grid_scenario_results.accuracy)
        best_structvector_grid_idx = argmax(structvector_grid_scenario_results.accuracy)
        best_structvector_grid_row = structvector_grid_scenario_results[best_structvector_grid_idx]

        structvector_speedup = round(grid_time / structvector_grid_time, digits = 1)
        structvector_reduction_pct = round(100 * (1 - structvector_grid_time / grid_time), digits = 1)

        log_both(styled"\nGrid Search (StructVector + FLoops, reference scenario):")
        log_both(styled"  Time: {yellow:$(round(structvector_grid_time, digits=2))}s")
        log_both(styled"  Speedup vs DataFrame: {green:$(structvector_speedup)}x")
        log_both(styled"  Time reduction: {green:$(structvector_reduction_pct)}%")
        log_both(styled"  Best accuracy: {green:$(round(best_structvector_grid_acc, digits=4))}")
        log_both(styled"  Δ Accuracy vs DataFrame Grid: {yellow:$(round(best_structvector_grid_acc - best_grid_acc, digits=4))}")
        log_both(styled"  Best parameters (optimized):")
        log_both(styled"    Threshold quantile: {cyan:$(round(best_structvector_grid_row.threshold_quantile, digits=3))}")
        log_both(styled"    Consecutive thresholds: {cyan:$(best_structvector_grid_row.consecutive_thresholds)}")
        log_both(styled"    Sensitivity: {green:$(round(best_structvector_grid_row.sensitivity, digits=4))}")
        log_both(styled"    Specificity: {green:$(round(best_structvector_grid_row.specificity, digits=4))}")

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
                    row.ews_enddate_type == reference_scenario.ews_enddate_type,
                results
            )

            if nrow(ms_scenario_results) == 0
                log_both(styled"\nMultistart ($name - $n_sobol Sobol points): {red:No results for reference scenario}")
                continue
            end

            best_ms_acc = maximum(ms_scenario_results.accuracy)
            acc_diff = best_ms_acc - best_grid_acc

            speedup = round(grid_time / ms_time, digits = 1)
            reduction_pct = round(100 * (1 - ms_time / grid_time), digits = 1)

            # Find best accuracy row in filtered results
            best_idx = argmax(ms_scenario_results.accuracy)
            best_row = ms_scenario_results[best_idx, :]

            log_both(styled"\nMultistart ($name - $n_sobol Sobol points, reference scenario):")
            log_both(styled"  Time: {yellow:$(round(ms_time, digits=2))}s")
            log_both(styled"  Speedup: {green:$(speedup)}x")
            log_both(styled"  Time reduction: {green:$(reduction_pct)}%")
            log_both(styled"  Best accuracy: {green:$(round(best_ms_acc, digits=4))}")
            log_both(styled"  Δ Accuracy vs Grid: {yellow:$(round(acc_diff, digits=4))}")

            # Print optimized parameters (only the ones that vary)
            log_both(styled"  Best parameters (optimized):")
            log_both(styled"    Threshold quantile: {cyan:$(round(best_row.threshold_quantile, digits=3))}")
            log_both(styled"    Consecutive thresholds: {cyan:$(best_row.consecutive_thresholds)}")
            log_both(styled"    Sensitivity: {green:$(round(best_row.sensitivity, digits=4))}")
            log_both(styled"    Specificity: {green:$(round(best_row.specificity, digits=4))}")
        end

        # NEW: Basic accuracy comparison summary
        log_both(styled"\n{bold blue:ACCURACY COMPARISON SUMMARY}")
        log_both("="^60)

        log_both(styled"DataFrame Grid Search:")
        log_both(styled"  Total scenarios: {cyan:$(nrow(grid_results))}")
        log_both(styled"  Best accuracy: {green:$(round(maximum(grid_results.accuracy), digits=4))}")
        log_both(styled"  Mean accuracy: {yellow:$(round(mean(grid_results.accuracy), digits=4))}")

        log_both(styled"\nStructVector Grid Search:")
        log_both(styled"  Total scenarios: {cyan:$(length(structvector_grid_results))}")
        log_both(styled"  Best accuracy: {green:$(round(maximum(structvector_grid_results.accuracy), digits=4))}")
        log_both(styled"  Mean accuracy: {yellow:$(round(mean(structvector_grid_results.accuracy), digits=4))}")
        log_both(styled"  Accuracy difference vs DataFrame: {yellow:$(round(mean(structvector_grid_results.accuracy) - mean(grid_results.accuracy), digits=6))}")

        for config in multistart_configs
            n_sobol = config.n_sobol
            name = config.name
            results = multistart_results[n_sobol][:res]

            log_both(styled"\nMultistart ($name - $n_sobol Sobol points):")
            log_both(styled"  Total scenarios: {cyan:$(nrow(results))}")
            log_both(styled"  Best accuracy: {green:$(round(maximum(results.accuracy), digits=4))}")
            log_both(styled"  Mean accuracy: {yellow:$(round(mean(results.accuracy), digits=4))}")
        end

        @info ""
        @info "---"
        @info ""
        @info "**Report saved to:** `$markdown_filename`"

        return grid_results, structvector_grid_results, multistart_results
    finally
        cleanup_logging()
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
