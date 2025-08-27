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
using StatsBase: StatsBase
using Distributions: Distributions
using ProgressMeter
using FLoops
using Statistics
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
        PoissonNoiseSpecification(1.0),
        PoissonNoiseSpecification(7.0),
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
        noise_specification = noise_specification_vec[1],  # PoissonNoiseSpecification(1.0)
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

"""
Compare accuracy values between grid search and multistart optimization
for all scenarios.
"""
function compare_all_scenarios(grid_results, multistart_results, specification_vecs)
    # Define scenario columns (non-optimized parameters)
    scenario_cols = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
    ]

    # Group grid results by scenario
    grid_grouped = groupby(grid_results, scenario_cols)

    comparison_data = []

    for grid_scenario in grid_grouped
        # Get best grid accuracy for this scenario
        best_grid_idx = argmax(grid_scenario.accuracy)
        best_grid_row = grid_scenario[best_grid_idx, :]
        best_grid_acc = best_grid_row.accuracy

        # Extract scenario identifier
        scenario_id = NamedTuple(
            col => best_grid_row[col] for col in scenario_cols
        )

        # Compare with each multistart configuration
        for (n_sobol, ms_data) in multistart_results
            ms_results = ms_data.res

            # Find matching scenario in multistart results
            ms_scenario = filter(
                row -> all(
                    row[col] == scenario_id[col] for col in scenario_cols
                ), ms_results
            )

            if nrow(ms_scenario) > 0
                best_ms_idx = argmax(ms_scenario.accuracy)
                best_ms_row = ms_scenario[best_ms_idx, :]
                best_ms_acc = best_ms_row.accuracy

                # Calculate difference
                acc_diff = best_ms_acc - best_grid_acc
                relative_diff = abs(acc_diff) / best_grid_acc * 100

                push!(
                    comparison_data, (
                        scenario_id = scenario_id,
                        grid_accuracy = best_grid_acc,
                        grid_percentile = best_grid_row.ews_threshold_percentile,
                        grid_consecutive = best_grid_row.ews_consecutive_thresholds,
                        grid_sensitivity = best_grid_row.sensitivity,
                        grid_specificity = best_grid_row.specificity,
                        ms_accuracy = best_ms_acc,
                        ms_percentile = best_ms_row.ews_threshold_percentile,
                        ms_consecutive = best_ms_row.ews_consecutive_thresholds,
                        ms_sensitivity = best_ms_row.sensitivity,
                        ms_specificity = best_ms_row.specificity,
                        n_sobol = n_sobol,
                        accuracy_diff = acc_diff,
                        relative_diff_pct = relative_diff,
                        matches = isapprox(best_grid_acc, best_ms_acc; atol = 1.0e-4),
                    )
                )
            end
        end
    end

    return DataFrame(comparison_data)
end

"""
Display summary statistics of accuracy comparison.
"""
function display_accuracy_comparison_summary(comparison_df)
    println(styled"\n{bold:Accuracy Comparison Summary}")
    println("-"^40)

    # Group by n_sobol configuration
    for n_sobol in sort(unique(comparison_df.n_sobol))
        config_df = filter(row -> row.n_sobol == n_sobol, comparison_df)

        n_scenarios = nrow(config_df)
        n_matching = sum(config_df.matches)
        match_rate = n_matching / n_scenarios * 100

        mean_diff = mean(config_df.accuracy_diff)
        std_diff = std(config_df.accuracy_diff)
        max_diff = maximum(abs.(config_df.accuracy_diff))
        mean_rel_diff = mean(config_df.relative_diff_pct)

        println(styled"\n{cyan:$n_sobol Sobol points}:")
        println(styled"  Scenarios compared: {yellow:$n_scenarios}")
        println(styled"  Matching accuracy (±1e-4): {green:$n_matching} ({green:$(round(match_rate, digits=1))%})")
        println(styled"  Mean accuracy difference: {yellow:$(round(mean_diff, digits=5))}")
        println(styled"  Std accuracy difference: {yellow:$(round(std_diff, digits=5))}")
        println(styled"  Max absolute difference: {yellow:$(round(max_diff, digits=5))}")
        println(styled"  Mean relative difference: {yellow:$(round(mean_rel_diff, digits=2))%}")

        # Identify scenarios with largest discrepancies
        if max_diff > 1.0e-3
            worst_scenarios = sort(config_df, :accuracy_diff)[1:min(3, nrow(config_df)), :]
            println(styled"  {red:Scenarios with largest negative differences}:")
            for row in eachrow(worst_scenarios)
                metric = row.scenario_id.ews_metric
                noise = get_noise_magnitude_description(row.scenario_id.noise_specification)
                diff = round(row.accuracy_diff, digits = 4)
                println(styled"    - $metric, $noise: {red:$diff}")
            end
        end

        # Show scenarios where multistart performs better
        better_scenarios = filter(row -> row.accuracy_diff > 1.0e-4, config_df)
        if nrow(better_scenarios) > 0
            println(styled"  {green:Scenarios where multistart performs better}: {green:$(nrow(better_scenarios))}")
        end
    end
    return
end

"""
Generate detailed accuracy verification report.
"""
function generate_accuracy_verification_report(
        comparison_df,
        grid_time,
        multistart_results
    )
    println(styled"\n{bold blue:DETAILED VERIFICATION REPORT}")
    println("="^60)

    # Performance summary
    println(styled"\n{bold:Performance Summary}")
    println(styled"Grid search time: {yellow:$(round(grid_time, digits=2))}s")

    for (n_sobol, ms_data) in sort(collect(multistart_results))
        ms_time = ms_data.time
        speedup = round(grid_time / ms_time, digits = 1)
        println(styled"Multistart ($n_sobol points): {yellow:$(round(ms_time, digits=2))}s (speedup: {green:$(speedup)x})")
    end

    # Accuracy verification by metric
    println(styled"\n{bold:Accuracy Verification by EWS Metric}")
    metrics = unique([row.scenario_id.ews_metric for row in eachrow(comparison_df)])

    for metric in metrics
        metric_df = filter(row -> row.scenario_id.ews_metric == metric, comparison_df)

        println(styled"\n{cyan:$metric}:")
        for n_sobol in sort(unique(metric_df.n_sobol))
            sobol_df = filter(row -> row.n_sobol == n_sobol, metric_df)
            n_match = sum(sobol_df.matches)
            n_total = nrow(sobol_df)
            match_pct = round(n_match / n_total * 100, digits = 1)

            println(styled"  $n_sobol Sobol points: {yellow:$n_match/$n_total matching ($match_pct%)}")
        end
    end

    # Parameter convergence analysis
    println(styled"\n{bold:Parameter Convergence Analysis}")

    for n_sobol in sort(unique(comparison_df.n_sobol))
        config_df = filter(row -> row.n_sobol == n_sobol, comparison_df)

        # Calculate parameter differences
        percentile_diffs = abs.(config_df.grid_percentile .- config_df.ms_percentile)
        consecutive_diffs = abs.(config_df.grid_consecutive .- config_df.ms_consecutive)

        mean_perc_diff = mean(percentile_diffs)
        mean_cons_diff = mean(consecutive_diffs)

        println(styled"\n{cyan:$n_sobol Sobol points}:")
        println(styled"  Mean percentile difference: {yellow:$(round(mean_perc_diff, digits=4))}")
        println(styled"  Mean consecutive threshold difference: {yellow:$(round(mean_cons_diff, digits=2))}")
    end

    # Save detailed comparison to CSV for further analysis
    output_dir = outdir("benchmark")
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    output_file = joinpath(
        output_dir,
        "accuracy_verification_$(Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")).csv"
    )

    # Flatten scenario_id for CSV export using NamedTuples
    export_data = map(eachrow(comparison_df)) do row
        # Extract scenario fields
        scenario_fields = Dict{Symbol, Any}()
        for (key, value) in pairs(row.scenario_id)
            if key == :noise_specification
                scenario_fields[:noise_description] = get_noise_magnitude_description(value)
            elseif key == :test_specification
                scenario_fields[:test_description] = get_test_description(value)
            else
                scenario_fields[key] = string(value)
            end
        end

        # Combine with metrics
        return merge(
            NamedTuple(scenario_fields),
            (
                grid_accuracy = row.grid_accuracy,
                grid_percentile = row.grid_percentile,
                grid_consecutive = row.grid_consecutive,
                grid_sensitivity = row.grid_sensitivity,
                grid_specificity = row.grid_specificity,
                ms_accuracy = row.ms_accuracy,
                ms_percentile = row.ms_percentile,
                ms_consecutive = row.ms_consecutive,
                ms_sensitivity = row.ms_sensitivity,
                ms_specificity = row.ms_specificity,
                n_sobol = row.n_sobol,
                accuracy_diff = row.accuracy_diff,
                relative_diff_pct = row.relative_diff_pct,
                matches = row.matches,
            )
        )
    end

    export_df = DataFrame(export_data)
    CSV.write(output_file, export_df)
    return println(styled"\n{green:Detailed comparison saved to: $output_file}")
end

function calculate_scenarios(spec_vecs)
    return length(spec_vecs.noise_specification_vec) *
        length(spec_vecs.test_specification_vec) *
        length(spec_vecs.percent_tested_vec) *
        length(spec_vecs.ews_metric_specification_vec) *
        length(spec_vecs.ews_enddate_type_vec) *
        length(spec_vecs.ews_threshold_window_vec) *
        length(spec_vecs.ews_metric_vec)
end

function calculate_grid_points(spec_vecs)
    return length(spec_vecs.ews_threshold_percentile_vec) *
        length(spec_vecs.ews_consecutive_thresholds_vec) *
        length(spec_vecs.ews_threshold_burnin_vec)
end

function generate_ensemble_data(ensemble_specification, null_specification, ensemble_outbreak_specification)
    println(styled"  Generating ensemble dynamics...")

    # Generate ensemble data by calling the core simulation function directly
    ensemble_data = generate_single_ensemble(ensemble_specification; seed = 42)

    println(styled"  Generating null dynamics...")

    # Generate null data
    null_data = generate_single_ensemble(null_specification; seed = 42)

    println(styled"  Processing outbreak data...")

    # Generate incidence arrays for the outbreak specification using the incidence vectors
    ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
        ensemble_data[:ensemble_inc_vecs], ensemble_outbreak_specification
    )

    null_inc_arr, null_thresholds_vec = create_inc_infec_arr(
        null_data[:ensemble_inc_vecs], ensemble_outbreak_specification
    )

    return (
        ensemble_specification = ensemble_specification,
        ensemble_single_incarr = ensemble_inc_arr,
        null_single_incarr = null_inc_arr,
        ensemble_single_Reff_thresholds_vec = ensemble_data[:ensemble_Reff_thresholds_vec],
        ensemble_single_periodsum_vecs = ensemble_thresholds_vec,
    )
end

function generate_single_ensemble(ensemble_spec; seed)
    @unpack state_parameters, dynamics_parameter_specification, time_parameters, nsims = ensemble_spec
    @unpack tstep, tlength, trange = time_parameters

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states), 2}(
        undef, tlength, nsims
    )

    ensemble_inc_vecs = Array{typeof(SVector(0)), 2}(
        undef, tlength, nsims
    )

    ensemble_beta_arr = zeros(Float64, tlength)
    ensemble_Reff_arr = zeros(Float64, tlength, nsims)
    ensemble_Reff_thresholds_vec = Vector{Array{Int64, 2}}(
        undef, size(ensemble_inc_vecs, 2)
    )

    dynamics_parameters = Vector{DynamicsParameters}(undef, nsims)

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        dynamics_parameters[sim] = DynamicsParameters(
            dynamics_parameter_specification; seed = run_seed
        )

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            dynamics_parameters[sim],
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_seir_arr = convert_svec_to_array(ensemble_seir_vecs)

    for sim in axes(ensemble_inc_vecs, 2)
        calculateReffective_t!(
            @view(ensemble_Reff_arr[:, sim]),
            ensemble_beta_arr,
            dynamics_parameters[sim],
            1,
            @view(ensemble_seir_arr[:, :, sim]),
        )

        ensemble_Reff_thresholds_vec[sim] = Reff_ge_than_one(
            @view(ensemble_Reff_arr[:, sim])
        )
    end

    return (
        ensemble_seir_arr = ensemble_seir_arr,
        ensemble_spec = ensemble_spec,
        dynamics_parameters = dynamics_parameters,
        ensemble_Reff_arr = ensemble_Reff_arr,
        ensemble_Reff_thresholds_vec = ensemble_Reff_thresholds_vec,
        ensemble_inc_vecs = ensemble_inc_vecs,
    )
end

# Run benchmark if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
