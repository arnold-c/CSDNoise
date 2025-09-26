# Common benchmarking functions that can be reused across scripts

using CSV: CSV

export generate_ensemble_data, generate_single_ensemble,
    calculate_scenarios, calculate_grid_points,
    compare_optimization_results, display_benchmark_summary,
    compare_all_scenarios, display_accuracy_comparison_summary,
    generate_accuracy_verification_report, create_ensemble_specs,
    save_benchmark_comparison_results, EnsembleSpecsParameters

"""
    generate_ensemble_data(ensemble_specification, null_specification, ensemble_outbreak_specification)

Generate ensemble and null data arrays for benchmarking optimization methods.
"""
function generate_ensemble_data(ensemble_specification, null_specification, ensemble_outbreak_specification)
    println(styled"  Generating ensemble dynamics...")

    # Generate ensemble data by calling the core simulation function directly
    ensemble_data = generate_single_ensemble(ensemble_specification; seed = 42)

    println(styled"  Generating null dynamics...")

    # Generate null data
    null_data = generate_single_ensemble(null_specification; seed = 42)

    println(styled"  Processing outbreak data...")

    # Generate incidence arrays for the outbreak specification using the incidence vectors
    emergent_incidence_arr, emergent_outbreak_threshold_vecs = create_inc_infec_arr(
        ensemble_data[:ensemble_inc_vecs], ensemble_outbreak_specification
    )

    null_incidence_arr, null_outbreak_threshold_vecs = create_inc_infec_arr(
        null_data[:ensemble_inc_vecs], ensemble_outbreak_specification
    )

    return (
        emergent_incidence_arr,
        null_incidence_arr,
        ensemble_single_Reff_thresholds_vec = ensemble_data[:ensemble_Reff_thresholds_vec],
        emergent_outbreak_threshold_vecs = emergent_outbreak_threshold_vecs,
    )
end

"""
    generate_single_ensemble(ensemble_spec; seed)

Generate a single ensemble simulation with the given specification and seed.
"""
function generate_single_ensemble(ensemble_spec::EnsembleSpecification; seed::Int64)
    @unpack state_parameters, dynamics_parameter_specification, time_parameters, nsims = ensemble_spec
    @unpack tstep, tlength, trange = time_parameters

    init_states_sv = SVector(state_parameters.init_states)

    # Get concrete type to avoid abstract element types
    init_state_type = typeof(init_states_sv)
    ensemble_seir_vecs = Array{init_state_type, 2}(
        undef, tlength, nsims
    )

    # Use concrete SVector type instead of typeof(SVector(0))
    ensemble_inc_vecs = Array{Int64, 2}(
        undef, tlength, nsims
    )

    beta_vec = zeros(Float64, tlength)
    ensemble_Reff_arr = zeros(Float64, tlength, nsims)
    ensemble_Reff_thresholds_vec = Vector{Array{Int64, 2}}(
        undef, size(ensemble_inc_vecs, 2)
    )

    dynamics_parameters = Vector{DynamicsParameters}(undef, nsims)

    # Use explicit loop instead of broadcasting to avoid runtime dispatch
    for i in eachindex(beta_vec)
        beta_vec[i] = calculate_beta_amp(
            dynamics_parameter_specification.beta_mean,
            dynamics_parameter_specification.beta_force,
            trange[i];
            seasonality = dynamics_parameter_specification.seasonality
        )
    end

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        local dynp = DynamicsParameters(
            dynamics_parameter_specification; seed = run_seed
        )
        dynamics_parameters[sim] = dynp

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            beta_vec,
            init_states_sv,
            dynp,
            time_parameters,
            run_seed,
        )

        calculateReffective_t!(
            @view(ensemble_Reff_arr[:, sim]),
            beta_vec,
            dynamics_parameters[sim],
            1,
            @view(ensemble_seir_vecs[:, sim])
        )

        ensemble_Reff_thresholds_vec[sim] = Reff_ge_than_one(
            @view(ensemble_Reff_arr[:, sim])
        )
    end

    return (
        ensemble_seir_vecs = ensemble_seir_vecs,
        ensemble_spec = ensemble_spec,
        dynamics_parameters = dynamics_parameters,
        ensemble_Reff_arr = ensemble_Reff_arr,
        ensemble_Reff_thresholds_vec = ensemble_Reff_thresholds_vec,
        ensemble_inc_vecs = ensemble_inc_vecs,
    )
end

"""
    calculate_scenarios(spec_vecs)

Calculate the total number of scenarios from specification vectors.
"""
function calculate_scenarios(spec_vecs)
    return length(spec_vecs.noise_specification_vec) *
        length(spec_vecs.test_specification_vec) *
        length(spec_vecs.percent_tested_vec) *
        length(spec_vecs.ews_metric_specification_vec) *
        length(spec_vecs.ews_enddate_type_vec) *
        length(spec_vecs.ews_threshold_window_vec) *
        length(spec_vecs.ews_metric_vec) *
        length(spec_vecs.ews_threshold_burnin_vec)
end

"""
    calculate_grid_points(spec_vecs)

Calculate the total number of grid points from specification vectors.
"""
function calculate_grid_points(spec_vecs)
    return length(spec_vecs.ews_threshold_quantile_vec) *
        length(spec_vecs.ews_consecutive_thresholds_vec)
end

"""
    compare_optimization_results(results1, results2, scenario_cols)

Compare two sets of optimization results for matching scenarios.
"""
function compare_optimization_results(results1, results2, scenario_cols)
    comparison_data = []

    # Group results by scenario
    grouped1 = groupby(results1, scenario_cols)
    grouped2 = groupby(results2, scenario_cols)

    for group1 in grouped1
        # Find matching scenario in results2
        scenario_key = NamedTuple(col => group1[1, col] for col in scenario_cols)

        matching_group2 = filter(
            group2 -> all(
                group2[1, col] == scenario_key[col] for col in scenario_cols
            ), grouped2
        )

        if length(matching_group2) > 0
            group2 = matching_group2[1]

            best_acc1 = maximum(group1.accuracy)
            best_acc2 = maximum(group2.accuracy)

            push!(
                comparison_data, (
                    scenario = scenario_key,
                    accuracy1 = best_acc1,
                    accuracy2 = best_acc2,
                    accuracy_diff = best_acc2 - best_acc1,
                    relative_diff = (best_acc2 - best_acc1) / best_acc1 * 100,
                )
            )
        end
    end

    return DataFrame(comparison_data)
end

"""
    display_benchmark_summary(results, label)

Display summary statistics for benchmark results.
"""
function display_benchmark_summary(results, label)
    println(styled"\n{bold:$label Summary}")
    println("-"^40)

    best_accuracy = maximum(results.accuracy)
    mean_accuracy = mean(results.accuracy)
    std_accuracy = std(results.accuracy)

    println(styled"Best accuracy: {green:$(round(best_accuracy, digits=4))}")
    println(styled"Mean accuracy: {yellow:$(round(mean_accuracy, digits=4))}")
    println(styled"Std accuracy: {yellow:$(round(std_accuracy, digits=4))}")

    # Find best parameters
    best_idx = argmax(results.accuracy)
    best_row = results[best_idx, :]

    println(styled"Best parameters:")
    println(styled"  Threshold quantile: {cyan:$(round(best_row.ews_threshold_quantile, digits=3))}")
    println(styled"  Consecutive thresholds: {cyan:$(best_row.ews_consecutive_thresholds)}")
    println(styled"  Sensitivity: {green:$(round(best_row.sensitivity, digits=4))}")
    return println(styled"  Specificity: {green:$(round(best_row.specificity, digits=4))}")
end

"""
    compare_all_scenarios(grid_results, multistart_results, specification_vecs)

Compare accuracy values between grid search and multistart optimization for all scenarios.
More sophisticated than compare_optimization_results - handles multiple multistart configurations.
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
                    # Don't match on optimization parameters
                    row[col] == scenario_id[col] for col in setdiff(scenario_cols, [:ews_threshold_window, :ews_threshold_burnin])
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
                        grid_quantile = best_grid_row.ews_threshold_quantile,
                        grid_consecutive = best_grid_row.ews_consecutive_thresholds,
                        grid_sensitivity = best_grid_row.sensitivity,
                        grid_specificity = best_grid_row.specificity,
                        ms_accuracy = best_ms_acc,
                        ms_quantile = best_ms_row.threshold_quantile,
                        ms_consecutive = best_ms_row.consecutive_thresholds,
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
    display_accuracy_comparison_summary(comparison_df)

Display summary statistics of accuracy comparison between grid search and multistart.
"""
function display_accuracy_comparison_summary(comparison_df)
    log_both(styled"\n{bold:Accuracy Comparison Summary}")
    log_both("-"^40)

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

        log_both(styled"\n{cyan:$n_sobol Sobol points}:")
        log_both(styled"  Scenarios compared: {yellow:$n_scenarios}")
        log_both(styled"  Matching accuracy (±1e-4): {green:$n_matching} ({green:$(round(match_rate, digits=1))%})")
        log_both(styled"  Mean accuracy difference: {yellow:$(round(mean_diff, digits=5))}")
        log_both(styled"  Std accuracy difference: {yellow:$(round(std_diff, digits=5))}")
        log_both(styled"  Max absolute difference: {yellow:$(round(max_diff, digits=5))}")
        log_both(styled"  Mean relative difference: {yellow:$(round(mean_rel_diff, digits=2))%}")

        # Identify scenarios with largest discrepancies
        if max_diff > 1.0e-3
            worst_scenarios = sort(config_df, :accuracy_diff)[1:min(3, nrow(config_df)), :]
            log_both(styled"  {red:Scenarios with largest negative differences}:")
            for row in eachrow(worst_scenarios)
                metric = row.scenario_id.ews_metric
                noise = get_noise_magnitude_description(row.scenario_id.noise_specification)
                diff = round(row.accuracy_diff, digits = 4)
                log_both(styled"    - $metric, $noise: {red:$diff}")
            end
        end

        # Show scenarios where multistart performs better
        better_scenarios = filter(row -> row.accuracy_diff > 1.0e-4, config_df)
        if nrow(better_scenarios) > 0
            log_both(styled"  {green:Scenarios where multistart performs better}: {green:$(nrow(better_scenarios))}")
        end
    end
    return
end

"""
    generate_accuracy_verification_report(comparison_df, grid_time, multistart_results)

Generate detailed accuracy verification report with CSV export.
"""
function generate_accuracy_verification_report(
        comparison_df,
        grid_time,
        multistart_results
    )
    log_both(styled"\n{bold blue:DETAILED VERIFICATION REPORT}")
    log_both("="^60)

    # Performance summary
    log_both(styled"\n{bold:Performance Summary}")
    log_both(styled"Grid search time: {yellow:$(round(grid_time, digits=2))}s")

    for (n_sobol, ms_data) in sort(collect(multistart_results))
        ms_time = ms_data.time
        speedup = round(grid_time / ms_time, digits = 1)
        log_both(styled"Multistart ($n_sobol points): {yellow:$(round(ms_time, digits=2))}s (speedup: {green:$(speedup)x})")
    end

    # Accuracy verification by metric
    log_both(styled"\n{bold:Accuracy Verification by EWS Metric}")
    metrics = unique([row.scenario_id.ews_metric for row in eachrow(comparison_df)])

    for metric in metrics
        metric_df = filter(row -> row.scenario_id.ews_metric == metric, comparison_df)

        log_both(styled"\n{cyan:$metric}:")
        for n_sobol in sort(unique(metric_df.n_sobol))
            sobol_df = filter(row -> row.n_sobol == n_sobol, metric_df)
            n_match = sum(sobol_df.matches)
            n_total = nrow(sobol_df)
            match_pct = round(n_match / n_total * 100, digits = 1)

            log_both(styled"  $n_sobol Sobol points: {yellow:$n_match/$n_total matching ($match_pct%)}")
        end
    end

    # Parameter convergence analysis
    log_both(styled"\n{bold:Parameter Convergence Analysis}")

    for n_sobol in sort(unique(comparison_df.n_sobol))
        config_df = filter(row -> row.n_sobol == n_sobol, comparison_df)

        # Calculate parameter differences
        quantile_diffs = abs.(config_df.grid_quantile .- config_df.ms_quantile)
        consecutive_diffs = abs.(config_df.grid_consecutive .- config_df.ms_consecutive)

        mean_perc_diff = mean(quantile_diffs)
        mean_cons_diff = mean(consecutive_diffs)

        log_both(styled"\n{cyan:$n_sobol Sobol points}:")
        log_both(styled"  Mean quantile difference: {yellow:$(round(mean_perc_diff, digits=4))}")
        log_both(styled"  Mean consecutive threshold difference: {yellow:$(round(mean_cons_diff, digits=2))}")
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
                grid_quantile = row.grid_quantile,
                grid_consecutive = row.grid_consecutive,
                grid_sensitivity = row.grid_sensitivity,
                grid_specificity = row.grid_specificity,
                ms_accuracy = row.ms_accuracy,
                ms_quantile = row.ms_quantile,
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
    return log_both(styled"\n{green:Detailed comparison saved to: $output_file}")
end

"""
    create_ensemble_specs(nsims)

Create ensemble specifications for a given number of simulations.
Returns ensemble_specification, null_specification, and outbreak_specification.
"""
function create_ensemble_specs(params::EnsembleSpecsParameters)
    # Base configuration
    ensemble_time_specification = SimTimeParameters(;
        burnin = 365.0 * params.burnin_years,
        tmin = 0.0,
        tmax = 365.0 * params.nyears,
        tstep = 1.0
    )

    # Calculate dynamics parameters
    mu = calculate_mu(params.annual_births_per_k)
    beta_mean = calculate_beta(params.R_0, params.gamma, mu, 1, params.ensemble_state_specification.init_states.N)
    epsilon = calculate_import_rate(mu, params.R_0, params.ensemble_state_specification.init_states.N)

    min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
        params.target_Reff,
        params.target_years,
        params.ensemble_state_specification.init_states.S,
        params.ensemble_state_specification.init_states.N,
        params.R_0,
        mu
    )

    @assert params.max_vaccination_coverage < min_burnin_vaccination_coverage "Set the maximum vaccination coverage for the emergent time series to $(params.max_vaccination_coverage), but it should be less than $min_burnin_vaccination_coverage"

    ensemble_dynamics_specification = DynamicsParameterSpecification(
        beta_mean,
        0.0,
        SeasonalityFunction(CosineSeasonality()),
        params.sigma,
        params.gamma,
        mu,
        27.0,
        epsilon,
        params.R_0,
        min_burnin_vaccination_coverage,
        1.0,
        params.min_vaccination_coverage,
        params.max_vaccination_coverage
    )

    null_dynamics_specification = DynamicsParameterSpecification(
        map(
            pn -> getproperty(ensemble_dynamics_specification, pn),
            filter(
                name -> name != :min_vaccination_coverage && name != :max_vaccination_coverage,
                propertynames(ensemble_dynamics_specification)
            )
        )...,
        nothing, nothing
    )

    # Create specifications with specified number of simulations
    ensemble_specification = EnsembleSpecification(
        params.ensemble_state_specification,
        ensemble_dynamics_specification,
        ensemble_time_specification,
        params.nsims
    )

    null_specification = EnsembleSpecification(
        params.ensemble_state_specification,
        null_dynamics_specification,
        ensemble_time_specification,
        params.nsims
    )

    outbreak_specification = OutbreakSpecification(5, 30, 500)

    return ensemble_specification, null_specification, outbreak_specification
end

"""
    save_benchmark_comparison_results(results_dict, filename_prefix; output_dir=outdir("benchmark"))

Save benchmark comparison results to CSV file. Generalized version that can handle
different types of benchmark comparisons.
"""
function save_benchmark_comparison_results(results_dict, filename_prefix; output_dir = outdir("benchmark"), full_filename = "")
    mkpath(output_dir)

    if isempty(full_filename)
        # Use default timestamped naming
        timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
        filename = "$(filename_prefix)_$(timestamp).csv"
    else
        # Use provided full filename
        filename = "$(full_filename).csv"
    end

    # Create summary DataFrame from results dictionary
    summary_data = []
    for (key, result) in results_dict
        # Handle different result structures
        if haskey(result, :best_accuracy) && haskey(result, :best_params)
            # Ensemble size comparison format
            push!(
                summary_data, (
                    comparison_key = key,
                    data_time = get(result, :data_time, missing),
                    opt_time = get(result, :opt_time, missing),
                    total_time = get(result, :total_time, missing),
                    best_accuracy = result.best_accuracy,
                    best_quantile = result.best_params.threshold_quantile,
                    best_consecutive = result.best_params.consecutive_thresholds,
                    sensitivity = result.best_params.sensitivity,
                    specificity = result.best_params.specificity,
                )
            )
        else
            # Generic format - just store what's available
            row_data = (comparison_key = key,)
            for (field, value) in pairs(result)
                if field != :results_df  # Skip large DataFrames
                    row_data = merge(row_data, NamedTuple{(field,)}((value,)))
                end
            end
            push!(summary_data, row_data)
        end
    end

    summary_df = DataFrame(summary_data)
    CSV.write(joinpath(output_dir, filename), summary_df)
    return println(styled"\n{green:Results saved to: $filename}")
end
