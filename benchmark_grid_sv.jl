#%%
using JET
using BenchmarkTools
using DrWatson
@quickactivate("CSDNoise")
using CSDNoise
using Dates
using NLopt
using StructArrays
using FLoops
using DataFrames
using StatsBase: mean


# Helper function to create test specification vectors for grid search
function create_gridsearch_test_specification_vectors()
    ensemble_spec, null_spec, _ = create_ensemble_specs(3)
    ensemble_specification_vec = [ensemble_spec]
    null_specification_vec = [null_spec]

    noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))]
    test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)]
    percent_tested_vec = [1.0]

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
    ews_threshold_quantile_vec = collect(0.5:0.02:0.99)
    ews_consecutive_thresholds_vec = collect(2:2:30)
    ews_threshold_burnin_vec = [Year(5)]

    return (;
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
end

#%%
specification_vecs = create_gridsearch_test_specification_vectors()

multivariate_optimization_vecs = specification_vecs[
    setdiff(
        propertynames(specification_vecs), (
            :ews_consecutive_thresholds_vec,
            :ews_threshold_quantile_vec,
        )
    ),
]

univariate_optimization_vecs = (
    specification_vecs[setdiff(propertynames(specification_vecs), (:ews_threshold_quantile_vec,))]...,
    ews_threshold_quantile_vec = [0.9],
)

all_scenarios = create_gridsearch_scenarios_structvector(specification_vecs)
univariate_optim_scenarios = create_gridsearch_scenarios_structvector(univariate_optimization_vecs)

existing_results = StructVector(OptimizationResult[])

missing_scenarios = find_missing_scenarios(all_scenarios, existing_results)
univariate_optim_missing_scenarios = find_missing_scenarios(univariate_optim_scenarios, existing_results)

@assert length(univariate_optim_missing_scenarios) == length(missing_scenarios) / length(specification_vecs.ews_threshold_quantile_vec)

ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(3)
data_arrs = generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)

#%%
# @code_warntype evaluate_gridsearch_scenarios(
#     missing_scenarios,
#     data_arrs;
#     executor = FLoops.SequentialEx(),
#     # executor = FLoops.ThreadedEx(),
#     save_results = false,
#     save_checkpoints = false,
#     checkpoint_dir = "",
#     verbose = false
# )
#
#%%
# @profview evaluate_gridsearch_scenarios(
#     missing_scenarios,
#     data_arrs;
#     executor = FLoops.SequentialEx(),
#     # executor = FLoops.ThreadedEx(),
#     save_results = false,
#     save_checkpoints = false,
#     checkpoint_dir = "",
#     verbose = false
# )
#

#%%
begin
    println("Original StructVector (Sequential)")
    a = @benchmark evaluate_gridsearch_scenarios(
        $missing_scenarios,
        $data_arrs;
        # executor = FLoops.ThreadedEx(),
        executor = FLoops.SequentialEx(),
        save_results = false,
        save_checkpoints = false,
        checkpoint_dir = "",
        verbose = false
    )
    display(a)

    println("\n\nOriginal StructVector (Multithreaded)")
    a = @benchmark evaluate_gridsearch_scenarios(
        $missing_scenarios,
        $data_arrs;
        executor = FLoops.ThreadedEx(),
        save_results = false,
        save_checkpoints = false,
        checkpoint_dir = "",
        verbose = false
    )
    display(a)


    println("\n\nStructVector without filters")
    a = @benchmark evaluate_gridsearch_scenarios2(
        $missing_scenarios,
        $data_arrs;
        # executor = FLoops.ThreadedEx(),
        executor = FLoops.SequentialEx(),
        save_results = false,
        save_checkpoints = false,
        checkpoint_dir = "",
        verbose = false
    )
    display(a)

    println("\n\nStructVector with cached data")
    a = @benchmark evaluate_gridsearch_scenarios_optimized(
        $missing_scenarios,
        $data_arrs;
        # executor = FLoops.ThreadedEx(),
        executor = FLoops.SequentialEx(),
        save_results = false,
        save_checkpoints = false,
        checkpoint_dir = "",
        verbose = false
    )
    display(a)

    println("\n\nStructVector with Bumper")
    a = @benchmark evaluate_gridsearch_scenarios_bumper(
        $missing_scenarios,
        $data_arrs;
        # executor = FLoops.ThreadedEx(),
        executor = FLoops.SequentialEx(),
        save_results = false,
        save_checkpoints = false,
        checkpoint_dir = "",
        verbose = false
    )
    display(a)

    println("\n\nDataFrames")
    a = @benchmark ews_hyperparam_gridsearch(
        $specification_vecs,
        $data_arrs;
        force = true,
        save_results = false,
        verbose = false,
        disable_time_check = true
    )
    display(a)

    for n in [20, 50, 100, 150]
        println("\n\nMultistart Multivariate Optimization ($n Sobol points)")
        local a = @benchmark ews_multistart_optimization(
            $multivariate_optimization_vecs,
            $data_arrs;
            quantile_bounds = (0.5, 0.99),
            consecutive_bounds = (2.0, 30.0),
            n_sobol_points = $n,
            maxeval = 1000,
            executor = FLoops.SequentialEx(),
            force = true,
            return_df = true,
            save_results = false,
            verbose = false,
            disable_time_check = true,
            xtol_rel = 1.0e-3,
            xtol_abs = 1.0e-3,
            ftol_rel = 1.0e-4,
        )
        display(a)

        println("\n\nMultistart Univariate Optimization ($n Sobol points)")
        local a = @benchmark CSDNoise.evaluate_gridsearch_scenarios_multistart(
            $univariate_optim_missing_scenarios,
            $data_arrs;
            quantile_bounds = (0.5, 0.99),
            save_results = true,
            save_checkpoints = false,
            verbose = false,
            # Optimization configuration
            n_sobol_points = $n,
            local_algorithm = NLopt.LN_BOBYQA,
            maxeval = 1000,
            xtol_rel = 1.0e-3,
            xtol_abs = 1.0e-3,
            ftol_rel = 1.0e-4,
        )
        display(a)
    end
end

#%%
original_res = evaluate_gridsearch_scenarios(
    missing_scenarios,
    data_arrs;
    # executor = floops.threadedex(),
    executor = FLoops.SequentialEx(),
    save_results = false,
    save_checkpoints = false,
    checkpoint_dir = "",
    verbose = false
);

original_res_threaded = evaluate_gridsearch_scenarios(
    missing_scenarios,
    data_arrs;
    executor = FLoops.ThreadedEx(),
    save_results = false,
    save_checkpoints = false,
    checkpoint_dir = "",
    verbose = false
);

#%%
original_res2 = evaluate_gridsearch_scenarios2(
    missing_scenarios,
    data_arrs;
    # executor = FLoops.ThreadedEx(),
    executor = FLoops.SequentialEx(),
    save_results = false,
    save_checkpoints = false,
    checkpoint_dir = "",
    verbose = false
);

#%%
optimized_res = evaluate_gridsearch_scenarios_optimized(
    missing_scenarios,
    data_arrs;
    # executor = FLoops.ThreadedEx(),
    executor = FLoops.SequentialEx(),
    save_results = false,
    save_checkpoints = false,
    checkpoint_dir = "",
    verbose = false
);

#%%
bumper_res = evaluate_gridsearch_scenarios_bumper(
    missing_scenarios,
    data_arrs;
    # executor = FLoops.ThreadedEx(),
    executor = FLoops.SequentialEx(),
    save_results = false,
    save_checkpoints = false,
    checkpoint_dir = "",
    verbose = false
);

#%%
df_results = ews_hyperparam_gridsearch(
    specification_vecs,
    data_arrs;
    force = true,
    save_results = false,
    return_df = false,
    verbose = false,
    disable_time_check = true
);

#%%
# Run multistart optimization with different configurations for comparison
multistart_fast_results = ews_multistart_optimization(
    multistart_specification_vecs,
    data_arrs;
    quantile_bounds = (0.5, 0.99),
    consecutive_bounds = (2.0, 30.0),
    n_sobol_points = 20,
    maxeval = 100,
    executor = FLoops.SequentialEx(),
    force = true,
    return_df = true,
    save_results = false,
    verbose = false,
    disable_time_check = true
);

multistart_balanced_results = ews_multistart_optimization(
    multistart_specification_vecs,
    data_arrs;
    quantile_bounds = (0.5, 0.99),
    consecutive_bounds = (2.0, 30.0),
    n_sobol_points = 50,
    maxeval = 200,
    executor = FLoops.SequentialEx(),
    force = true,
    return_df = true,
    save_results = false,
    verbose = false,
    disable_time_check = true
);

multistart_thorough_results = ews_multistart_optimization(
    multistart_specification_vecs,
    data_arrs;
    quantile_bounds = (0.5, 0.99),
    consecutive_bounds = (2.0, 30.0),
    n_sobol_points = 100,
    maxeval = 500,
    executor = FLoops.SequentialEx(),
    force = true,
    return_df = true,
    save_results = false,
    verbose = false,
    disable_time_check = true
);


#%%
original_res == original_res2
original_res == original_res_threaded
original_res == optimized_res
original_res == bumper_res

#%%
# Compare results - check that optimal results are equivalent
function compare_all_results(df_results, sv_results)
    # Convert StructVector to DataFrame for comparison
    sv_df = DataFrame(sv_results)

    # Rename StructVector columns to match DataFrame if needed
    sv_df_renamed = rename(
        sv_df, :threshold_quantile => :ews_threshold_quantile,
        :consecutive_thresholds => :ews_consecutive_thresholds
    )

    # Standardize burnin periods: convert all Years to Days in DataFrame results
    df_results_standardized = copy(df_results)
    for i in 1:nrow(df_results_standardized)
        burnin = df_results_standardized[i, :ews_threshold_burnin]
        if isa(burnin, Dates.Year)
            # Convert Year to Days (assuming 365 days per year)
            df_results_standardized[i, :ews_threshold_burnin] = Dates.Day(round(Dates.days(burnin)))
        end
    end

    # Method: Compare by creating unique scenario keys using string representations
    function create_scenario_key(row)
        return (
            string(row.noise_specification),
            string(row.test_specification),
            row.percent_tested,
            string(row.ews_metric_specification),
            string(row.ews_enddate_type),
            row.ews_metric,
            string(row.ews_threshold_window),
            string(row.ews_threshold_burnin),
            row.ews_threshold_quantile,
            row.ews_consecutive_thresholds,
        )
    end

    # Create dictionaries mapping scenario keys to results and original rows
    df_dict = Dict()
    for row in eachrow(df_results_standardized)
        key = create_scenario_key(row)
        df_dict[key] = (
            accuracy = row.accuracy,
            sensitivity = row.sensitivity,
            specificity = row.specificity,
            row = row,
        )
    end

    sv_dict = Dict()
    for row in eachrow(sv_df_renamed)
        key = create_scenario_key(row)
        sv_dict[key] = (
            accuracy = row.accuracy,
            sensitivity = row.sensitivity,
            specificity = row.specificity,
            row = row,
        )
    end

    # Get all unique scenario keys from both datasets
    df_keys = Set(keys(df_dict))
    sv_keys = Set(keys(sv_dict))
    all_keys = union(df_keys, sv_keys)

    # Debug: Check if keys are the same
    println("Debug: DataFrame rows: $(nrow(df_results_standardized)), unique keys: $(length(df_keys))")
    println("Debug: StructVector rows: $(nrow(sv_df_renamed)), unique keys: $(length(sv_keys))")
    println("Debug: Union keys count: $(length(all_keys))")
    println("Debug: Keys only in DataFrame: $(length(setdiff(df_keys, sv_keys)))")
    println("Debug: Keys only in StructVector: $(length(setdiff(sv_keys, df_keys)))")

    # Check for duplicate keys within each dataset
    if length(df_keys) != nrow(df_results_standardized)
        println("Warning: DataFrame has duplicate keys! $(nrow(df_results_standardized)) rows but $(length(df_keys)) unique keys")
    end
    if length(sv_keys) != nrow(sv_df_renamed)
        println("Warning: StructVector has duplicate keys! $(nrow(sv_df_renamed)) rows but $(length(sv_keys)) unique keys")
    end

    if length(df_keys) != length(sv_keys) || length(all_keys) != length(df_keys)
        println("Debug: Key mismatch detected!")

        # Show a few examples of differing keys
        only_in_df = setdiff(df_keys, sv_keys)
        only_in_sv = setdiff(sv_keys, df_keys)

        if !isempty(only_in_df)
            println("Debug: First key only in DataFrame:")
            println(first(only_in_df))
        end

        if !isempty(only_in_sv)
            println("Debug: First key only in StructVector:")
            println(first(only_in_sv))
        end

        # Let's also check if the column names are different
        println("Debug: DataFrame columns: $(names(df_results_standardized))")
        println("Debug: StructVector columns: $(names(sv_df_renamed))")

        # Show first few keys from each to see the difference
        println("Debug: First DataFrame key:")
        if !isempty(df_keys)
            first_df_key = first(df_keys)
            println(first_df_key)
        end

        println("Debug: First StructVector key:")
        if !isempty(sv_keys)
            first_sv_key = first(sv_keys)
            println(first_sv_key)
        end
    end

    # Create comparison DataFrame
    comparison_rows = []

    for key in all_keys
        # Extract scenario information from key
        noise_spec_str, test_spec_str, percent_tested, ews_metric_spec_str,
            ews_enddate_type_str, ews_metric, ews_threshold_window_str,
            ews_threshold_burnin_str, ews_threshold_quantile, ews_consecutive_thresholds = key

        # Check if scenario exists in both datasets
        in_dataframe = haskey(df_dict, key)
        in_structvector = haskey(sv_dict, key)

        # Get results for each method (if available)
        df_accuracy = in_dataframe ? df_dict[key].accuracy : missing
        df_sensitivity = in_dataframe ? df_dict[key].sensitivity : missing
        df_specificity = in_dataframe ? df_dict[key].specificity : missing

        sv_accuracy = in_structvector ? sv_dict[key].accuracy : missing
        sv_sensitivity = in_structvector ? sv_dict[key].sensitivity : missing
        sv_specificity = in_structvector ? sv_dict[key].specificity : missing

        # Calculate differences and matches (only if both values exist)
        accuracy_diff = (in_dataframe && in_structvector) ? abs(df_accuracy - sv_accuracy) : missing
        sensitivity_diff = (in_dataframe && in_structvector) ? abs(df_sensitivity - sv_sensitivity) : missing
        specificity_diff = (in_dataframe && in_structvector) ? abs(df_specificity - sv_specificity) : missing

        accuracy_match = ismissing(accuracy_diff) ? missing : accuracy_diff < 1.0e-10
        sensitivity_match = ismissing(sensitivity_diff) ? missing : sensitivity_diff < 1.0e-10
        specificity_match = ismissing(specificity_diff) ? missing : specificity_diff < 1.0e-10

        # Overall match for this scenario
        scenario_match = if in_dataframe && in_structvector
            !ismissing(accuracy_match) && !ismissing(sensitivity_match) && !ismissing(specificity_match) &&
                accuracy_match && sensitivity_match && specificity_match
        else
            false  # Can't match if not in both datasets
        end

        push!(
            comparison_rows, (
                # Scenario identifiers
                noise_specification = noise_spec_str,
                test_specification = test_spec_str,
                percent_tested = percent_tested,
                ews_metric_specification = ews_metric_spec_str,
                ews_enddate_type = ews_enddate_type_str,
                ews_metric = ews_metric,
                ews_threshold_window = ews_threshold_window_str,
                ews_threshold_burnin = ews_threshold_burnin_str,
                ews_threshold_quantile = ews_threshold_quantile,
                ews_consecutive_thresholds = ews_consecutive_thresholds,

                # Presence indicators
                in_dataframe = in_dataframe,
                in_structvector = in_structvector,

                # DataFrame results
                df_accuracy = df_accuracy,
                df_sensitivity = df_sensitivity,
                df_specificity = df_specificity,

                # StructVector results
                sv_accuracy = sv_accuracy,
                sv_sensitivity = sv_sensitivity,
                sv_specificity = sv_specificity,

                # Differences
                accuracy_diff = accuracy_diff,
                sensitivity_diff = sensitivity_diff,
                specificity_diff = specificity_diff,

                # Match indicators
                accuracy_match = accuracy_match,
                sensitivity_match = sensitivity_match,
                specificity_match = specificity_match,
                scenario_match = scenario_match,
            )
        )
    end
    comparison_df = DataFrame(comparison_rows)
    println("Number of scenarios that don't match accuracy values: $(nrow(comparison_df) - sum(comparison_df.accuracy_match))")
    println("\tMean accuracy diff: $(mean(comparison_df.accuracy_diff))")
    println("Number of scenarios that don't match sensitivity values: $(nrow(comparison_df) - sum(comparison_df.sensitivity_match))")
    println("\tMean sensitivity diff: $(mean(comparison_df.sensitivity_diff))")
    println("Number of scenarios that don't match specificity values: $(nrow(comparison_df) - sum(comparison_df.specificity_match))")
    println("\tMean specificity diff: $(mean(comparison_df.specificity_diff))")

    return comparison_df
end

# Debug: Check sizes before comparison
println("Debug: Full DataFrame results: $(nrow(df_results)) rows")
println("Debug: Full StructVector results: $(length(original_res)) rows")

comparison_df = compare_all_results(df_results, original_res);

#%%
# Compare multistart optimization results with grid search
println("\n=== MULTISTART OPTIMIZATION RESULTS COMPARISON ===")
println("Grid Search (DataFrame) - Best accuracy: $(round(maximum(df_results.accuracy), digits = 4))")
println("Grid Search (StructVector) - Best accuracy: $(round(maximum(original_res.accuracy), digits = 4))")
println("Multistart Fast (20 Sobol) - Best accuracy: $(round(maximum(multistart_fast_results.accuracy), digits = 4))")
println("Multistart Balanced (50 Sobol) - Best accuracy: $(round(maximum(multistart_balanced_results.accuracy), digits = 4))")
println("Multistart Thorough (100 Sobol) - Best accuracy: $(round(maximum(multistart_thorough_results.accuracy), digits = 4))")

println("\nMean accuracies:")
println("Grid Search (DataFrame): $(round(mean(df_results.accuracy), digits = 4))")
println("Grid Search (StructVector): $(round(mean(original_res.accuracy), digits = 4))")
println("Multistart Fast: $(round(mean(multistart_fast_results.accuracy), digits = 4))")
println("Multistart Balanced: $(round(mean(multistart_balanced_results.accuracy), digits = 4))")
println("Multistart Thorough: $(round(mean(multistart_thorough_results.accuracy), digits = 4))")

println("\nOptimal parameters found:")
best_grid_idx = argmax(df_results.accuracy)
best_grid_row = df_results[best_grid_idx, :]
println("Grid Search - Quantile: $(round(best_grid_row.ews_threshold_quantile, digits = 3)), Consecutive: $(best_grid_row.ews_consecutive_thresholds)")

best_fast_idx = argmax(multistart_fast_results.accuracy)
best_fast_row = multistart_fast_results[best_fast_idx, :]
println("Multistart Fast - Quantile: $(round(best_fast_row.threshold_quantile, digits = 3)), Consecutive: $(best_fast_row.consecutive_thresholds)")

best_balanced_idx = argmax(multistart_balanced_results.accuracy)
best_balanced_row = multistart_balanced_results[best_balanced_idx, :]
println("Multistart Balanced - Quantile: $(round(best_balanced_row.threshold_quantile, digits = 3)), Consecutive: $(best_balanced_row.consecutive_thresholds)")

best_thorough_idx = argmax(multistart_thorough_results.accuracy)
best_thorough_row = multistart_thorough_results[best_thorough_idx, :]
println("Multistart Thorough - Quantile: $(round(best_thorough_row.threshold_quantile, digits = 3)), Consecutive: $(best_thorough_row.consecutive_thresholds)")
