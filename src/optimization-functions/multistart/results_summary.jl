using DataFrames: DataFrame

export summarize_optimization_results,
    create_results_dataframe

"""
    summarize_optimization_results(results_df)

Provide a summary of optimization results for user feedback.
"""
function summarize_optimization_results(results_df::DataFrame)
    if nrow(results_df) == 0
        println(styled"{yellow:No results to summarize}")
        return
    end

    println(styled"{green:=== Optimization Results Summary ===}")
    println(styled"Total scenarios: {yellow:$(nrow(results_df))}")

    # Accuracy statistics
    acc_stats = (
        mean = mean(results_df.accuracy),
        median = median(results_df.accuracy),
        min = minimum(results_df.accuracy),
        max = maximum(results_df.accuracy),
    )

    println(
        styled"Accuracy - Mean: {blue:$(round(acc_stats.mean, digits=3))}, " *
            "Median: {blue:$(round(acc_stats.median, digits = 3))}, " *
            "Range: {blue:$(round(acc_stats.min, digits = 3))} - {blue:$(round(acc_stats.max, digits = 3))}"
    )

    # Count by major categories
    noise_counts = combine(groupby(results_df, :noise_specification), nrow => :count)
    println(styled"Noise specifications: {cyan:$(nrow(noise_counts))} types")

    metric_counts = combine(groupby(results_df, :ews_metric), nrow => :count)
    println(styled"EWS metrics: {cyan:$(nrow(metric_counts))} types")

    # Best performing scenarios
    best_idx = argmax(results_df.accuracy)
    best_accuracy = results_df.accuracy[best_idx]
    return println(
        styled"Best accuracy: {green:$(round(best_accuracy, digits=4))} " *
            "with metric: {yellow:$(results_df.ews_metric[best_idx])}"
    )
end

"""
    create_results_dataframe(results, scenarios)

Convert optimization results to DataFrame format.
"""
function create_results_dataframe(
        scenarios::StructVector{OptimizationScenario},
        results::StructVector{OptimizedValues},
    )
    @assert length(scenarios) == length(results)
    scenarios_df = DataFrame(scenarios)
    results_df = DataFrame(results)
    return hcat(scenarios_df, results_df)
end
