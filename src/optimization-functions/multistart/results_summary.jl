export summarize_optimization_results,
    create_results_dataframe

"""
    summarize_optimization_results(results_df)

Provide a summary of optimization results for user feedback.
"""
function summarize_optimization_results(results_df::DF.DataFrame)
    if DF.nrow(results_df) == 0
        println(StyledStrings.styled"{yellow:No results to summarize}")
        return
    end

    println(StyledStrings.styled"{green:=== Optimization Results Summary ===}")
    println(StyledStrings.styled"Total scenarios: {yellow:$(nrow(results_df))}")

    # Accuracy statistics
    acc_stats = (
        mean = StatsBase.mean(results_df.accuracy),
        median = StatsBase.median(results_df.accuracy),
        min = minimum(results_df.accuracy),
        max = maximum(results_df.accuracy),
    )

    println(
        StyledStrings.styled"Accuracy - Mean: {blue:$(round(acc_stats.mean, digits=3))}, " *
            "Median: {blue:$(round(acc_stats.median, digits = 3))}, " *
            "Range: {blue:$(round(acc_stats.min, digits = 3))} - {blue:$(round(acc_stats.max, digits = 3))}"
    )

    # Count by major categories
    noise_counts = DF.combine(
        DF.groupby(results_df, :noise_specification),
        nrow => :count
    )
    println(StyledStrings.styled"Noise specifications: {cyan:$(nrow(noise_counts))} types")

    metric_counts = DF.combine(
        DF.groupby(results_df, :ews_metric),
        nrow => :count
    )
    println(StyledStrings.styled"EWS metrics: {cyan:$(nrow(metric_counts))} types")

    # Best performing scenarios
    best_idx = argmax(results_df.accuracy)
    best_accuracy = results_df.accuracy[best_idx]
    return println(
        StyledStrings.styled"Best accuracy: {green:$(round(best_accuracy, digits=4))} " *
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
    scenarios_df = DF.DataFrame(scenarios)
    results_df = DF.DataFrame(results)
    return DF.hcat(scenarios_df, results_df)
end
