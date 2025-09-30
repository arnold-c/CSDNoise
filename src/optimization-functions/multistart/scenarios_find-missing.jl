using DataFrames: DataFrame, nrow, antijoin

export find_missing_scenarios

"""
    find_missing_scenarios(all_scenarios_df, completed_results_df)

Find scenarios that haven't been computed yet using antijoin.
Returns DataFrame of missing scenarios.
"""
function find_missing_scenarios(all_scenarios_df::DataFrame, completed_results_df::DataFrame)
    if nrow(completed_results_df) == 0
        return all_scenarios_df
    end

    # Define scenario columns for comparison
    scenario_cols = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
    ]

    # Find scenarios that are in all_scenarios but not in completed_results
    missing_scenarios_df = antijoin(
        all_scenarios_df,
        completed_results_df;
        on = scenario_cols
    )

    return missing_scenarios_df
end
