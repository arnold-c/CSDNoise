export merge_results_safely

"""
    merge_results_safely(df1, df2)

Safely merge two results DataFrames, removing duplicates based on scenario parameters.
"""
function merge_results_safely(df1::DF.DataFrame, df2::DF.DataFrame)
    if DF.nrow(df1) == 0
        return df2
    elseif DF.nrow(df2) == 0
        return df1
    end

    # Define scenario columns for duplicate detection
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

    # Combine DataFrames
    combined_df = DF.vcat(df1, df2; cols = :union)

    # Remove duplicates, keeping first occurrence
    # (assumes df1 contains more recent/authoritative results)
    unique_df = unique(combined_df, scenario_cols)

    return unique_df
end
