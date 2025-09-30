using DataFrames: DataFrame, nrow

export dataframe_row_to_scenario

"""
    dataframe_row_to_scenario(row)

Convert a DataFrame row to an OptimizationScenario struct.
"""
function dataframe_row_to_scenario(row::DataFrameRow)
    return OptimizationScenario(
        row.ensemble_specification,
        row.null_specification,
        row.noise_specification,
        row.test_specification,
        row.percent_tested,
        row.ews_metric_specification,
        row.ews_enddate_type,
        row.ews_threshold_window,
        row.ews_threshold_burnin,
        row.ews_metric,
    )
end
