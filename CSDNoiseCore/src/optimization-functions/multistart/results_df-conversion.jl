export df_to_structvector

"""
    df_to_structvector(df, ::Type{OptimizationResult})

Convert DataFrame to StructVector{OptimizationResult} for migration from old format.
"""
function df_to_structvector(df::DF.DataFrame, ::Type{OptimizationResult})
    return StructVector(
        map(eachrow(df)) do row
            OptimizationResult(
                row.noise_specification,
                row.test_specification,
                row.percent_tested,
                row.ews_metric_specification,
                row.ews_enddate_type,
                row.ews_threshold_window,
                row.ews_threshold_burnin,
                row.ews_metric,
                row.threshold_quantile,
                row.consecutive_thresholds,
                row.accuracy,
                row.sensitivity,
                row.specificity
            )
        end
    )
end
