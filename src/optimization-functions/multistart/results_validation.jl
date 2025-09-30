using DataFrames: DataFrame, names

export validate_results_integrity

"""
    validate_results_integrity(results_df)

Validate the integrity of results DataFrame.
"""
function validate_results_integrity(results_df::DataFrame)
    required_cols = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
        :ews_threshold_quantile,
        :ews_consecutive_thresholds,
        :accuracy,
        :sensitivity,
        :specificity,
    ]

    # Check required columns exist
    missing_cols = setdiff(required_cols, names(results_df))
    if !isempty(missing_cols)
        @warn "Missing required columns: $missing_cols"
        return false
    end

    # Check for missing values in critical columns
    for col in [:accuracy, :sensitivity, :specificity]
        if any(ismissing, results_df[!, col])
            @warn "Missing values found in column: $col"
            return false
        end
    end

    # Check accuracy bounds
    if any(x -> x < 0 || x > 1, results_df.accuracy)
        @warn "Accuracy values outside [0,1] range found"
        return false
    end

    @info "Results DataFrame validation passed"
    return true
end
