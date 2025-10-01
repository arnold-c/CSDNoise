export ews_as_df

"""
    ews_as_df(ews::EWSMetrics)

Convert EWSMetrics struct to a DataFrame for analysis and visualization.

# Arguments
- `ews::EWSMetrics`: EWS metrics struct to convert

# Returns
- `DataFrame`: DataFrame with columns for each EWS metric (excluding specification and tau values)

# Notes
- Excludes the ews_specification field and tau correlation coefficients
- Useful for data analysis, plotting, and export to other formats
- Column names correspond to the EWS metric names
"""
function ews_as_df(ews::EWSMetrics)
    metrics = filter(
        x -> x != :ews_specification && !contains(string(x), "tau"),
        propertynames(ews),
    )
    df =
        reduce(
        hcat,
        map(metric -> getproperty(ews, metric), metrics),
    ) |>
        array -> DF.DataFrame(array, [metrics...])
    return df
end
