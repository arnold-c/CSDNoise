using DataFrames: DataFrame
using Dates: Dates

export create_empty_results_dataframe

"""
    create_empty_results_dataframe()

Create an empty DataFrame with the correct structure for multistart optimization results.
"""
function create_empty_results_dataframe()
    return DataFrame(
        noise_specification = NoiseSpecification[],
        test_specification = eltype(IndividualTestSpecification)[],
        percent_tested = Float64[],
        ews_metric_specification = EWSMetricSpecification[],
        ews_enddate_type = EWSEndDateType[],
        ews_threshold_window = Union{Type{ExpandingThresholdWindow}, Type{RollingThresholdWindow}}[],
        ews_metric = String[],
        ews_threshold_quantile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_threshold_burnin = Union{Dates.Day, Dates.Year}[],
        accuracy = Float64[],
        sensitivity = Float64[],
        specificity = Float64[]
    )
end
