export prepare_line_plot_df!

using DataFrames: DataFrames
using DrWatson: DrWatson
using StatsBase: StatsBase

function prepare_line_plot_df!(
        output_df,
        gdf::T1,
        ewsmetric = "mean",
        outcome = :accuracy,
        tests = [
            IndividualTestSpecification(0.8, 0.8, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
            IndividualTestSpecification(0.95, 0.95, 0),
            IndividualTestSpecification(0.96, 0.96, 0),
            IndividualTestSpecification(0.97, 0.97, 0),
            IndividualTestSpecification(0.98, 0.98, 0),
            IndividualTestSpecification(0.99, 0.99, 0),
            IndividualTestSpecification(1.0, 1.0, 0),
        ];
        tiebreaker_preference = "specificity",
    ) where {T1 <: Union{<:DataFrames.DataFrame, <:DataFrames.DataFrames.SubDataFrame}}
    tiebreaker_args = if tiebreaker_preference == "speed"
        (:ews_consecutive_thresholds, false)
    elseif tiebreaker_preference == "specificity"
        (:specificity, true)
    else
        error(
            "Invalid preference: $tiebreaker_preference. Please choose either \"speed\" or \"specificity\"."
        )
    end

    metric_df = subset(
        gdf,
        :ews_metric => ByRow(==(ewsmetric)),
        :test_specification => ByRow(in(tests)),
    )

    filtered_test_df =
        map(collect(groupby(metric_df, :test_specification))) do df
        sort(df, order(tiebreaker_args[1]; rev = tiebreaker_args[2]))[
            1, :,
        ]
    end |>
        x -> vcat(DataFrame.(x)...; cols = :union)

    filtered_test_df[!, :test_sens] =
        getproperty.(filtered_test_df.test_specification, :sensitivity)
    filtered_test_df[!, :test_spec] =
        getproperty.(filtered_test_df.test_specification, :sensitivity)
    filtered_test_df[!, :test_result_lag] =
        getproperty.(filtered_test_df.test_specification, :test_result_lag)
    sort!(
        filtered_test_df,
        [:test_sens, :test_spec, :test_result_lag];
        rev = [true, true, false],
    )

    append!(output_df, filtered_test_df; cols = :union)

    return nothing
end
