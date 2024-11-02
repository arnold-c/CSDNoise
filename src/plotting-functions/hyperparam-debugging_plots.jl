using UnPack: @unpack

function hyperparam_debugging_Reff_plot(
    inc_vec,
    null_inc_vec,
    Reff_vec,
    null_Reff_vec,
    Reff_thresholds,
    null_Reff_thresholds,
    outbreak_bounds,
    null_outbreak_bounds,
    ewsmetric,
    test_positive_vec,
    null_test_positive_vec,
    ews_vals_vec,
    null_ews_vals_vec,
    exceed_thresholds_vec,
    null_exceed_thresholds_vec,
    threshold_percentiles_vec,
    null_threshold_percentiles_vec,
    detection_index_vec,
    null_detection_index_vec,
    ensemble_time_specification;
    plottitle = "",
    xlims = (0, 12),
    kwargs...,
)
    ews_specification = ews_vals_vec.ews_specification
    @unpack aggregation = ews_specification
    aggregation_int = Dates.days(aggregation)

    aggregated_inc_vec = aggregate_timeseries(inc_vec, aggregation)
    aggregated_null_inc_vec = aggregate_timeseries(null_inc_vec, aggregation)

    aggregated_Reff_vec = aggregate_Reff_vec(Reff_vec, aggregation)
    aggregated_Reff_thresholds = Reff_thresholds .÷ aggregation_int

    aggregated_null_Reff_vec = aggregate_Reff_vec(
        null_Reff_vec, aggregation
    )
    aggregated_null_Reff_thresholds =
        null_Reff_thresholds .÷ aggregation_int

    aggregated_outbreak_bounds =
        outbreak_bounds[(outbreak_bounds[:, 4] .== 1), [1, 2]] .÷
        aggregation_int

    aggregated_null_outbreak_bounds =
        null_outbreak_bounds[(null_outbreak_bounds[:, 4] .== 1), [1, 2]] .÷
        aggregation_int

    aggregated_test_vec = aggregate_timeseries(test_positive_vec, aggregation)

    aggregated_test_movingavg_vec = zeros(
        Int64, size(aggregated_test_vec)
    )

    calculate_movingavg!(
        aggregated_test_movingavg_vec,
        aggregated_test_vec,
        7,
    )

    aggregated_null_test_vec = aggregate_timeseries(
        null_test_positive_vec, aggregation
    )

    aggregated_null_test_movingavg_vec = zeros(
        Int64, size(aggregated_null_test_vec)
    )

    calculate_movingavg!(
        aggregated_null_test_movingavg_vec,
        aggregated_null_test_vec,
        7,
    )

    return Reff_ews_plot(
        aggregated_inc_vec,
        aggregated_null_inc_vec,
        aggregated_Reff_vec,
        aggregated_null_Reff_vec,
        aggregated_Reff_thresholds,
        aggregated_null_Reff_thresholds,
        ews_vals_vec,
        null_ews_vals_vec,
        ewsmetric,
        aggregated_outbreak_bounds,
        aggregated_null_outbreak_bounds,
        exceed_thresholds_vec,
        null_exceed_thresholds_vec,
        threshold_percentiles_vec,
        null_threshold_percentiles_vec,
        detection_index_vec,
        null_detection_index_vec,
        ensemble_time_specification;
        plottitle = plottitle,
        xlims = xlims,
        kwargs...,
    )
end
