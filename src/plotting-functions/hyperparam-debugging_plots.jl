function hyperparam_debugging_Reff_plot(
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_arr,
    null_single_Reff_arr,
    ensemble_single_Reff_thresholds_vec,
    null_single_Reff_thresholds_vec,
    ensemble_single_periodsum_vecs,
    null_single_periodsum_vecs,
    ewsmetric,
    vec_of_testarr,
    vec_of_null_testarr,
    vec_of_ews_vals_vec,
    vec_of_null_ews_vals_vec,
    vec_of_exceed_thresholds,
    vec_of_null_exceed_thresholds,
    vec_of_detection_index_vec,
    vec_of_null_detection_index_vec,
    ensemble_time_specification;
    selected_sim = 1,
    test_index = 1,
    xlims = (0, 12),
    ylims_metric = (nothing, nothing),
    ylims_inc = (nothing, nothing),
)
    aggregated_inc_vec = aggregate_timeseries(
        @view(ensemble_single_incarr[:, 1, selected_sim]),
        7,
    )

    aggregated_outbreak_status_vec = aggregate_thresholds_vec(
        @view(ensemble_single_incarr[:, 3, selected_sim]),
        7,
    )

    aggregated_null_inc_vec = aggregate_timeseries(
        @view(null_single_incarr[:, 1, selected_sim]),
        7,
    )

    aggregated_null_outbreak_status_vec = aggregate_thresholds_vec(
        @view(null_single_incarr[:, 3, selected_sim]),
        7,
    )

    aggregated_Reff_vec = aggregate_Reff_vec(
        @view(ensemble_single_Reff_arr[:, selected_sim]),
        7,
    )

    aggregated_Reff_thresholds_arr =
        ensemble_single_Reff_thresholds_vec[selected_sim] .÷ 7

    aggregated_null_Reff_vec = aggregate_Reff_vec(
        @view(null_single_Reff_arr[:, selected_sim]),
        7,
    )

    aggregated_null_Reff_thresholds_arr =
        null_single_Reff_thresholds_vec[selected_sim] .÷ 7

    aggregated_outbreak_thresholds_arr =
        ensemble_single_periodsum_vecs[selected_sim][
            (ensemble_single_periodsum_vecs[selected_sim][:, 4] .== 1),
            [1, 2],
        ] .÷ 7

    aggregated_null_outbreak_thresholds_arr =
        null_single_periodsum_vecs[selected_sim][
            (null_single_periodsum_vecs[selected_sim][:, 4] .== 1),
            [1, 2],
        ] .÷ 7

    aggregated_test_vec = aggregate_timeseries(
        @view(vec_of_testarr[test_index][:, 5, selected_sim]),
        7,
    )

    aggregated_test_movingavg_vec = zeros(
        Int64, size(aggregated_test_vec)
    )

    calculate_movingavg!(
        aggregated_test_movingavg_vec,
        aggregated_test_vec,
        7,
    )

    aggregated_null_test_vec = aggregate_timeseries(
        @view(vec_of_null_testarr[test_index][:, 5, selected_sim]),
        7,
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
        aggregated_Reff_thresholds_arr,
        aggregated_null_Reff_thresholds_arr,
        vec_of_ews_vals_vec[test_index][selected_sim],
        vec_of_null_ews_vals_vec[test_index][selected_sim],
        ewsmetric,
        aggregated_outbreak_thresholds_arr,
        aggregated_null_outbreak_thresholds_arr,
        vec(vec_of_exceed_thresholds[test_index][selected_sim, 1]),
        vec(vec_of_null_exceed_thresholds[test_index][selected_sim, 1]),
        vec_of_detection_index_vec[test_index][selected_sim],
        vec_of_null_detection_index_vec[test_index][selected_sim],
        ensemble_time_specification;
        xlims = xlims,
        ylims_metric = ylims_metric,
        ylims_inc = ylims_inc,
    )
end
