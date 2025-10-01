using DataFrames: DataFrame, ByRow, subset
using UnPack: @unpack

export simulate_and_plot_ews_survival

function simulate_and_plot_ews_survival(
        optimal_heatmap_df,
        ews_metric_specification,
        ews_enddate_type,
        ews_threshold_burnin,
        ews_threshold_window,
        noise_specification,
        ensemble_specification,
        individual_test_specification,
        emergent_incidence_arr,
        null_incidence_arr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = "mean",
        plottitle = "Survival",
        subtitle = "",
    )
    subset_df = subset(
        optimal_heatmap_df,
        :ews_metric_specification => ByRow(==(ews_metric_specification)),
        :ews_enddate_type => ByRow(==(ews_enddate_type)),
        :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
        :ews_threshold_window => ByRow(==(ews_threshold_window)),
        :noise_specification => ByRow(==(noise_specification)),
    )

    return simulate_and_plot_ews_survival(
        subset_df,
        ews_metric_specification,
        ews_threshold_burnin,
        ensemble_specification,
        individual_test_specification,
        emergent_incidence_arr,
        null_incidence_arr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = ews_metric,
        plottitle = plottitle,
        subtitle = subtitle,
    )
end

function simulate_and_plot_ews_survival(
        subset_df,
        ews_metric_specification,
        ews_threshold_burnin,
        ensemble_specification,
        individual_test_specification,
        emergent_incidence_arr,
        null_incidence_arr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = "mean",
        plottitle = "Survival",
        subtitle = "",
    )
    test_subset_df = subset(
        subset_df,
        :test_specification => ByRow(==(individual_test_specification)),
    )

    survival_df = simulate_ews_survival_data(
        test_subset_df,
        ensemble_specification,
        emergent_incidence_arr,
        null_incidence_arr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = ews_metric,
    )[1]

    detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
        survival_df
    )

    @unpack aggregation = ews_metric_specification

    return ews_survival_plot(
        detection_survival_vecs,
        null_survival_vecs,
        survival_df.enddate;
        ews_aggregation = aggregation,
        burnin = ews_threshold_burnin,
        plottitle = plottitle,
        subtitle = subtitle,
    )
end
