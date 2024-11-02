using DataFrames: DataFrames

function create_optimal_ews_plots(
    optimal_ews_df,
    ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec,
    tiebreaker_preference_vec;
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
    ],
    ews_metrics = [
        "autocorrelation",
        "autocovariance",
        "coefficient_of_variation",
        "index_of_dispersion",
        "kurtosis",
        "mean",
        "skewness",
        "variance",
    ],
    force = false,
    base_plotpath = joinpath(
        plotsdir(), "ensemble", "ews-hyperparam-optimization"
    ),
)
    grouped_optimal_ews_df = DataFrames.groupby(
        optimal_ews_df,
        [
            :noise_specification,
            :percent_tested,
            :ews_metric_specification,
            :ews_enddate_type,
            :ews_threshold_burnin,
        ],
    )

    for (gdf, tiebreaker_preference) in
        Iterators.product(grouped_optimal_ews_df, tiebreaker_preference_vec)
        optimal_heatmap_df = optimal_ews_heatmap_df(
            gdf;
            tiebreaker_preference = tiebreaker_preference,
            optimal_grouping_parameters = optimal_grouping_parameters,
        )

        @assert nrow(optimal_heatmap_df) ==
            length(ews_metrics) *
                unique(grouped_optimal_ews_df.test_specification)

        noise_plotdir = getdirpath(noise_specification)
        percent_tested = optimal_ews_df[1, :percent_tested]
        noise_plotpath = joinpath(
            base_plotpath,
            noise_plotdir,
            "percent-tested-$(percent_tested)",
        )
        mkpath(noise_plotpath)

        noise_descripton = get_noise_magnitude_description(noise_specification)

        plotpath = joinpath(
            noise_plotpath,
            "ews-heatmap_tiebreaker-$(tiebreaker_preference).png",
        )

        if !isfile(plotpath) || force
            ews_metric_specification = optimal_heatmap_df.ews_metric_specification[1]
            @unpack method, aggregation, bandwidth, lag =
                ews_metric_specification
            ews_enddate_type = optimal_heatmap_df.ews_enddate_type[1]

            heatmap_subtitle =
                "Noise: $(noise_descripton), Percent Tested: $(percent_tested), $(ews_enddate_type)" *
                "\n$(get_ews_metric_specification_description(ews_metric_specification))" *
                "\nP = Percentile Threshold, C = Consecutive Thresholds, S = Specificity"

            optimal_heatmap_plot = optimal_ews_heatmap_plot(
                optimal_heatmap_df; subtitle = heatmap_subtitle
            )

            save(
                plotpath,
                optimal_heatmap_plot,
            )
        end

        # test_description = get_test_description(individual_test_specification)
        #
        # survival_plotdir = joinpath(
        #     noise_plotdir,
        #     "survival",
        # )
        #
        # for ews_metric in ews_metrics
        #     plotpath = joinpath(
        #         survival_plotdir,
        #         "ews_survival_$(ews_metric)_ews_agg-$(aggregation).png",
        #     )
        #     if !isfile(plotpath) || force
        #         survival_plot = simulate_and_plot_ews_survival(
        #             subset_df,
        #             ews_metric_specification,
        #             ews_threshold_burnin,
        #             ensemble_specification,
        #             individual_test_specification,
        #             ensemble_single_incarr,
        #             null_single_incarr,
        #             ensemble_single_Reff_thresholds_vec;
        #             ews_metric = ews_metric,
        #         )
        #         save(
        #             plotpath,
        #             survival_plot,
        #         )
        #     end
        # end
    end
    return nothing
end
