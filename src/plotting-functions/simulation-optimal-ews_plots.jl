using DataFrames: DataFrames
using ProgressMeter

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
    base_plotpath = joinpath(plotsdir(), "ensemble"),
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

    @showprogress for (gdf, tiebreaker_preference) in
                      Iterators.product(
        grouped_optimal_ews_df, tiebreaker_preference_vec
    )
        optimal_heatmap_df = optimal_ews_heatmap_df(
            gdf;
            tiebreaker_preference = tiebreaker_preference,
            optimal_grouping_parameters = optimal_grouping_parameters,
        )

        @assert nrow(optimal_heatmap_df) ==
            length(ews_metrics) *
                length(unique(optimal_ews_df.test_specification))

        burnin_time = ensemble_specification.time_parameters.burnin
        @unpack min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_vaccination_coverage,
        max_vaccination_coverage =
            ensemble_specification.dynamics_parameter_specification

        ensemble_vax_plotpath = joinpath(
            base_plotpath,
            "burnin-time-$(burnin_time)",
            "min-burnin-vax_$(min_burnin_vaccination_coverage)",
            "max-burnin-vax_$(max_burnin_vaccination_coverage)",
            "min-vax_$(min_vaccination_coverage)",
            "max-vax_$(max_vaccination_coverage)",
        )

        percent_tested = optimal_ews_df[1, :percent_tested]
        noise_specification = optimal_heatmap_df.noise_specification[1]
        noise_dirpath = getdirpath(noise_specification)

        noise_plotdir = joinpath(
            ensemble_vax_plotpath,
            noise_dirpath,
            "percent-tested_$(percent_tested)",
            # "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)",
        )
        ews_metric_specification = optimal_heatmap_df.ews_metric_specification[1]
        ews_enddate_type = optimal_heatmap_df.ews_enddate_type[1]
        ews_enddate_type_str = split(string(ews_enddate_type), "::")[1]
        ews_threshold_burnin = optimal_heatmap_df.ews_threshold_burnin[1]

        ews_plotdir = joinpath(
            noise_plotdir,
            "optimal-heatmaps",
            ews_metric_specification.dirpath,
            ews_enddate_type_str,
        )
        mkpath(ews_plotdir)

        noise_descripton = get_noise_magnitude_description(noise_specification)

        plotpath = joinpath(
            ews_plotdir,
            "ews-heatmap_tiebreaker-$(tiebreaker_preference).png",
        )

        if !isfile(plotpath) || force
            @unpack method, aggregation, bandwidth, lag =
                ews_metric_specification

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

        for test_specification in unique(optimal_heatmap_df.test_specification)
            test_description = get_test_description(test_specification)

            test_plotdir = joinpath(
                noise_plotdir,
                "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)",
            )

            survival_plotdir = joinpath(
                test_plotdir,
                "survival",
                ews_metric_specification.dirpath,
                ews_enddate_type_str,
                "tiebreaker-$(tiebreaker_preference)",
            )
            mkpath(survival_plotdir)

            for ews_metric in ews_metrics
                plotpath = joinpath(
                    survival_plotdir,
                    "ews_survival_$(ews_metric).png",
                )
                if !isfile(plotpath) || force
                    survival_plot = simulate_and_plot_ews_survival(
                        optimal_heatmap_df,
                        ews_metric_specification,
                        ews_threshold_burnin,
                        ensemble_specification,
                        test_specification,
                        ensemble_single_incarr,
                        null_single_incarr,
                        ensemble_single_Reff_thresholds_vec;
                        ews_metric = ews_metric,
                    )
                    save(
                        plotpath,
                        survival_plot,
                    )
                end
            end
        end
    end
    return nothing
end
