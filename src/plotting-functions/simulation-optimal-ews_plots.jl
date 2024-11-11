using DataFrames: DataFrames
using ProgressMeter
using Match: Match

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
    output_format = "png",
    pt_per_unit = 0.75,
    px_per_unit = 2,
)
    Match.@match output_format begin
        "pdf" || "svg" => include(srcdir("cairomakie-plotting-setup.jl"))
        "png" || "jpg" || "jpeg" => include(srcdir("makie-plotting-setup.jl"))
    end

    Match.@match output_format begin
        "pdf" => CairoMakie.activate!(; type = "pdf", pt_per_unit = pt_per_unit)
        "svg" => CairoMakie.activate!(; type = "svg", pt_per_unit = pt_per_unit)
        "png" || "jpg" || "jpeg" =>
            GLMakie.activate!(; px_per_unit = px_per_unit)
    end

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

    ngroups = length(grouped_optimal_ews_df) * length(tiebreaker_preference_vec)

    prog = Progress(ngroups)
    for (gdf, tiebreaker_preference) in
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
            "ews-threshold-burnin_$(ews_threshold_burnin)",
        )
        mkpath(ews_plotdir)

        noise_descripton = get_noise_magnitude_description(noise_specification)

        plotpath = joinpath(
            ews_plotdir,
            "ews-heatmap_tiebreaker-$(tiebreaker_preference).$(output_format)",
        )

        if !isfile(plotpath) || force
            @unpack method, aggregation, bandwidth, lag =
                ews_metric_specification

            heatmap_subtitle =
                "Noise: $(noise_descripton), Percent Tested: $(percent_tested), $(ews_enddate_type)" *
                "\n$(get_ews_metric_specification_description(ews_metric_specification)), Threshold Burnin: $(ews_threshold_burnin), Tiebreaker: $(tiebreaker_preference)" *
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

            for ews_metric in ews_metrics
                survival_plottitle = "Metric: $(ews_metric)"
                survival_subtitle =
                    "\nNoise: $(noise_descripton), Percent Tested: $(percent_tested), $(ews_enddate_type)" *
                    "\n$(get_ews_metric_specification_description(ews_metric_specification)), Threshold Burnin: $(ews_threshold_burnin), Tiebreaker: $(tiebreaker_preference), $(test_description)"

                test_plotdir = joinpath(
                    noise_plotdir,
                    "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)",
                )

                survival_plotdir = joinpath(
                    test_plotdir,
                    "survival",
                    ews_metric_specification.dirpath,
                    ews_enddate_type_str,
                    "ews-threshold-burnin_$(ews_threshold_burnin)",
                    "tiebreaker-$(tiebreaker_preference)",
                )
                mkpath(survival_plotdir)

                plotpath = joinpath(
                    survival_plotdir,
                    "ews_survival_$(ews_metric).$(output_format)",
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
                        plottitle = survival_plottitle,
                        subtitle = survival_subtitle,
                    )
                    save(
                        plotpath,
                        survival_plot,
                    )
                end
            end
        end
        next!(prog)
    end
    return nothing
end
