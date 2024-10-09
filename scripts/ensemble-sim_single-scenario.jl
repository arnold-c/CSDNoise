#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes
using Dates

using CSDNoise
using StructArrays
using SumTypes
using Try
using DataFrames
using StatsBase: StatsBase

using ProgressMeter

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec = IndividualTestSpecification(
    1.0, 1.0, 0
)
ensemble_single_outbreak_detection_spec = OutbreakDetectionSpecification(
    10, 7, 1.0, 1.0, "movingavg"
)

#%%
ensemble_single_seir_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_seir_arr"]

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

ensemble_single_Reff_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_arr"]

ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

#%%
nsims_plot = 5
ensemble_incarr_Reff_plot(
    ensemble_single_incarr[:, :, 1:nsims_plot],
    ensemble_single_Reff_arr[:, 1:nsims_plot],
    ensemble_single_Reff_thresholds_vec[1:nsims_plot],
    ensemble_single_periodsum_vecs[1:nsims_plot], ;
    outbreak_alpha = 0.1,
    Reff_alpha = 1,
)

#%%
ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_incarr,
    ensemble_single_seir_arr,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    threshold = ensemble_outbreak_specification.outbreak_threshold,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_incidence_prevalence.png"
    ),
    ensemble_single_scenario_incidence_prevalence_plot,
)

#%%
# Open a textfile for writing
io = open(scriptsdir("ensemble-sim_single-scenario.log.txt"), "a")
write(io, "============================================\n")

force = false

test_specification_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

sims = (
    1,
    4,
)

ews_method_vec = [
    Centered,
    Backward,
]
ews_aggregation_days_vec = [7, 14, 28]
ews_bandwidth_days_vec = [52 * 7]
ews_lag_days_vec = [1]

ews_spec_vec = create_combinations_vec(
    calculate_bandwidth_and_return_ews_metric_spec,
    (
        ews_method_vec,
        ews_aggregation_days_vec,
        ews_bandwidth_days_vec,
        ews_lag_days_vec,
    ),
)

ews_metrics = [
    "autocorrelation",
    "autocovariance",
    "coefficient_of_variation",
    "index_of_dispersion",
    "kurtosis",
    "mean",
    "skewness",
    "variance",
]

ews_df = DataFrame(
    "ews_metric" => String[],
    "test_specification" => IndividualTestSpecification[],
    "ews_enddate_type" => EWSEndDateType[],
    "ews_metric_value" => Float64[],
    "ews_metric_vector" => Vector{Float64}[],
)

@showprogress for (noise_specification, ews_metric_specification) in
                  Iterators.product(
    [ensemble_noise_specification_vec[1]],
    [ews_spec_vec[1]],
)
    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )[1]
    noisedir = getdirpath(noise_specification)

    for test_specification in test_specification_vec
        testarr = create_testing_arrs(
            ensemble_single_incarr,
            noisearr,
            ensemble_single_outbreak_detection_spec.percent_tested,
            test_specification,
        )

        for ews_enddate_type in
            (
            Main.Reff_start,
            Main.Reff_end,
            Main.Outbreak_start,
            Main.Outbreak_end,
            Main.Outbreak_middle,
        )
            thresholds = SumTypes.@cases ews_enddate_type begin
                [Reff_start, Reff_end] =>
                    ensemble_single_Reff_thresholds_vec
                [Outbreak_start, Outbreak_end, Outbreak_middle] =>
                    ensemble_single_periodsum_vecs
            end

            enddate_vec = zeros(Int64, size(testarr, 3))
            ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
                undef, size(testarr, 3)
            )
            inc_ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
                undef, size(testarr, 3)
            )
            fill!(ews_vals_vec, missing)
            fill!(inc_ews_vals_vec, missing)

            for sim in axes(testarr, 3)
                enddate = calculate_ews_enddate(
                    thresholds[sim],
                    ews_enddate_type,
                )

                if Try.isok(enddate)
                    enddate_vec[sim] = Try.unwrap(enddate)

                    inc_ews_vals_vec[sim] = EWSMetrics(
                        ews_metric_specification,
                        @view(
                            ensemble_single_incarr[1:enddate_vec[sim], 1, sim]
                        )
                    )

                    ews_vals_vec[sim] = EWSMetrics(
                        ews_metric_specification,
                        @view(testarr[1:enddate_vec[sim], 5, sim])
                    )
                else
                    write(
                        io,
                        "$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))\n",
                    )
                    write(io, "Error:\t$(Try.unwrap_err(enddate))\n")
                    write(io, "Simulation:\t$(sim)\n\n")
                end
            end

            filter!(x -> !ismissing(x), ews_vals_vec)
            filter!(x -> !ismissing(x), inc_ews_vals_vec)

            @assert length(ews_vals_vec) == length(inc_ews_vals_vec)

            ews_vals_sa = StructArray(convert(Vector{EWSMetrics}, ews_vals_vec))
            inc_ews_vals_sa = StructArray(
                convert(Vector{EWSMetrics}, inc_ews_vals_vec)
            )

            ensemble_noise_plotpath = joinpath(
                plotsdir(),
                "ensemble",
                noisedir,
                "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)",
                ews_metric_specification.dirpath,
                split(string(ews_enddate_type), "::")[1],
            )
            mkpath(ensemble_noise_plotpath)

            for ews_metric in ews_metrics
                plotpath = joinpath(
                    ensemble_noise_plotpath,
                    "ensemble-sim_single-scenario_ews-$(ews_metric)-tau-distribution.png",
                )

                if !isfile(plotpath) || force
                    plot = simulation_tau_distribution(
                        ews_vals_sa,
                        inc_ews_vals_sa,
                        ews_metric;
                        plottitle = "$(get_test_description(test_specification)), $(get_noise_magnitude_description(noise_specification)): $(method_string(ews_metric_specification.method)) $(split(string(ews_enddate_type), "::")[1]) EWS $(ews_metric) Tau Distribution",
                    )

                    save(
                        plotpath,
                        plot;
                        size = (2200, 1600),
                    )
                end

                simulation_tau_heatmap_df!(
                    ews_df,
                    ews_vals_sa,
                    ews_metric;
                    individual_test_specification = test_specification,
                    ews_enddate_type = ews_enddate_type,
                    statistic_function = StatsBase.mean,
                )
            end

            GC.gc(true)
            @info "Finished plotting the single scenario for $(noisedir), $(ews_metric_specification.dirpath), $(ews_enddate_type), $(get_test_description(test_specification))"
            println(
                "================================================================="
            )
        end
    end

    for ews_enddate_type in (
        Main.Reff_start,
        Main.Reff_end,
        Main.Outbreak_start,
        Main.Outbreak_end,
        Main.Outbreak_middle,
    )
        plotpath = joinpath(
            plotsdir(),
            "ensemble",
            noisedir,
            "tau-heatmaps",
            ews_metric_specification.dirpath,
            split(string(ews_enddate_type), "::")[1],
        )
        mkpath(plotpath)
        plotpath = joinpath(
            plotpath, "ews-tau-heatmap_mean.png"
        )

        tau_heatmap = tycho_tau_heatmap_plot(
            subset(
                ews_df, :ews_enddate_type => ByRow(==(ews_enddate_type))
            );
            statistic_function = titlecase("mean"),
            plottitle = "Kendall's Tau Heatmap (Mean)\n$(method_string(ews_metric_specification.method)), $(split(string(ews_enddate_type), "::")[1]), $(get_noise_magnitude_description(noise_specification))",
        )

        save(
            plotpath,
            tau_heatmap;
            size = (2200, 1600),
        )
    end
end

close(io)

#%%
#     for sim in sims
#         enddate = enddate_vec[sim]
#
#         ews_vals = ews_vals_vec[sim]
#
#         aggregated_noise_vec = aggregate_timeseries(
#             @view(noisearr[:, sim]),
#             ews_metric_specification.aggregation,
#         )
#
#         aggregated_inc_vec = aggregate_timeseries(
#             @view(ensemble_single_incarr[:, 1, sim]),
#             ews_metric_specification.aggregation,
#         )
#         aggregated_outbreak_status_vec = aggregate_thresholds_vec(
#             @view(ensemble_single_incarr[:, 3, sim]),
#             ews_metric_specification.aggregation,
#         )
#
#         aggregated_test_vec = aggregate_timeseries(
#             @view(testarr[:, 5, sim]),
#             ews_metric_specification.aggregation,
#         )
#
#         aggregated_test_movingavg_vec = zeros(
#             Int64, size(aggregated_test_vec)
#         )
#
#         calculate_movingavg!(
#             aggregated_test_movingavg_vec,
#             aggregated_test_vec,
#             ensemble_single_outbreak_detection_spec.moving_average_lag,
#         )
#
#         aggregated_Reff_vec = aggregate_Reff_vec(
#             @view(ensemble_single_Reff_arr[:, sim]),
#             ews_metric_specification.aggregation,
#         )
#
#         aggregated_Reff_thresholds_arr =
#             ensemble_single_Reff_thresholds_vec[sim] .÷
#             ews_metric_specification.aggregation
#
#         aggregated_outbreak_thresholds_arr =
#             ensemble_single_periodsum_vecs[sim][
#                 (ensemble_single_periodsum_vecs[sim][:, 4] .== 1), [1, 2]
#             ] .÷ ews_metric_specification.aggregation
#
#         plot_all_single_scenarios(
#             aggregated_noise_vec,
#             noisedir,
#             aggregated_inc_vec,
#             aggregated_outbreak_status_vec,
#             aggregated_test_vec,
#             aggregated_test_movingavg_vec,
#             aggregated_Reff_vec,
#             aggregated_Reff_thresholds_arr,
#             aggregated_outbreak_thresholds_arr,
#             ews_vals,
#             ews_metric_specification.dirpath,
#             split(string(ews_enddate_type), "::")[1],
#             ensemble_single_individual_test_spec,
#             ensemble_single_outbreak_detection_spec,
#             ensemble_time_specification;
#             aggregation = ews_metric_specification.aggregation,
#             sim = sim,
#             force = true,
#         )
#     end
#
#     GC.gc(true)
#     @info "Finished plotting the single scenario for $(noisedir), $(ews_metric_specification.dirpath)"
#     println("=================================================================")
# end
