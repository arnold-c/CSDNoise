#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes

using CSDNoise
using StructArrays
using SumTypes

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

for (noise_specification, ews_metric_specification) in Iterators.product(
    ensemble_noise_specification_vec[1:2],
    ews_spec_vec,
)
    scenario_specification = ScenarioSpecification(
        ensemble_specification,
        ensemble_outbreak_specification,
        noise_specification,
        ensemble_single_outbreak_detection_spec,
        ensemble_single_individual_test_spec,
        ews_metric_specification,
    )

    noisearr, poisson_noise_prop = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )
    noisedir = getdirpath(noise_specification)

    testarr = create_testing_arrs(
        ensemble_single_incarr,
        noisearr,
        ensemble_single_outbreak_detection_spec.percent_tested,
        ensemble_single_individual_test_spec,
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

        for sim in sims
            enddate = calculate_ews_enddate(
                thresholds[sim],
                ews_enddate_type,
            )

            Debugger.@bp

            ews_vals = EWSMetrics(
                ews_metric_specification,
                @view(testarr[1:enddate, 5, sim]),
            )

            aggregated_noise_vec = aggregate_timeseries(
                @view(noisearr[:, sim]),
                ews_metric_specification.aggregation,
            )

            aggregated_inc_vec = aggregate_timeseries(
                @view(ensemble_single_incarr[:, 1, sim]),
                ews_metric_specification.aggregation,
            )
            aggregated_outbreak_status_vec = aggregate_thresholds_vec(
                @view(ensemble_single_incarr[:, 3, sim]),
                ews_metric_specification.aggregation,
            )

            aggregated_test_vec = aggregate_timeseries(
                @view(testarr[:, 5, sim]),
                ews_metric_specification.aggregation,
            )

            aggregated_test_movingavg_vec = zeros(
                Int64, size(aggregated_test_vec)
            )
            calculate_movingavg!(
                aggregated_test_movingavg_vec,
                aggregated_test_vec,
                ensemble_single_outbreak_detection_spec.moving_average_lag,
            )

            aggregated_Reff_vec = aggregate_Reff_vec(
                @view(ensemble_single_Reff_arr[:, sim]),
                ews_metric_specification.aggregation,
            )

            aggregated_Reff_thresholds_arr =
                ensemble_single_Reff_thresholds_vec[sim] .÷
                ews_metric_specification.aggregation

            aggregated_outbreak_thresholds_arr =
                ensemble_single_periodsum_vecs[sim][
                    (ensemble_single_periodsum_vecs[sim][:, 4] .== 1), [1, 2]
                ] .÷ ews_metric_specification.aggregation

            plot_all_single_scenarios(
                aggregated_noise_vec,
                noisedir,
                aggregated_inc_vec,
                aggregated_outbreak_status_vec,
                aggregated_test_vec,
                aggregated_test_movingavg_vec,
                aggregated_Reff_vec,
                aggregated_Reff_thresholds_arr,
                aggregated_outbreak_thresholds_arr,
                ews_vals,
                ews_metric_specification.dirpath,
                split(string(ews_enddate_type), "::")[1],
                ensemble_single_individual_test_spec,
                ensemble_single_outbreak_detection_spec,
                ensemble_time_specification;
                aggregation = ews_metric_specification.aggregation,
                sim = sim,
                force = true,
            )
        end
    end

    GC.gc(true)
    @info "Finished plotting the single scenario for $(noisedir), $(ews_metric_specification.dirpath)"
    println("=================================================================")
end
