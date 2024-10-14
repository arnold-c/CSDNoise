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
using Random: Random
using Distributions: Distributions
using StyledStrings

using ProgressMeter

include(srcdir("makie-plotting-setup.jl"))
include(srcdir("ensemble-parameters.jl"))

#%%
ensemble_model_type = ("seasonal-infectivity-import", "tau-leaping")

burnin_years = 5
nyears = 20
burnin_time = 365.0 * burnin_years
ensemble_time_specification = SimTimeParameters(;
    burnin = 365.0 * burnin_years, tmin = 0.0, tmax = 365.0 * nyears,
    tstep = 1.0,
)

ensemble_state_specification = StateParameters(
    500_000,
    Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
)

mu = calculate_mu(27)
beta_mean = calculate_beta(
    R0, GAMMA, mu, 1, ensemble_state_specification.init_states.N
)
epsilon = calculate_import_rate(
    mu, R0, ensemble_state_specification.init_states.N
)

min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
    0.9,
    ensemble_state_specification.init_states.S,
    ensemble_state_specification.init_states.N,
    R0,
    mu,
    burnin_years * 2,
)

max_burnin_vaccination_coverage = 1.0
min_vaccination_coverage = 0.0
max_vaccination_coverage = 0.8

ensemble_dynamics_specification = DynamicsParameterSpecification(
    beta_mean,
    0.0,
    cos,
    SIGMA,
    GAMMA,
    mu,
    27,
    epsilon,
    R0,
    min_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage,
    min_vaccination_coverage,
    max_vaccination_coverage,
)

#%%
ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
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
ensemble_vax_plotpath = joinpath(
    plotsdir(),
    "ensemble",
    "burnin-time-$(burnin_time)",
    "min-burnin-vax_$(min_burnin_vaccination_coverage)",
    "max-burnin-vax_$(max_burnin_vaccination_coverage)",
    "min-vax_$(min_vaccination_coverage)",
    "max-vax_$(max_vaccination_coverage)",
)
mkpath(ensemble_vax_plotpath)

#%%
nsims_plot = 10
Random.seed!(1234)
selected_sims = StatsBase.sample(
    collect(1:size(ensemble_single_incarr, 3)),
    nsims_plot;
    replace = false,
)

ensemble_single_scenario_incidence_Reff_plot = ensemble_incarr_Reff_plot(
    ensemble_single_incarr[:, :, selected_sims],
    ensemble_single_Reff_arr[:, selected_sims],
    ensemble_single_Reff_thresholds_vec[selected_sims],
    ensemble_single_periodsum_vecs[selected_sims], ;
    outbreak_alpha = 0.1,
    Reff_alpha = 1,
)

plotpath = joinpath(
    ensemble_vax_plotpath,
    "ensemble-sim_single-scenario_incidence_Reff.png",
)

save(
    plotpath,
    ensemble_single_scenario_incidence_Reff_plot;
    size = (2200, 1600),
)

#%%
calculate_vaccination_rate_to_achieve_Reff(
    0.9, ensemble_state_specification.init_states,
    ensemble_dynamics_specification, 10,
)

#%%
ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_incarr,
    ensemble_single_seir_arr,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    sim = 1,
    threshold = ensemble_outbreak_specification.outbreak_threshold,
)

plotpath = joinpath(
    ensemble_vax_plotpath,
    "ensemble-sim_single-scenario_incidence_prevalence.png",
)

save(
    plotpath,
    ensemble_single_scenario_incidence_prevalence_plot;
    size = (2200, 1600),
)

#%%
# Open a textfile for writing
io = open(scriptsdir("ensemble-sim_ews-visualization.log.txt"), "a")

force = false

test_specification_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

percent_tested_vec = [1.0]

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

ews_threshold_window = Main.Expanding
ews_threshold_percentile = 0.95
consecutive_thresholds = 2

ews_df = DataFrame(
    "ews_metric" => String[],
    "test_specification" => IndividualTestSpecification[],
    "ews_metric_specification" => EWSMetricSpecification[],
    "ews_enddate_type" => EWSEndDateType[],
    "ews_metric_value" => Float64[],
    "ews_metric_vector" => Vector{Float64}[],
)

@showprogress for (
    noise_specification, percent_tested, ews_metric_specification
) in
                  Iterators.product(
    [ensemble_noise_specification_vec[1]],
    percent_tested_vec,
    ews_spec_vec,
)
    println(
        styled"{green:\n=================================================================}"
    )
    println(
        styled"Creating EWS visualization for:\n{green,inverse: $(getdirpath(noise_specification))}, {red,inverse: $(percent_tested)}, {blue,inverse: $(ews_metric_specification.dirpath)}"
    )

    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )[1]
    noisedir = getdirpath(noise_specification)

    for ews_enddate_type in
        (
        Main.Reff_start,
        Main.Reff_end,
        Main.Outbreak_start,
        # Main.Outbreak_end,
        # Main.Outbreak_middle,
    )
        ews_enddate_type_str = split(string(ews_enddate_type), "::")[1]
        println(
            styled"\t\tEWS end date type: {magenta: $(ews_enddate_type_str)}"
        )

        thresholds = SumTypes.@cases ews_enddate_type begin
            [Reff_start, Reff_end] =>
                ensemble_single_Reff_thresholds_vec
            [Outbreak_start, Outbreak_end, Outbreak_middle] =>
                ensemble_single_periodsum_vecs
        end

        for test_specification in test_specification_vec
            println(
                styled"\t\t\t\t-> Test specification tau distribution: {blue: $(get_test_description(test_specification))}"
            )

            testarr = create_testing_arrs(
                ensemble_single_incarr,
                noisearr,
                percent_tested,
                test_specification,
            )
            enddate_vec = zeros(Int64, size(testarr, 3))
            failed_sims = zeros(Int64, size(testarr, 3))
            ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
                undef, size(testarr, 3)
            )
            inc_ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
                undef, size(testarr, 3)
            )
            fill!(ews_vals_vec, missing)
            fill!(inc_ews_vals_vec, missing)

            exceeds_threshold_arr = Array{Matrix{Bool},2}(
                undef, size(testarr, 3), length(ews_metrics)
            )
            detection_index_arr = Array{Union{Nothing,Int64},2}(
                undef, size(testarr, 3), length(ews_metrics)
            )
            fill!(detection_index_arr, nothing)

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

                    for (j, ews_metric) in pairs(ews_metrics)
                        exceeds_threshold_arr[sim, j] = expanding_ews_thresholds(
                            ews_vals_vec[sim],
                            Symbol(ews_metric),
                            ews_threshold_window;
                            percentiles = ews_threshold_percentile,
                        )[2]

                        detection_index_arr[sim, j] = calculate_ews_trigger_index(
                            exceeds_threshold_arr[sim, j];
                            consecutive_thresholds = consecutive_thresholds,
                        )
                    end

                else
                    failed_sims[sim] = sim
                    write(
                        io,
                        "$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))\n",
                    )
                    write(io, "Error:\t$(Try.unwrap_err(enddate))\n")
                    write(io, "Simulation:\t$(sim)\n\n")
                end
            end

            filtered_ews_vals_vec = filter(!ismissing, ews_vals_vec)
            filtered_inc_ews_vals_vec = filter(!ismissing, inc_ews_vals_vec)
            filter!(x -> x != 0, failed_sims)

            # filter!(!isnan, ews_lead_time)
            # percentile_tail = (1 - lead_time_percentile) / 2

            @assert length(ews_vals_vec) == length(inc_ews_vals_vec)

            ews_vals_sa = StructArray(
                convert(Vector{EWSMetrics}, filtered_ews_vals_vec)
            )
            inc_ews_vals_sa = StructArray(
                convert(Vector{EWSMetrics}, filtered_inc_ews_vals_vec)
            )

            ensemble_noise_plotpath = joinpath(
                ensemble_vax_plotpath,
                noisedir,
                "percent-tested_$(percent_tested)",
                "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)",
            )
            ensemble_ews_plotpath = joinpath(
                ensemble_noise_plotpath,
                ews_metric_specification.dirpath,
                ews_enddate_type_str,
            )
            mkpath(ensemble_ews_plotpath)

            for ews_metric in ews_metrics
                plotdir = joinpath(
                    ensemble_ews_plotpath,
                    "tau-distributions",
                )
                mkpath(plotdir)
                plotpath = joinpath(
                    plotdir,
                    "ensemble-sim_single-scenario_ews-$(ews_metric)-tau-distribution.png",
                )

                if !isfile(plotpath) || force
                    plot = simulation_tau_distribution(
                        ews_vals_sa,
                        inc_ews_vals_sa,
                        ews_metric;
                        plottitle = "$(get_test_description(test_specification)) ($(percent_tested*100)% tested), $(min_vaccination_coverage)-$(max_vaccination_coverage) Vaccination Coverage, $(get_noise_magnitude_description(noise_specification)): $(method_string(ews_metric_specification.method)) $(ews_enddate_type_str) EWS $(ews_metric) Tau Distribution",
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
                    ews_metric_specification = ews_metric_specification,
                    ews_enddate_type = ews_enddate_type,
                    statistic_function = StatsBase.mean,
                )
            end

            println(
                styled"\t\t\t\t-> Single scenario plots"
            )

            for sim in [selected_sims[1]]
                if sim in failed_sims
                    write(
                        io,
                        "Tried to plot single EWS for simulation $(sim), but failed as no end date was found\n\n",
                    )
                    continue
                end
                enddate = enddate_vec[sim]

                ews_vals = ews_vals_vec[sim]

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
                    7,
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
                        (ensemble_single_periodsum_vecs[sim][:, 4] .== 1),
                        [1, 2],
                    ] .÷ ews_metric_specification.aggregation

                plotdir = joinpath(
                    ensemble_ews_plotpath, "single-scenario", "sim-$(sim)"
                )
                mkpath(plotdir)

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
                    exceeds_threshold_arr[sim, :],
                    detection_index_arr[sim, :],
                    ews_metric_specification.dirpath,
                    ews_enddate_type_str,
                    test_specification,
                    percent_tested,
                    ensemble_time_specification;
                    aggregation = ews_metric_specification.aggregation,
                    sim = sim,
                    force = true,
                    base_plotpath = plotdir,
                )
            end
        end

        plotpath = joinpath(
            plotsdir(),
            "ensemble",
            "burnin-time-$(burnin_time)",
            "min-burnin-vax_$(min_burnin_vaccination_coverage)",
            "max-burnin-vax_$(max_burnin_vaccination_coverage)",
            "min-vax_$(min_vaccination_coverage)",
            "max-vax_$(max_vaccination_coverage)",
            noisedir,
            "percent-tested_$(percent_tested)",
            "tau-heatmaps",
            ews_metric_specification.dirpath,
            ews_enddate_type_str,
        )
        mkpath(plotpath)
        plotpath = joinpath(
            plotpath, "ews-tau-heatmap_mean.png"
        )

        println(
            styled"\t\t\t\t-> Tau heatmap"
        )
        tau_heatmap = tycho_tau_heatmap_plot(
            subset(
                ews_df,
                :ews_enddate_type => ByRow(==(ews_enddate_type)),
                :ews_metric_specification =>
                    ByRow(==(ews_metric_specification)),
            );
            statistic_function = titlecase("mean"),
            plottitle = "Kendall's Tau Heatmap (Mean)\n($(percent_tested*100)% tested), $(min_vaccination_coverage)-$(max_vaccination_coverage) Vaccination Coverage, $(ews_metric_specification.dirpath), $(ews_enddate_type_str), $(get_noise_magnitude_description(noise_specification))",
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
