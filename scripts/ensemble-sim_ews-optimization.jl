using DataFrames: IteratorSize
#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack

using CSDNoise
using StructArrays
using SumTypes
using Try
using DataFrames
using StatsBase: StatsBase
using Random: Random
using Distributions: Distributions
using StyledStrings
using Dates

using ProgressMeter

include(srcdir("makie-plotting-setup.jl"))

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
    burnin_years * 2,
    ensemble_state_specification.init_states.S,
    ensemble_state_specification.init_states.N,
    R0,
    mu,
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
    0.6,
    0.8,
)

null_dynamics_specification = DynamicsParameterSpecification(
    map(
        pn -> getproperty(ensemble_dynamics_specification, pn),
        filter(
            name ->
                name != :min_vaccination_coverage &&
                    name != :max_vaccination_coverage,
            propertynames(ensemble_dynamics_specification),
        ),
    )...,
    nothing,
    nothing,
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

null_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    null_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

ensemble_outbreak_specification = OutbreakSpecification(
    5, 30, 500
)

#%%
ensemble_dynamics_parameters_sa = StructVector(
    get_ensemble_file(
        ensemble_specification
    )["dynamics_parameters"],
)

null_dynamics_parameters_sa = StructVector(
    get_ensemble_file(
        null_specification
    )["dynamics_parameters"],
)

@assert ensemble_dynamics_parameters_sa.burnin_vaccination_coverage ==
    null_dynamics_parameters_sa.burnin_vaccination_coverage ==
    null_dynamics_parameters_sa.vaccination_coverage

#%%
ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

null_single_scenario_inc_file = get_ensemble_file(
    null_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

null_single_incarr = null_single_scenario_inc_file["ensemble_inc_arr"]
null_single_periodsum_vecs = null_single_scenario_inc_file["ensemble_thresholds_vec"]

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
# Open a textfile for writing
io = open(scriptsdir("ensemble-sim_ews-optimization.log.txt"), "a")

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
ews_threshold_burnin = 10
ews_threshold_percentile = 0.95
consecutive_thresholds = 2

ews_df = DataFrame(
    "test_specification" => IndividualTestSpecification[],
    "ews_enddate_type" => EWSEndDateType[],
    "noise_specification" => NoiseSpecification[],
    "ews_metric" => String[],
    "ews_threshold_percentile" => Float64[],
    "consecutive_thresholds" => Int[],
    "ews_threshold_burnin" => Int[],
    "true_positives" => Int64[],
    "true_negatives" => Int64[],
    "accuracy" => Float64[],
    "sensitivity" => Float64[],
    "specificity" => Float64[],
)

@showprogress for (
    noise_specification, percent_tested, ews_metric_specification
) in
                  Iterators.product(
    [PoissonNoiseSpecification(8.0)],
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

    for (
        ews_enddate_type, ews_burnin, ews_percentile, ews_consecutive_thresholds
    ) in
        Iterators.product(
        (
            Main.Reff_start,
            Main.Reff_end,
        ),
        collect(10:1:15),
        collect(0.9:0.02:0.98),
        collect(2:1:10),
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

            null_testarr = create_testing_arrs(
                null_single_incarr,
                noisearr,
                percent_tested,
                test_specification,
            )

            enddate_vec = zeros(Int64, size(testarr, 3))
            failed_sims = zeros(Int64, size(testarr, 3))
            ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
                undef, size(testarr, 3)
            )
            null_ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
                undef, size(testarr, 3)
            )
            fill!(ews_vals_vec, missing)
            fill!(null_ews_vals_vec, missing)

            exceeds_threshold_arr = Array{Matrix{Bool},2}(
                undef, size(testarr, 3), length(ews_metrics)
            )
            null_exceeds_threshold_arr = Array{Matrix{Bool},2}(
                undef, size(testarr, 3), length(ews_metrics)
            )
            detection_index_arr = Array{Union{Nothing,Int64},2}(
                undef, size(testarr, 3), length(ews_metrics)
            )
            null_detection_index_arr = Array{Union{Nothing,Int64},2}(
                undef, size(testarr, 3), length(ews_metrics)
            )
            fill!(detection_index_arr, nothing)
            fill!(null_detection_index_arr, nothing)

            for sim in axes(testarr, 3)
                enddate = calculate_ews_enddate(
                    thresholds[sim],
                    ews_enddate_type,
                )

                if Try.isok(enddate)
                    enddate_vec[sim] = Try.unwrap(enddate)

                    ews_vals_vec[sim] = EWSMetrics(
                        ews_metric_specification,
                        @view(testarr[1:enddate_vec[sim], 5, sim])
                    )

                    null_ews_vals_vec[sim] = EWSMetrics(
                        ews_metric_specification,
                        @view(null_testarr[1:enddate_vec[sim], 5, sim])
                    )

                    for (j, ews_metric) in pairs(ews_metrics)
                        exceeds_threshold_arr[sim, j] = expanding_ews_thresholds(
                            ews_vals_vec[sim],
                            Symbol(ews_metric),
                            ews_threshold_window;
                            percentiles = ews_percentile,
                            burn_in = ews_burnin,
                        )[2]

                        detection_index_arr[sim, j] = calculate_ews_trigger_index(
                            exceeds_threshold_arr[sim, j];
                            consecutive_thresholds = ews_consecutive_thresholds,
                        )

                        null_exceeds_threshold_arr[sim, j] = expanding_ews_thresholds(
                            null_ews_vals_vec[sim],
                            Symbol(ews_metric),
                            ews_threshold_window;
                            percentiles = ews_percentile,
                            burn_in = ews_burnin,
                        )[2]

                        null_detection_index_arr[sim, j] = calculate_ews_trigger_index(
                            null_exceeds_threshold_arr[sim, j];
                            consecutive_thresholds = ews_consecutive_thresholds,
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

            for (j, ews_metric) in pairs(ews_metrics)
                true_positives = length(
                    filter(!ismissing, detection_index_arr[:, j])
                )
                true_negatives = length(
                    filter(ismissing, null_detection_index_arr[:, j])
                )
                sensitivity = true_positives / ensemble_nsims
                specificity = true_negatives / ensemble_nsims
                accuracy = (sensitivity + specificity) / 2

                push!(
                    ews_df,
                    (
                        test_specification,
                        ews_enddate_type,
                        noise_specification,
                        ews_metric,
                        ews_percentile,
                        ews_consecutive_thresholds,
                        ews_burnin,
                        true_positives,
                        true_negatives,
                        accuracy,
                        sensitivity,
                        specificity,
                    ),
                )
            end

            # filtered_ews_vals_vec = filter(!ismissing, ews_vals_vec)
            # filtered_null_ews_vals_vec = filter(!ismissing, null_ews_vals_vec)
            # filter!(x -> x != 0, failed_sims)
        end
    end
end

close(io)
