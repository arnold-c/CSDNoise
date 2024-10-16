using DataFrames
using StyledStrings
using ProgressMeter
using UnPack: @unpack

function ews_hyperparam_optimization(
    specification_vecs,
    data_arrs;
    filepath = outdir("ews_hyperparam_optimization.jld2"),
    io_file = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
    force = false,
    specification_vec_tuples = (
        noise_specification = NoiseSpecification[],
        test_specification = IndividualTestSpecification[],
        percent_tested = Float64[],
        ews_metric_specification = EWSMetricSpecification[],
        ews_enddate_type = EWSEndDateType[],
        ews_threshold_window = EWSThresholdWindow[],
        ews_threshold_burnin = Int[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_metric = String[],
    ),
)
    ews_df =
        if isfile(filepath) && !force
            load(filepath)["ews_df"]
        else
            DataFrame(
            (
                specification_vec_tuples...,
                true_positives = Int64[],
                true_negatives = Int64[],
                accuracy = Float64[],
                sensitivity = Float64[],
                specificity = Float64[],
            )
    )
        end

    return ews_hyperparam_optimization!(
        ews_df,
        specification_vecs,
        data_arrs;
        io_file = io_file,
        specification_vec_tuples = specification_vec_tuples,
    )
end

function ews_hyperparam_optimization!(
    ews_df,
    specification_vecs,
    data_arrs;
    io_file = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
    specification_vec_tuples = (
        noise_specification = NoiseSpecification[],
        test_specification = IndividualTestSpecification[],
        percent_tested = Float64[],
        ews_metric_specification = EWSMetricSpecification[],
        ews_enddate_type = EWSEndDateType[],
        ews_threshold_window = EWSThresholdWindow[],
        ews_threshold_burnin = Int[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_metric = String[],
    ),
)
    @assert map(propertynames(specification_vecs)) do pn
        Symbol(match(r"(.*)(_vec)$", string(pn)).captures[1])
    end ==
        propertynames(specification_vec_tuples)

    missing_specification_vecs = check_missing_ews_hyperparameter_simulations(
        ews_df,
        specification_vecs;
        specification_vec_tuples = specification_vec_tuples,
    )

    @assert map(propertynames(missing_specification_vecs)) do pn
        Symbol(match(r"(missing_)(.*)$", string(pn)).captures[2])
    end ==
        propertynames(specification_vecs)

    @unpack missing_noise_specification_vec,
    missing_test_specification_vec,
    missing_percent_tested_vec,
    missing_ews_metric_specification_vec,
    missing_ews_enddate_type_vec,
    missing_ews_threshold_window_vec,
    missing_ews_threshold_burnin_vec,
    missing_ews_threshold_percentile_vec,
    missing_ews_consecutive_thresholds_vec,
    missing_ews_metric_vec = missing_specification_vecs

    @unpack ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec,
    ensemble_single_periodsum_vecs =
        data_arrs

    ensemble_nsims = size(ensemble_single_incarr, 3)

    println("All worked as expected")
    println(missing_specification_vecs)
    return nothing

    io = open(io_file, "a")

    @showprogress for noise_specification in missing_noise_specification_vec
        println(
            styled"{green:\n=================================================================}"
        )
        println(
            styled"Noise type: {green,inverse: $(getdirpath(noise_specification))}"
        )

        noisearr = create_noise_arr(
            noise_specification,
            ensemble_single_incarr;
            ensemble_specification = ensemble_specification,
            seed = 1234,
        )[1]

        for (test_specification, percent_tested) in
            Iterators.product(
            missing_test_specification_vec, missing_percent_tested_vec
        )
            println(
                styled"\t\t\t\t-> Test specification: {blue: $(get_test_description(test_specification))}, Percent tested: {red,inverse: $(percent_tested)}"
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

            for (
                ews_metric_specification,
                ews_enddate_type,
                ews_threshold_window,
                ews_burnin,
                ews_percentile,
                ews_consecutive_thresholds,
            ) in
                Iterators.product(
                missing_ews_metric_specification_vec,
                missing_ews_enddate_type_vec,
                missing_ews_threshold_window_vec,
                missing_ews_threshold_burnin_vec,
                missing_ews_threshold_percentile_vec,
                missing_ews_consecutive_thresholds_vec,
            )
                ews_enddate_type_str = split(string(ews_enddate_type), "::")[1]
                println(
                    styled"\t\tEWS hyperparameters\n\t\tEWS metric specification: {blue,inverse: $(ews_metric_specification.dirpath)}, End date type: {magenta: $(ews_enddate_type_str)}, EWS burn-in: {yellow: $(ews_burnin)}, EWS percentile: {magenta,inverse: $(ews_percentile)}, EWS consecutive thresholds: {yellow,inverse: $(ews_consecutive_thresholds)}"
                )

                thresholds = SumTypes.@cases ews_enddate_type begin
                    [Reff_start, Reff_end] =>
                        ensemble_single_Reff_thresholds_vec
                    [Outbreak_start, Outbreak_end, Outbreak_middle] =>
                        ensemble_single_periodsum_vecs
                end

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
                    undef, size(testarr, 3), length(missing_ews_metric_vec)
                )
                null_exceeds_threshold_arr = Array{Matrix{Bool},2}(
                    undef, size(testarr, 3), length(missing_ews_metric_vec)
                )
                detection_index_arr = Array{Union{Nothing,Int64},2}(
                    undef, size(testarr, 3), length(missing_ews_metric_vec)
                )
                null_detection_index_arr = Array{Union{Nothing,Int64},2}(
                    undef, size(testarr, 3), length(missing_ews_metric_vec)
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

                        for (j, ews_metric) in pairs(missing_ews_metric_vec)
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

                for (j, ews_metric) in pairs(missing_ews_metric_vec)
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
                            percent_tested,
                            ews_metric_specification,
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
            end
        end
    end
end

function check_missing_ews_hyperparameter_simulations(
    ews_df,
    specification_vecs;
    specification_vec_tuples = (
        noise_specification = NoiseSpecification[],
        test_specification = IndividualTestSpecification[],
        percent_tested = Float64[],
        ews_metric_specification = EWSMetricSpecification[],
        ews_enddate_type = EWSEndDateType[],
        ews_threshold_window = EWSThresholdWindow[],
        ews_threshold_burnin = Int[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_metric = String[],
    ),
)
    run_params_df = DataFrame(
        specification_vec_tuples
    )

    for (
        noise_specification,
        test_specification,
        percent_tested,
        ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window,
        ews_threshold_burnin,
        ews_threshold_percentile,
        ews_consecutive_thresholds,
        ews_metric,
    ) in Iterators.product(specification_vecs...)
        push!(
            run_params_df,
            (
                noise_specification,
                test_specification,
                percent_tested,
                ews_metric_specification,
                ews_enddate_type,
                ews_threshold_window,
                ews_threshold_burnin,
                ews_threshold_percentile,
                ews_consecutive_thresholds,
                ews_metric,
            ),
        )
    end

    missing_run_params_df = antijoin(
        run_params_df,
        ews_df;
        on = [propertynames(specification_vec_tuples)...],
    )

    missing_params_nt = create_missing_run_params_nt(missing_run_params_df)

    return missing_params_nt
end

function create_missing_run_params_nt(missing_run_params_df)
    pns = names(missing_run_params_df)
    output_pns = Symbol.("missing_" .* pns .* "_vec")
    vals = map(pns) do pn
        unique(missing_run_params_df[!, pn])
    end

    return (; zip(output_pns, vals)...)
end
