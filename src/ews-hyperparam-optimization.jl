using Dates: Dates
using DataFrames
using StyledStrings
using ProgressMeter
using UnPack: @unpack
using REPL: REPL
using REPL.TerminalMenus: RadioMenu, request
using Try: Try

function ews_hyperparam_optimization(
    specification_vecs,
    data_arrs;
    filedir = outdir("ensemble", "ews-hyperparam-optimization"),
    gridsearch_filename_base = "ews-hyperparam-gridsearch.jld2",
    gridsearch_output_filepath = joinpath(
        filedir,
        string(Dates.now()) * "_" * gridsearch_filename_base,
    ),
    optimization_filename_base = "ews-hyperparam-optimization.jld2",
    optimization_output_filepath = joinpath(
        filedir,
        string(Dates.now()) * "_" * optimization_filename_base,
    ),
    logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
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
    time_per_run_s = 0.11,
)
    ews_df = ews_hyperparam_gridsearch(
        specification_vecs,
        data_arrs;
        filedir = filedir,
        filename_base = gridsearch_filename_base,
        output_filepath = gridsearch_output_filepath,
        logfilepath = logfilepath,
        force = force,
        specification_vec_tuples = specification_vec_tuples,
        time_per_run_s = time_per_run_s,
        return_df = true,
    )

    load_filepath = load_most_recent_hyperparam_gridsearch_file(
        optimization_filename_base,
        filedir,
    )
    optimal_ews_df = filter_optimal_ews_hyperparam_gridsearch(ews_df)

    if Try.isok(load_filepath) && !force
        previous_optimal_ews_df = load(Try.unwrap(load_filepath))["optimal_ews_df"]

        if optimal_ews_df == previous_optimal_ews_df
            @info "🟢 Previous optimal ews_df is the same as the current. No need to save a new version. 🟢"
            return optimal_ews_df
        end
    end

    @tagsave(
        optimization_output_filepath,
        Dict("optimal_ews_df" => optimal_ews_df)
    )
    @info "🟢 Saved optimal ews_df to $(optimization_output_filepath) 🟢"

    return optimal_ews_df
end

function filter_optimal_ews_hyperparam_gridsearch(ews_df)
    return map(
        collect(
            groupby(
                ews_df,
                [
                    :noise_specification,
                    :test_specification,
                    :percent_tested,
                    :ews_metric_specification,
                    :ews_enddate_type,
                    :ews_metric,
                ],
            ),
        ),
    ) do df
        max_accuracy = maximum(df[!, :accuracy])
        return subset(df, :accuracy => ByRow(==(max_accuracy)))
    end |>
           x -> vcat(x...; cols = :union)
end

function ews_hyperparam_gridsearch(
    specification_vecs,
    data_arrs;
    filedir = outdir("ensemble", "ews-hyperparam-optimization"),
    filename_base = "ews-hyperparam-gridsearch.jld2",
    output_filepath = joinpath(
        filedir,
        string(Dates.now()) * "_" * filename_base,
    ),
    logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
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
    time_per_run_s = 0.11,
    return_df = true,
)
    load_filepath = load_most_recent_hyperparam_gridsearch_file(
        filename_base,
        filedir,
    )

    if Try.isok(load_filepath) && !force
        ews_df = load(Try.unwrap(load_filepath))["ews_df"]
    else
        ews_df = DataFrame(
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

    val = ews_hyperparam_gridsearch!(
        ews_df,
        specification_vecs,
        data_arrs;
        logfilepath = logfilepath,
        specification_vec_tuples = specification_vec_tuples,
        time_per_run_s = time_per_run_s,
    )

    if Try.iserr(val)
        println(Try.unwrap_err(val))
        if return_df
            println("Returning previously computed ews_df")
            return ews_df
        end

        return nothing
    end

    time_taken_s, missing_runs = Try.unwrap(val)
    time_per_run_s = time_taken_s / missing_runs

    time_taken_message = if time_taken_s < 60
        "$(round(time_taken_s; digits = 0)) seconds"
    else
        "$(round(time_taken_s / 60; digits = 2)) minutes"
    end

    @info "🟢 ews_hyperparam_optimization.jl completed $(missing_runs) missing simulations in $(time_taken_message) ($(round(time_per_run_s; digits = 4)) seconds per run). 🟢"

    @tagsave(output_filepath, Dict("ews_df" => ews_df))

    if return_df
        return ews_df
    end

    return nothing
end

function load_most_recent_hyperparam_gridsearch_file(
    filename_base,
    filedir,
)
    @assert isdir(filedir)
    optimization_files = readdir(filedir)

    if length(optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filter_regex = Regex("(.*)$(filename_base)\$")

    filtered_optimization_files = filter(
        f -> contains(f, filter_regex),
        optimization_files,
    )

    if length(filtered_optimization_files) == 0
        return Try.Err("No optimization files found.")
    end

    filtered_optimization_datetimes = Vector{Union{Try.Ok,Try.Err}}(
        undef, length(filtered_optimization_files)
    )

    for (i, f) in pairs(filtered_optimization_files)
        matches = match(filter_regex, f)
        if isnothing(matches)
            filtered_optimization_datetimes[i] = Try.Err(
                "No matches for filename $(f)"
            )
            continue
        end

        filtered_optimization_datetimes[i] = Try.Ok(
            tryparse(
                Dates.DateTime,
                strip(
                    matches[1],
                    '_',
                ),
            ),
        )
    end

    filtered_optimization_datetimes = filter(
        Try.isok, filtered_optimization_datetimes
    )

    if length(filtered_optimization_datetimes) == 0
        return Try.Err("No optimization files found.")
    end

    most_recent_optimization_datetime = sort(
        Try.unwrap.(filtered_optimization_datetimes)
    )[end]

    most_recent_filepath = joinpath(
        filedir,
        string(most_recent_optimization_datetime) *
        "_$(filename_base)",
    )
    return Try.Ok(most_recent_filepath)
end

function ews_hyperparam_gridsearch!(
    ews_df,
    specification_vecs,
    data_arrs;
    logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
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
    time_per_run_s = 0.11,
)
    @assert map(propertynames(specification_vecs)) do pn
        Symbol(match(r"(.*)(_vec)$", string(pn)).captures[1])
    end ==
        propertynames(specification_vec_tuples)

    missing_specification_vecs = check_missing_ews_hyperparameter_simulations(
        ews_df,
        specification_vecs;
        specification_vec_tuples = specification_vec_tuples,
        time_per_run_s = time_per_run_s,
    )

    if Try.iserr(missing_specification_vecs)
        return missing_specification_vecs
    end

    @assert map(propertynames(Try.unwrap(missing_specification_vecs))) do pn
        Symbol(match(r"(missing_)(.*)$", string(pn)).captures[2])
    end ==
        (propertynames(specification_vecs)..., :runs)

    @unpack missing_noise_specification_vec,
    missing_test_specification_vec,
    missing_percent_tested_vec,
    missing_ews_metric_specification_vec,
    missing_ews_enddate_type_vec,
    missing_ews_threshold_window_vec,
    missing_ews_threshold_burnin_vec,
    missing_ews_threshold_percentile_vec,
    missing_ews_consecutive_thresholds_vec,
    missing_ews_metric_vec,
    missing_runs = Try.unwrap(missing_specification_vecs)

    @unpack ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec,
    ensemble_single_periodsum_vecs =
        data_arrs

    ensemble_nsims = size(ensemble_single_incarr, 3)

    logfile = open(logfilepath, "a")

    start_time = time()
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
                ews_threshold_burnin,
                ews_threshold_percentile,
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
                    styled"\t\tEWS hyperparameters\n\t\tEWS metric specification: {blue,inverse: $(ews_metric_specification.dirpath)}, End date type: {magenta: $(ews_enddate_type_str)}, EWS window: $(ews_threshold_window), EWS burn-in: {yellow: $(ews_threshold_burnin)}, EWS percentile: {magenta,inverse: $(ews_threshold_percentile)}, EWS consecutive thresholds: {yellow,inverse: $(ews_consecutive_thresholds)}"
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
                                percentiles = ews_threshold_percentile,
                                burn_in = ews_threshold_burnin,
                            )[2]

                            detection_index_arr[sim, j] = calculate_ews_trigger_index(
                                exceeds_threshold_arr[sim, j];
                                consecutive_thresholds = ews_consecutive_thresholds,
                            )

                            null_exceeds_threshold_arr[sim, j] = expanding_ews_thresholds(
                                null_ews_vals_vec[sim],
                                Symbol(ews_metric),
                                ews_threshold_window;
                                percentiles = ews_threshold_percentile,
                                burn_in = ews_threshold_burnin,
                            )[2]

                            null_detection_index_arr[sim, j] = calculate_ews_trigger_index(
                                null_exceeds_threshold_arr[sim, j];
                                consecutive_thresholds = ews_consecutive_thresholds,
                            )
                        end

                    else
                        failed_sims[sim] = sim
                        write(
                            logfile,
                            "$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))\n",
                        )
                        write(logfile, "Error:\t$(Try.unwrap_err(enddate))\n")
                        write(logfile, "Simulation:\t$(sim)\n\n")
                    end
                end

                for (j, ews_metric) in pairs(missing_ews_metric_vec)
                    true_positives = length(
                        filter(!isnothing, detection_index_arr[:, j])
                    )
                    true_negatives = length(
                        filter(isnothing, null_detection_index_arr[:, j])
                    )
                    sensitivity = true_positives / ensemble_nsims
                    specificity = true_negatives / ensemble_nsims
                    accuracy = (sensitivity + specificity) / 2

                    push!(
                        ews_df,
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
    end_time = time()
    return Try.Ok((; time = end_time - start_time, missing_runs = missing_runs))
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
    time_per_run_s = 0.11,
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

    missing_runs = nrow(missing_run_params_df)

    if missing_runs == 0
        return Try.Err("No missing simulations")
    end

    if missing_runs > 0
        nrun_time_s = missing_runs * time_per_run_s
        nrun_time_minutes = round(nrun_time_s / 60; digits = 2)
        nrun_time_message = if nrun_time_s < 10
            "less than 10 seconds"
        elseif nrun_time_s < 60
            "approximately $(round(nrun_time_s; digits = 0)) seconds"
        else
            "approximately $(nrun_time_minutes) minutes"
        end
        choice = request(
            "There are $(missing_runs) missing simulations. This is estimated to take $(nrun_time_message). Do you want to continue?",
            RadioMenu(["No", "Yes"]; ctrl_c_interrupt = false),
        )

        if choice != 2
            return Try.Err("User aborted")
        end

        println("Continuing ...")
    end

    missing_params_nt = create_missing_run_params_nt(missing_run_params_df)

    return Try.Ok(missing_params_nt)
end

function create_missing_run_params_nt(missing_run_params_df)
    pns = names(missing_run_params_df)
    output_pns = Symbol.("missing_" .* pns .* "_vec")
    vals = map(pns) do pn
        unique(missing_run_params_df[!, pn])
    end

    missing_runs = nrow(missing_run_params_df)

    return (; zip(output_pns, vals)..., missing_runs)
end
