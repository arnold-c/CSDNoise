using Dates: Dates
using DataFrames
using StyledStrings
using ProgressMeter
using UnPack: @unpack
using REPL: REPL
using REPL.TerminalMenus: RadioMenu, request
using Try: Try
using Match: Match

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
    return_df = true,
    specification_vec_tuples = (
        noise_specification = NoiseSpecification[],
        test_specification = IndividualTestSpecification[],
        percent_tested = Float64[],
        ews_metric_specification = EWSMetricSpecification[],
        ews_enddate_type = EWSEndDateType[],
        ews_threshold_window = Union{
            Type{ExpandingThresholdWindow},Type{RollingThresholdWindow}
        }[],
        ews_threshold_burnin = Union{Dates.Day,Dates.Year}[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_metric = String[],
    ),
    subset_optimal_parameters = [],
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
    ],
    disable_time_check = false,
    time_per_run_s = 0.08,
)
    if !isdir(filedir)
        mkpath(filedir)
    end

    ews_df = ews_hyperparam_gridsearch(
        specification_vecs,
        data_arrs;
        filedir = filedir,
        filename_base = gridsearch_filename_base,
        output_filepath = gridsearch_output_filepath,
        logfilepath = logfilepath,
        force = force,
        specification_vec_tuples = specification_vec_tuples,
        disable_time_check = disable_time_check,
        time_per_run_s = time_per_run_s,
        return_df = true,
    )

    load_filepath = get_most_recent_hyperparam_filepath(
        optimization_filename_base,
        filedir,
    )
    optimal_ews_df = filter_optimal_ews_hyperparam_gridsearch(
        ews_df;
        subset_optimal_parameters = subset_optimal_parameters,
        optimal_grouping_parameters = optimal_grouping_parameters,
    )

    if Try.isok(load_filepath) && !force
        previous_optimal_ews_df = load(Try.unwrap(load_filepath))["optimal_ews_df"]

        if optimal_ews_df == previous_optimal_ews_df
            @info "🟨 Previous optimal ews_df is the same as the current. No need to save a new version. 🟨"
            return optimal_ews_df
        end
    end

    @tagsave(
        optimization_output_filepath,
        Dict("optimal_ews_df" => optimal_ews_df)
    )
    @info "🟢 Saved optimal ews_df to $(optimization_output_filepath) 🟢"

    if return_df
        return optimal_ews_df
    end

    return nothing
end

function filter_optimal_ews_hyperparam_gridsearch(
    ews_df;
    subset_optimal_parameters = [],
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
    ],
)
    return map(
        collect(
            groupby(
                subset(ews_df, subset_optimal_parameters),
                optimal_grouping_parameters,
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
        ews_threshold_window = Union{
            Type{ExpandingThresholdWindow},Type{RollingThresholdWindow}
        }[],
        ews_threshold_burnin = Union{Dates.Day,Dates.Year}[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_metric = String[],
    ),
    disable_time_check = false,
    time_per_run_s = 0.08,
    return_df = true,
)
    load_filepath = get_most_recent_hyperparam_filepath(
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
        disable_time_check = disable_time_check,
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

function load_most_recent_hyperparam_file(
    filename_base,
    filedir,
)
    filepath = get_most_recent_hyperparam_filepath(
        filename_base,
        filedir,
    )

    if Try.iserr(filepath)
        return filepath
    end

    return load(Try.unwrap(filepath))
end

function get_most_recent_hyperparam_filepath(
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
        ews_threshold_window = Union{
            Type{ExpandingThresholdWindow},Type{RollingThresholdWindow}
        }[],
        ews_threshold_burnin = Union{Dates.Day,Dates.Year}[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_metric = String[],
    ),
    disable_time_check = false,
    time_per_run_s = 0.08,
)
    specification_vec_names = [
        Symbol(match(r"(.*)(_vec)$", string(pn)).captures[1]) for
        pn in propertynames(specification_vecs)
    ]

    @assert specification_vec_names ==
        [propertynames(specification_vec_tuples)...]

    missing_specification_vecs = check_missing_ews_hyperparameter_simulations(
        ews_df,
        specification_vecs;
        specification_vec_names = specification_vec_names,
        disable_time_check = disable_time_check,
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
    for noise_specification in missing_noise_specification_vec
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
    specification_vec_names = (
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_threshold_percentile,
        :ews_consecutive_thresholds,
        :ews_metric,
    ),
    disable_time_check = false,
    time_per_run_s = 0.08,
)
    # run_params_df = DataFrame(
    #     specification_vec_tuples
    # )
    #
    # for (
    #     noise_specification,
    #     test_specification,
    #     percent_tested,
    #     ews_metric_specification,
    #     ews_enddate_type,
    #     ews_threshold_window,
    #     ews_threshold_burnin,
    #     ews_threshold_percentile,
    #     ews_consecutive_thresholds,
    #     ews_metric,
    # ) in
    #     push!(
    #         run_params_df,
    #         (
    #             noise_specification,
    #             test_specification,
    #             percent_tested,
    #             ews_metric_specification,
    #             ews_enddate_type,
    #             ews_threshold_window,
    #             ews_threshold_burnin,
    #             ews_threshold_percentile,
    #             ews_consecutive_thresholds,
    #             ews_metric,
    #         ),
    #     )
    # end

    run_params_df = DataFrame(Iterators.product(specification_vecs...))
    rename!(run_params_df, specification_vec_names)

    missing_run_params_df = antijoin(
        run_params_df,
        ews_df;
        on = specification_vec_names,
    )

    missing_runs = nrow(missing_run_params_df)

    if missing_runs == 0
        return Try.Err("No missing simulations")
    end

    if missing_runs > 0 && !disable_time_check
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

function optimal_ews_heatmap_df(
    optimal_ews_df;
    tiebreaker_preference = "speed",
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
    ],
)
    tiebreaker_args = Match.@match tiebreaker_preference begin
        "speed" => (:ews_consecutive_thresholds, false)
        "specificity" => (:specificity, true)
        _ => error(
            "Invalid preference: $tiebreaker_preference. Please choose either \"speed\" or \"specificity\"."
        )
    end

    return map(
        collect(
            groupby(
                optimal_ews_df,
                optimal_grouping_parameters,
            ),
        ),
    ) do df
        sort(df, order(tiebreaker_args[1]; rev = tiebreaker_args[2]))[1, :]
    end |>
           x -> vcat(DataFrame.(x)...; cols = :union)
end

function optimal_ews_heatmap_plot(
    df;
    outcome = :accuracy,
    baseline_test = IndividualTestSpecification(1.0, 1.0, 0),
    colormap = :Blues,
    textcolorthreshold = 0.6,
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    outcome_str = titlecase(string(outcome))

    if !haskey(kwargs_dict, :plottitle)
        plottitle =
            "Heatmap: " * titlecase(string(outcome))
    end

    if !haskey(kwargs_dict, :subtitle)
        ews_enddate_type_str = split(string(df[1, :ews_enddate_type]), "::")[1]

        subtitle =
            "Noise: $(get_noise_magnitude_description(df[1, :noise_specification])), $(ews_enddate_type_str)" *
            "\nP = Percentile Threshold, C = Consecutive Thresholds, S = Specificity"
    end

    df[!, :test_sens] =
        getproperty.(df.test_specification, :sensitivity)
    df[!, :test_spec] = getproperty.(df.test_specification, :sensitivity)
    df[!, :test_result_lag] =
        getproperty.(df.test_specification, :test_result_lag)
    sort!(
        df,
        [:test_sens, :test_spec, :test_result_lag];
        rev = [true, true, false],
    )

    unique_tests = unique(df.test_specification)
    if mapreduce(x -> x == baseline_test, +, unique_tests) == 0
        error("default_test $baseline_test must be in unique_tests")
    end

    default_test_metric_order, mat = create_ews_heatmap_matrix(
        df,
        outcome,
    )

    threshold_percentile_matrix = create_ews_heatmap_matrix(
        df,
        :ews_threshold_percentile,
        default_test_metric_order,
    )

    consecutive_thresholds_df = create_ews_heatmap_matrix(
        df,
        :ews_consecutive_thresholds,
        default_test_metric_order,
    )

    specificity_df = create_ews_heatmap_matrix(
        df,
        :specificity,
        default_test_metric_order,
    )

    function test_axis_label(test)
        return "($(test.sensitivity), $(test.specificity), $(test.test_result_lag))"
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        subtitle = subtitle,
        xlabel = "Test Specification (Sensitivity, Specificity, Lag)",
        ylabel = "EWS Metric",
        xticks = (1:length(unique_tests), test_axis_label.(unique_tests)),
        yticks = (
            1:length(default_test_metric_order), default_test_metric_order
        ),
    )

    hmap = heatmap!(
        ax,
        mat;
        colormap = colormap,
        colorrange = (0, 1),
    )

    for j in axes(mat, 2), i in axes(mat, 1)
        val = mat[i, j]
        textcolor = abs(val) < textcolorthreshold ? :black : :white
        acc = round(mat[i, j]; digits = 2)
        perc = round(threshold_percentile_matrix[i, j]; digits = 2)
        consec = consecutive_thresholds_df[i, j]
        spec = round(specificity_df[i, j]; digits = 2)
        text!(
            ax,
            "Accuracy: $(acc)\nP: $(perc) C: $(consec) S: $(spec)";
            position = (i, j),
            color = textcolor,
            align = (:center, :center),
        )
    end

    limits!(
        ax,
        (0, length(unique_tests) + 1),
        (0, length(default_test_metric_order) + 1),
    )
    #
    Colorbar(fig[1, 2], hmap; label = outcome_str, width = 15, ticksize = 15)
    return fig
end

function create_ews_heatmap_matrix(
    df, outcome::Symbol, ews_metric_order
)
    ordered_df =
        unstack(
            select(df, [:ews_metric, :test_specification, outcome]),
            :test_specification,
            outcome,
        ) |>
        df -> df[indexin(ews_metric_order, df.ews_metric), :]

    return Matrix(ordered_df[:, 2:end])'
end

function create_ews_heatmap_matrix(
    df, outcome::Symbol
)
    ordered_df =
        unstack(
            select(df, [:ews_metric, :test_specification, outcome]),
            :test_specification,
            outcome,
        ) |>
        df -> sort(df, order(2; rev = false))

    default_test_metric_order = ordered_df.ews_metric

    return default_test_metric_order, Matrix(ordered_df[:, 2:end])'
end

function ews_survival_plot(
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    plottitle = "Survival",
    subtitle = "",
    ews_aggregation = Day(7),
    burnin = Year(5),
    endpoint_aggregation = Day(30),
    alpha = 1.0,
)
    @unpack detection_survival_vec, detection_indices_vec =
        detection_survival_vecs
    @unpack null_survival_vec, null_indices_vec = null_survival_vecs

    filtered_enddate_vec = filter(isinteger, enddate_vec)
    times =
        collect(1:Dates.days(ews_aggregation):maximum(filtered_enddate_vec)) ./
        365

    detection_survival_times = vcat(
        0,
        times[detection_indices_vec[1]],
        times[detection_indices_vec]...,
        times[end],
    )
    null_survival_times = vcat(
        0,
        times[null_indices_vec[1]],
        times[null_indices_vec]...,
        times[end],
    )

    detection_survival_vec = vcat(
        detection_survival_vec[1], detection_survival_vec[1],
        detection_survival_vec..., detection_survival_vec[end],
    )
    null_survival_vec = vcat(
        null_survival_vec[1],
        null_survival_vec[1],
        null_survival_vec...,
        null_survival_vec[end],
    )

    enddate_vec = div.(enddate_vec, Dates.days(endpoint_aggregation))
    unique_enddate_vec = sort(unique(enddate_vec))
    enddate_times =
        (unique_enddate_vec .* Dates.days(endpoint_aggregation)) ./ 365

    enddate_counts = map(unique_enddate_vec) do enddate
        sum(enddate_vec .== enddate)
    end

    fig = Figure()
    hist_ax = Axis(
        fig[1, 1];
        limits = (0, maximum(enddate_times), nothing, nothing),
    )
    surv_ax = Axis(
        fig[1, 1];
        title = plottitle,
        subtitle = subtitle,
        xlabel = "Time (Years)",
        ylabel = "Survival Numbers",
        xticks = 1:1:10,
        limits = (0, maximum(enddate_times), 0, nothing),
    )

    hidexdecorations!(hist_ax)
    hideydecorations!(hist_ax)

    barplot!(
        hist_ax,
        enddate_times,
        enddate_counts;
        color = (:darkgray, 1.0),
    )

    lines!(
        surv_ax,
        detection_survival_times,
        detection_survival_vec;
        color = (:blue, alpha),
        label = "Detection",
    )

    lines!(
        surv_ax,
        null_survival_times,
        null_survival_vec;
        color = (:red, alpha),
        label = "Null",
    )

    vlines!(
        surv_ax, Int64(Dates.value(burnin)); color = :black, linestyle = :dash
    )

    Legend(fig[1, 2], surv_ax; orientation = :vertical)

    return fig
end

function create_ews_survival_data(
    ews_optimal_simulation_df
)
    nsims = nrow(ews_optimal_simulation_df)
    unique_detection_indices = sort(
        filter(
            !isnothing, unique(ews_optimal_simulation_df.detection_index)
        ),
    )

    unique_null_detection_indices = sort(
        filter(
            !isnothing,
            unique(
                ews_optimal_simulation_df.null_detection_index
            ),
        ),
    )

    detection_indices_counts = map(
        ((i, v),) -> sum(ews_optimal_simulation_df.detection_index .== v),
        enumerate(unique_detection_indices),
    )

    null_detection_indices_counts = map(
        ((i, v),) -> sum(
            ews_optimal_simulation_df.null_detection_index .== v
        ),
        enumerate(unique_null_detection_indices),
    )

    @assert nsims ==
        sum(detection_indices_counts) +
            sum(ews_optimal_simulation_df.detection_index .== nothing)

    @assert nsims ==
        sum(null_detection_indices_counts) +
            sum(ews_optimal_simulation_df.null_detection_index .== nothing)

    detection_survival_vec = nsims .- cumsum(detection_indices_counts)
    null_survival_vec = nsims .- cumsum(null_detection_indices_counts)

    return (
        (;
            detection_survival_vec,
            detection_indices_vec = Int64.(unique_detection_indices),
        ),
        (;
            null_survival_vec,
            null_indices_vec = Int64.(unique_null_detection_indices),
        ),
    )
end

function simulate_ews_survival_data(
    optimal_ews_df,
    ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    thresholds;
    ews_metric = "mean",
    logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
)
    logfile = open(logfilepath, "a")

    survival_df = DataFrame(
    (
        test_specification = IndividualTestSpecification[],
        ews_metric = String[],
        enddate = Vector{Union{Int64,Try.Err}}(),
        detection_index = Vector{Union{Nothing,Int64}}(),
        null_detection_index = Vector{Union{Nothing,Int64}}(),
    )
)

    subset_optimal_ews_df = subset(
        optimal_ews_df,
        :ews_metric => ByRow(==(ews_metric)),
    )

    for n in [
        :noise_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
    ]
        @assert length(unique(subset_optimal_ews_df[:, n])) == 1
    end

    noise_specification = subset_optimal_ews_df[1, :noise_specification]
    percent_tested = subset_optimal_ews_df[1, :percent_tested]
    ews_metric_specification = subset_optimal_ews_df[
        1, :ews_metric_specification
    ]
    ews_enddate_type = subset_optimal_ews_df[1, :ews_enddate_type]
    ews_threshold_window = subset_optimal_ews_df[1, :ews_threshold_window]
    ews_threshold_burnin = subset_optimal_ews_df[1, :ews_threshold_burnin]

    ews_metric = subset_optimal_ews_df[1, :ews_metric]

    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )[1]

    enddate_vec = Vector{Union{Try.Ok,Try.Err}}(
        undef, size(ensemble_single_incarr, 3)
    )

    for sim in axes(ensemble_single_incarr, 3)
        enddate_vec[sim] = calculate_ews_enddate(
            thresholds[sim],
            ews_enddate_type,
        )
    end
    failed_sims = sum(Try.iserr.(enddate_vec))

    vec_of_testarr = Vector{Array{Float64,3}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_testarr = Vector{Array{Float64,3}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_ews_vals_vec = Vector{Vector{Union{Missing,EWSMetrics}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_ews_vals_vec = Vector{Vector{Union{Missing,EWSMetrics}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_exceed_thresholds = Vector{Vector{Matrix{Bool}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_exceed_thresholds = Vector{Vector{Matrix{Bool}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_threshold_percentiles = Vector{Vector{Matrix{Float64}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_threshold_percentiles = Vector{Vector{Matrix{Float64}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_detection_index_vec = Vector{Vector{Union{Nothing,Int64}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_detection_index_vec = Vector{Vector{Union{Nothing,Int64}}}(
        undef, nrow(subset_optimal_ews_df)
    )

    for (i, df_row) in pairs(eachrow(subset_optimal_ews_df))
        test_specification = df_row[:test_specification]
        ews_threshold_percentile = df_row[:ews_threshold_percentile]
        ews_consecutive_thresholds = df_row[:ews_consecutive_thresholds]

        testarr = create_testing_arrs(
            ensemble_single_incarr,
            noisearr,
            percent_tested,
            test_specification,
        )

        vec_of_testarr[i] = testarr

        null_testarr = create_testing_arrs(
            null_single_incarr,
            noisearr,
            percent_tested,
            test_specification,
        )
        vec_of_null_testarr[i] = null_testarr

        ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
            undef, size(testarr, 3)
        )
        null_ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
            undef, size(testarr, 3)
        )
        fill!(ews_vals_vec, missing)
        fill!(null_ews_vals_vec, missing)

        exceeds_threshold_vec = Vector{Matrix{Bool}}(
            undef, size(testarr, 3)
        )
        null_exceeds_threshold_vec = Vector{Matrix{Bool}}(
            undef, size(testarr, 3)
        )
        threshold_percentiles_vec = Vector{Matrix{Float64}}(
            undef, size(testarr, 3)
        )
        null_threshold_percentiles_vec = Vector{Matrix{Float64}}(
            undef, size(testarr, 3)
        )

        detection_index_vec = Vector{Union{Nothing,Int64}}(
            undef, size(testarr, 3)
        )
        null_detection_index_vec = Vector{Union{Nothing,Int64}}(
            undef, size(testarr, 3)
        )
        fill!(detection_index_vec, nothing)
        fill!(null_detection_index_vec, nothing)

        for sim in axes(testarr, 3)
            if Try.isok(enddate_vec[sim])
                enddate = Try.unwrap(enddate_vec[sim])
                ews_vals_vec[sim] = EWSMetrics(
                    ews_metric_specification,
                    @view(testarr[1:enddate, 5, sim])
                )

                null_ews_vals_vec[sim] = EWSMetrics(
                    ews_metric_specification,
                    @view(null_testarr[1:enddate, 5, sim])
                )

                exceeds_threshold_vec[sim] = expanding_ews_thresholds(
                    ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[2]

                threshold_percentiles_vec[sim] = expanding_ews_thresholds(
                    ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[1]

                detection_index_vec[sim] = calculate_ews_trigger_index(
                    exceeds_threshold_vec[sim];
                    consecutive_thresholds = ews_consecutive_thresholds,
                )

                null_exceeds_threshold_vec[sim] = expanding_ews_thresholds(
                    null_ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[2]

                null_threshold_percentiles_vec[sim] = expanding_ews_thresholds(
                    null_ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[1]

                null_detection_index_vec[sim] = calculate_ews_trigger_index(
                    null_exceeds_threshold_vec[sim];
                    consecutive_thresholds = ews_consecutive_thresholds,
                )

            else
                write(
                    logfile,
                    "$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))\n",
                )
                write(logfile, "Error:\t$(Try.unwrap_err(enddate_vec[sim]))\n")
                write(logfile, "Simulation:\t$(sim)\n\n")
            end
        end

        vec_of_ews_vals_vec[i] = ews_vals_vec
        vec_of_null_ews_vals_vec[i] = null_ews_vals_vec
        vec_of_exceed_thresholds[i] = exceeds_threshold_vec
        vec_of_null_exceed_thresholds[i] = null_exceeds_threshold_vec
        vec_of_threshold_percentiles[i] = threshold_percentiles_vec
        vec_of_null_threshold_percentiles[i] = null_threshold_percentiles_vec
        vec_of_detection_index_vec[i] = detection_index_vec
        vec_of_null_detection_index_vec[i] = null_detection_index_vec

        append!(
            survival_df,
            DataFrame(;
                test_specification = repeat(
                    [test_specification], length(detection_index_vec)
                ),
                enddate = Try.unwrap_or_else.(Try.Err, enddate_vec),
                ews_metric = repeat([ews_metric], length(detection_index_vec)),
                detection_index = detection_index_vec,
                null_detection_index = null_detection_index_vec,
            );
            promote = true,
        )
    end

    return survival_df,
    (
        vec_of_testarr,
        vec_of_null_testarr,
        vec_of_ews_vals_vec,
        vec_of_null_ews_vals_vec,
        vec_of_exceed_thresholds,
        vec_of_null_exceed_thresholds,
        vec_of_threshold_percentiles,
        vec_of_null_threshold_percentiles,
        vec_of_detection_index_vec,
        vec_of_null_detection_index_vec,
    )
end
