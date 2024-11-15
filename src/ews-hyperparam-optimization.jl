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
        :ews_threshold_window,
        :ews_threshold_burnin,
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
        :ews_threshold_window,
        :ews_threshold_burnin,
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
        :ews_threshold_window,
        :ews_threshold_burnin,
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
    colormap = :RdBu,
    colorrange = [0.2, 0.8],
    textcolorthreshold = (0.4, 0.68),
    kwargs...,
)
    kwargs_dict = Dict{Symbol,Any}(kwargs)

    outcome_str = titlecase(string(outcome))

    if !haskey(kwargs_dict, :plottitle)
        plottitle =
            "Heatmap: " * titlecase(string(outcome))
    else
        plottitle = kwargs_dict[:plottitle]
    end

    if !haskey(kwargs_dict, :subtitle)
        ews_enddate_type_str = split(string(df[1, :ews_enddate_type]), "::")[1]

        subtitle =
            "Noise: $(get_noise_magnitude_description(df[1, :noise_specification])), $(ews_enddate_type_str)" *
            "\nP = Percentile Threshold, C = Consecutive Thresholds, S = Specificity"
    else
        subtitle = kwargs_dict[:subtitle]
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

    if minimum(mat) < colorrange[1]
        colorrange[1] = floor(minimum(mat); digits = 1)
    end
    if maximum(mat) > colorrange[2]
        colorrange[2] = ceil(maximum(mat); digits = 1)
    end

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
        colorrange = colorrange,
    )

    for j in axes(mat, 2), i in axes(mat, 1)
        val = mat[i, j]
        if length(textcolorthreshold) == 1
            textcolor = val <= textcolorthreshold ? :black : :white
        elseif length(textcolorthreshold) == 2
            textcolor =
                if val >= textcolorthreshold[1] && val <= textcolorthreshold[2]
                    :black
                else
                    :white
                end
        else
            error(
                "variable `textcolorthreshold` should be length 1 or 2. Instead received length $(length(textcolorthreshold))"
            )
        end
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

function simulate_and_plot_ews_survival(
    optimal_heatmap_df,
    ews_metric_specification,
    ews_enddate_type,
    ews_threshold_burnin,
    ews_threshold_window,
    noise_specification,
    ensemble_specification,
    individual_test_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec;
    ews_metric = "mean",
    plottitle = "Survival",
    subtitle = "",
)
    subset_df = subset(
        optimal_heatmap_df,
        :ews_metric_specification => ByRow(==(ews_metric_specification)),
        :ews_enddate_type => ByRow(==(ews_enddate_type)),
        :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
        :ews_threshold_window => ByRow(==(ews_threshold_window)),
        :noise_specification => ByRow(==(noise_specification)),
    )

    return simulate_and_plot_ews_survival(
        subset_df,
        ews_metric_specification,
        ews_threshold_burnin,
        ensemble_specification,
        individual_test_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = ews_metric,
        plottitle = plottitle,
        subtitle = subtitle,
    )
end

function simulate_and_plot_ews_survival(
    subset_df,
    ews_metric_specification,
    ews_threshold_burnin,
    ensemble_specification,
    individual_test_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec;
    ews_metric = "mean",
    plottitle = "Survival",
    subtitle = "",
)
    test_subset_df = subset(
        subset_df,
        :test_specification => ByRow(==(individual_test_specification)),
    )

    survival_df = simulate_ews_survival_data(
        test_subset_df,
        ensemble_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = ews_metric,
    )[1]

    detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
        survival_df
    )

    @unpack aggregation = ews_metric_specification

    return ews_survival_plot(
        detection_survival_vecs,
        null_survival_vecs,
        survival_df.enddate;
        ews_aggregation = aggregation,
        burnin = ews_threshold_burnin,
        plottitle = plottitle,
        subtitle = subtitle,
    )
end

function ews_survival_plot(
    survival_df::T1;
    test_specification_vec = [
        IndividualTestSpecification(1.0, 1.0, 0),
        IndividualTestSpecification(0.9, 0.9, 0),
    ],
    noise_specification_vec = [
        PoissonNoiseSpecification(1.0),
        PoissonNoiseSpecification(7.0),
        DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.8734),
        DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.102),
    ],
    endpoint_aggregation = Dates.Day(30),
    linestyle_vec = [:solid, :dot],
    alpha = 1.0,
    trim_burnin = true,
    plottitle = "Survival",
    subtitle = "",
    nbanks = length(test_specification_vec) * length(linestyle_vec),
    legend_rowsize = Makie.Relative(0.05),
) where {T1<:DataFrames.DataFrame}
    if sum(
        filter(
            !isnothing,
            indexin(
                unique(survival_df.test_specification),
                test_specification_vec,
            ),
        ),
    ) != sum(eachindex(test_specification_vec))
        error(
            "Not all test specified are present in the dataframe.\nTrying to plot survival curves for:\n$(test_specification_vec).\nFound these in the dataframe:\n$(unique(survival_df.test_specification))"
        )
    end

    if length(test_specification_vec) > length(linestyle_vec)
        error(
            "There are more tests specified ($(length(test_specification_vec))) than line styles ($(length(linestyle_vec)))."
        )
    end

    if sum(
        filter(
            !isnothing,
            indexin(
                unique(survival_df.noise_specification),
                noise_specification_vec,
            ),
        ),
    ) != sum(eachindex(noise_specification_vec))
        error(
            "Not all noise structures specified are present in the dataframe.\nTrying to plot survival curves for:\n$(noise_specification_vec).\nFound these in the dataframe:\n$(unique(survival_df.noise_specification))"
        )
    end

    @assert length(unique(survival_df.ews_metric_specification)) == 1
    @assert length(unique(survival_df.ews_threshold_burnin)) == 1
    ews_aggregation = survival_df.ews_metric_specification[1].aggregation
    burnin = survival_df.ews_threshold_burnin[1]

    subset_survival_df = subset(
        survival_df,
        :test_specification => ByRow(in(test_specification_vec)),
        :noise_specification => ByRow(in(noise_specification_vec)),
    )

    noise_grouped_dfs = groupby(subset_survival_df, :noise_specification)

    num_noise = length(noise_specification_vec)

    if num_noise > 4
        error(
            "Trying to plot too many metric facets $(length(ewsmetric_grouped_dfs)). Max allowed is 4"
        )
    end

    fig = Figure()

    for (noise_num, noise_gdf) in enumerate(noise_grouped_dfs)
        noise_description = noise_table_description(
            noise_gdf.noise_specification[1]
        )
        ax_position = Match.@match noise_num begin
            1 => (1, 1)
            2 => (1, 2)
            3 => (2, 1)
            4 => (2, 2)
        end
        gl = fig[ax_position...] = GridLayout()
        surv_ax = nothing
        for (i, test_specification) in pairs(test_specification_vec)
            detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
                subset(
                    noise_gdf,
                    :test_specification => ByRow(==(test_specification)),
                ),
            )

            times,
            enddate_times,
            enddate_counts,
            detection_survival_times,
            detection_survival_vec,
            null_survival_times,
            null_survival_vec,
            nsims = prepare_survival_facet_params(
                detection_survival_vecs,
                null_survival_vecs,
                noise_gdf.enddate;
                ews_aggregation = ews_aggregation,
                endpoint_aggregation = endpoint_aggregation,
            )

            if i == 1
                survival_plot_histogram!(
                    gl,
                    enddate_times,
                    enddate_counts;
                    trim_burnin = trim_burnin,
                    burnin = burnin,
                )
                surv_ax = Axis(
                    gl[1, 1];
                    title = noise_description,
                    xlabel = "Time (Years)",
                    ylabel = "Survival Numbers",
                    xticks = 1:1:10,
                    limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
                )
            end

            survival_plot_lines!(
                surv_ax,
                times,
                detection_survival_times,
                detection_survival_vec,
                null_survival_times,
                null_survival_vec;
                alpha = alpha,
                nsims = nsims,
                trim_burnin = trim_burnin,
                burnin = burnin,
                linestyle = linestyle_vec[i],
            )
        end
        if ax_position[2] != 1
            hideydecorations!(surv_ax)
        end

        nrows = ceil(num_noise / 2)
        if ax_position[1] != nrows && ax_position[2] + 2 <= num_noise
            hidexdecorations!(surv_ax)
        end
    end

    Legend(
        fig[0, :],
        vcat(
            [
                LineElement(; linestyle = style) for
                style in linestyle_vec[1:length(test_specification_vec)]
            ],
            [
                PolyElement(; color = col) for
                col in [:red, :blue]
            ],
        ),
        vcat(
            get_test_description.(test_specification_vec),
            ["Null", "Emergent"],
        ),
        "";
        nbanks = nbanks,
    )

    rowsize!(fig.layout, 0, legend_rowsize)
    if num_noise == 1
        colsize!(fig.layout, 1, Makie.Relative(1))
    end

    return fig
end

# function ews_survival_facet!(
#     gl,
#     gdf::T1;
#     facet_title = "Survival",
#     subtitle = "",
#     ews_aggregation = Day(7),
#     burnin = Year(5),
#     endpoint_aggregation = Day(30),
#     alpha = 1.0,
#     trim_burnin = true,
# ) where {T1<:DataFrames.DataFrame}
#     return nothing
# end

function ews_survival_plot(
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    facet_title = "Survival",
    ews_aggregation = Day(7),
    burnin = Year(5),
    endpoint_aggregation = Day(30),
    alpha = 1.0,
    trim_burnin = true,
)
    fig = Figure()

    gl = fig[1, 1] = GridLayout()

    ews_survival_facet!(
        gl,
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        facet_title = facet_title,
        ews_aggregation = ews_aggregation,
        burnin = burnin,
        endpoint_aggregation = endpoint_aggregation,
        alpha = alpha,
        trim_burnin = trim_burnin,
    )

    Legend(
        fig[1, 2],
        contents(gl)[2];
        orientation = :vertical,
    )

    return fig
end

function ews_survival_facet!(
    gl,
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    facet_title = "Survival",
    ews_aggregation = Day(7),
    burnin = Year(5),
    endpoint_aggregation = Day(30),
    alpha = 1.0,
    trim_burnin = true,
)
    times,
    enddate_times,
    enddate_counts,
    detection_survival_times,
    detection_survival_vec,
    null_survival_times,
    null_survival_vec,
    nsims = prepare_survival_facet_params(
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        ews_aggregation = ews_aggregation,
        endpoint_aggregation = endpoint_aggregation,
    )

    survival_plot_histogram!(
        gl,
        enddate_times,
        enddate_counts;
        trim_burnin = trim_burnin,
        burnin = burnin,
    )

    survival_plot_lines!(
        gl,
        enddate_times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        alpha = alpha,
        nsims = nsims,
        facet_title = facet_title,
        trim_burnin = trim_burnin,
        burnin = burnin,
    )

    return nothing
end

function prepare_survival_facet_params(
    detection_survival_vecs,
    null_survival_vecs,
    enddate_vec;
    ews_aggregation = Day(7),
    endpoint_aggregation = Day(30),
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
        times[detection_indices_vec]...,
        times[end],
    )
    null_survival_times = vcat(
        0,
        times[null_indices_vec]...,
        times[end],
    )

    nsims = detection_survival_vec[1]

    enddate_vec = div.(enddate_vec, Dates.days(endpoint_aggregation))
    unique_enddate_vec = sort(unique(enddate_vec))
    enddate_times =
        (unique_enddate_vec .* Dates.days(endpoint_aggregation)) ./ 365

    enddate_counts = map(unique_enddate_vec) do enddate
        sum(enddate_vec .== enddate)
    end

    return times, enddate_times,
    enddate_counts,
    detection_survival_times,
    detection_survival_vec,
    null_survival_times,
    null_survival_vec,
    nsims
end

function survival_plot_histogram!(
    gl,
    enddate_times,
    enddate_counts;
    trim_burnin = true,
    burnin = Dates.Year(5),
)
    hist_ax = Axis(
        gl[1, 1];
        limits = (0, maximum(enddate_times), nothing, nothing),
    )

    hidexdecorations!(hist_ax)
    hideydecorations!(hist_ax)

    barplot!(
        hist_ax,
        enddate_times,
        enddate_counts;
        color = (:darkgray, 1.0),
    )

    if trim_burnin
        xlims!(hist_ax, Dates.value(burnin), enddate_times[end])
    end

    return nothing
end

function survival_plot_lines!(
    gl::T1,
    times,
    detection_survival_times,
    detection_survival_vec,
    null_survival_times,
    null_survival_vec;
    linestyle = :solid,
    alpha = 1.0,
    nsims = 100,
    facet_title = "Survival",
    trim_burnin = true,
    burnin = Dates.Year(5),
) where {T1<:Makie.GridLayout}
    surv_ax = Axis(
        gl[1, 1];
        title = facet_title,
        xlabel = "Time (Years)",
        ylabel = "Survival Numbers",
        xticks = 1:1:10,
        limits = (0, maximum(times), 0, ceil(1.1 * nsims)),
    )
    survival_plot_lines!(
        surv_ax,
        times,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec;
        linestyle = linestyle,
        alpha = alpha,
        nsims = nsims,
        trim_burnin = trim_burnin,
        burnin = burnin,
    )
    return nothing
end

function survival_plot_lines!(
    surv_ax::T1,
    times,
    detection_survival_times,
    detection_survival_vec,
    null_survival_times,
    null_survival_vec;
    linestyle = :solid,
    alpha = 1.0,
    nsims = 100,
    facet_title = "Survival",
    trim_burnin = true,
    burnin = Dates.Year(5),
) where {T1<:Makie.Axis}
    lines!(
        surv_ax,
        detection_survival_times,
        detection_survival_vec;
        color = (:blue, alpha),
        linestyle = linestyle,
        label = "Detection",
    )

    lines!(
        surv_ax,
        null_survival_times,
        null_survival_vec;
        color = (:red, alpha),
        linestyle = linestyle,
        label = "Null",
    )

    if trim_burnin
        xlims!(surv_ax, Dates.value(burnin), times[end])
    else
        vlines!(
            surv_ax, Int64(Dates.value(burnin)); color = :black,
            linestyle = :dash,
        )
    end
    return nothing
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

    detection_survival_vec, detection_indices_vec = calculate_detection_indices_and_survival(
        ews_optimal_simulation_df.detection_index,
        unique_detection_indices,
        nsims,
    )

    null_survival_vec, null_indices_vec = calculate_detection_indices_and_survival(
        ews_optimal_simulation_df.null_detection_index,
        unique_null_detection_indices,
        nsims,
    )

    return (
        (; detection_survival_vec, detection_indices_vec),
        (; null_survival_vec, null_indices_vec),
    )
end

function calculate_detection_indices_and_survival(
    detection_vec,
    detection_indices,
    nsims,
)
    if isempty(detection_indices)
        @assert nsims ==
            sum(detection_vec .== nothing)

        survival_vec = [nsims, nsims]
        indices_vec = Int64[]
    else
        detection_indices_counts = map(
            ((i, v),) -> sum(
                detection_vec .== v
            ),
            enumerate(detection_indices),
        )

        @assert nsims ==
            sum(detection_indices_counts) +
                sum(detection_vec .== nothing)

        survival_vec = nsims .- cumsum(detection_indices_counts)

        survival_vec = vcat(
            nsims, nsims, survival_vec..., survival_vec[end]
        )
        detection_indices = Int64.(detection_indices)

        indices_vec = vcat(
            detection_indices[1],
            detection_indices,
        )
    end

    return survival_vec, indices_vec
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
        ews_metric_specification = EWSMetricSpecification[],
        ews_threshold_burnin = Union{<:Dates.Day,<:Dates.Year}[],
        noise_specification = Union{
            <:PoissonNoiseSpecification,<:DynamicalNoiseSpecification
        }[],
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
                ews_metric_specification = repeat(
                    [ews_metric_specification], length(detection_index_vec)
                ),
                ews_threshold_burnin = repeat(
                    [ews_threshold_burnin], length(detection_index_vec)
                ),
                noise_specification = repeat(
                    [noise_specification], length(detection_index_vec)
                ),
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
    ),
    noisearr
end
