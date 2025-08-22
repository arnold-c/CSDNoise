using Dates: Dates, Day
using DataFrames
using ProgressMeter
using UnPack: @unpack
using MultistartOptimization
using NLopt
using FLoops: FLoops
using StyledStrings
using Try: Try
using DrWatson: @tagsave
using JLD2

"""
    OptimizationTracker

Mutable struct to track the best solution and its metrics during optimization.
"""
mutable struct OptimizationTracker
    best_loss::Float64
    best_accuracy::Float64
    best_sensitivity::Float64
    best_specificity::Float64
    best_params::Vector{Float64}

    function OptimizationTracker()
        return new(Inf, 0.0, 0.0, 0.0, Float64[])
    end
end

"""
    ews_multistart_optimization(specification_vecs, data_arrs; kwargs...)

Optimize EWS hyperparameters using multistart optimization with Sobol sequences.
Alternative to grid search that scales better with parameter dimensionality.

# Arguments
- `specification_vecs`: Named tuple containing specification vectors
- `data_arrs`: Named tuple containing ensemble data arrays

# Keyword Arguments
- `filedir`: Directory for saving results
- `optimization_filename_base`: Base filename for results
- `n_sobol_points`: Number of Sobol sequence points per scenario (default: 100)
- `local_algorithm`: NLopt local optimization algorithm (default: NLopt.LN_BOBYQA)
- `maxeval`: Maximum evaluations per local optimization (default: 1000)
- `xtol_rel`: Relative tolerance for parameter convergence (default: 1e-3)
- `xtol_abs`: Absolute tolerance for parameter convergence (default: 1e-3)
- `ftol_rel`: Relative tolerance for function convergence (default: 1e-4)
- `executor`: Executor for loops (default: `FLoops.SequentialEx()`)
- `percentile_bounds`: Bounds for threshold percentile (default: (0.5, 0.99))
- `consecutive_bounds`: Bounds for consecutive thresholds (default: (1.0, 10.0))
- `force`: Force recomputation even if results exist (default: false)
- `return_df`: Return DataFrame of results (default: true)
- `save_results`: Save results to file (default: true)
- `verbose`: Print progress information (default: true)

# Returns
- DataFrame with optimal parameters for each scenario
"""
function ews_multistart_optimization(
        specification_vecs,
        data_arrs;
        # File management
        filedir = outdir("ensemble", "ews-multistart-optimization"),
        optimization_filename_base = "ews-multistart-optimization.jld2",
        optimization_output_filepath = joinpath(
            filedir,
            string(Dates.now()) * "_" * optimization_filename_base,
        ),
        # Optimization configuration
        n_sobol_points = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
        executor = FLoops.SequentialEx(),
        # Parameter bounds
        percentile_bounds = (0.5, 0.99),
        consecutive_bounds = (1.0, 10.0),
        # Control options
        force = false,
        return_df = true,
        save_results = true,
        verbose = true,
    )
    if !isdir(filedir)
        mkpath(filedir)
    end

    # Check for previous results
    if !force
        previous_results = load_previous_multistart_results(
            filedir, optimization_filename_base
        )
        if !isnothing(previous_results)
            @info "Loading previous optimization results"
            return previous_results
        end
    end

    # Create scenario specifications
    scenarios = create_optimization_scenarios(specification_vecs)
    n_scenarios = length(scenarios)

    if verbose
        println(styled"{green:Starting Multistart Optimization}")
        println(styled"Scenarios to optimize: {yellow:$n_scenarios}")
        println(styled"Sobol points per scenario: {blue:$n_sobol_points}")
    end

    # Initialize results storage
    results = Vector{NamedTuple}(undef, n_scenarios)

    # Setup progress tracking
    prog = Progress(n_scenarios; desc = "Optimizing scenarios: ")

    # Define parameter bounds (only percentile and consecutive thresholds)
    bounds = (
        lowers = [percentile_bounds[1], consecutive_bounds[1]],
        uppers = [percentile_bounds[2], consecutive_bounds[2]],
    )

    # Optimization configuration
    optim_config = (
        n_sobol_points = n_sobol_points,
        local_algorithm = local_algorithm,
        maxeval = maxeval,
        xtol_rel = xtol_rel,
        xtol_abs = xtol_abs,
        ftol_rel = ftol_rel,
    )

    # Run optimization for each scenario
    FLoops.@floop executor for (idx, scenario) in enumerate(scenarios)
        result = optimize_single_scenario(
            scenario,
            data_arrs,
            bounds,
            optim_config
        )
        results[idx] = result
        next!(prog)
    end

    # Convert results to DataFrame
    results_df = create_results_dataframe(results, scenarios)

    # Save results
    if save_results
        @tagsave(
            optimization_output_filepath,
            Dict("multistart_optimization_df" => results_df)
        )
        @info "Saved optimization results to $optimization_output_filepath"
    end

    return return_df ? results_df : nothing
end

"""
    optimize_single_scenario(scenario, data_arrs, bounds, config)

Optimize EWS parameters for a single scenario using multistart optimization.
"""
function optimize_single_scenario(
        scenario::NamedTuple,
        data_arrs::NamedTuple,
        bounds::NamedTuple,
        config::NamedTuple
    )
    @unpack ews_metric = scenario

    simulated_data = simulate_timeseries_arrs(scenario, data_arrs)

    # Create tracker instance for this scenario
    tracker = OptimizationTracker()

    # Create objective function closure that updates tracker
    objective = params -> ews_objective_function_with_tracking(
        params,
        scenario,
        simulated_data,
        tracker
    )

    # Setup multistart problem
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        bounds.lowers,
        bounds.uppers
    )

    # Configure local optimization method
    local_method = MultistartOptimization.NLopt_local_method(
        config.local_algorithm;
        xtol_rel = config.xtol_rel,
        xtol_abs = config.xtol_abs,
        # ftol_rel = config.ftol_rel,
        maxeval = config.maxeval,
    )

    # Configure multistart method (TikTak uses Sobol sequences)
    multistart_method = MultistartOptimization.TikTak(config.n_sobol_points)

    # Run optimization
    result = MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    # Extract optimal parameters from tracker (which has the best metrics)
    optimal_params = map_continuous_to_ews_parameters(tracker.best_params)

    return (
        optimal_threshold_percentile = optimal_params.threshold_percentile,
        optimal_consecutive_thresholds = optimal_params.consecutive_thresholds,
        accuracy = tracker.best_accuracy,
        sensitivity = tracker.best_sensitivity,
        specificity = tracker.best_specificity,
        # n_evaluations = result.n_calls,
        # convergence_status = string(result.ret_code),
    )
end

"""
    ews_objective_function_with_tracking(params_vec, scenario, data_arrs, tracker)

Objective function for EWS parameter optimization that tracks the best solution.
Returns 1 - accuracy for minimization and updates tracker with best metrics.
"""
function ews_objective_function_with_tracking(
        params_vec::Vector{Float64},
        scenario::NamedTuple,
        data_arrs::NamedTuple,
        tracker::OptimizationTracker
    )
    @unpack ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window, ews_threshold_burnin, ews_metric = scenario

    @unpack ensemble_single_incarr,
        testarr,
        null_testarr,
        thresholds = data_arrs

    # Map continuous parameters to EWS parameters
    ews_params = map_continuous_to_ews_parameters(params_vec)

    # Calculate accuracy
    ensemble_nsims = size(ensemble_single_incarr, 3)
    true_positives = 0
    true_negatives = 0

    for sim in 1:ensemble_nsims
        enddate = calculate_ews_enddate(thresholds[sim], ews_enddate_type)

        if Try.isok(enddate)
            enddate_val = Try.unwrap(enddate)

            # Calculate EWS metrics
            ews_vals = EWSMetrics(
                ews_metric_specification,
                @view(testarr[1:enddate_val, 5, sim])
            )

            null_ews_vals = EWSMetrics(
                ews_metric_specification,
                @view(null_testarr[1:enddate_val, 5, sim])
            )

            # Check threshold exceedances
            exceeds_threshold = expanding_ews_thresholds(
                ews_vals,
                Symbol(ews_metric),
                ews_threshold_window;
                percentiles = ews_params.threshold_percentile,
                burn_in = ews_threshold_burnin,
            )[2]

            detection_index = calculate_ews_trigger_index(
                exceeds_threshold;
                consecutive_thresholds = ews_params.consecutive_thresholds,
            )

            null_exceeds_threshold = expanding_ews_thresholds(
                null_ews_vals,
                Symbol(ews_metric),
                ews_threshold_window;
                percentiles = ews_params.threshold_percentile,
                burn_in = ews_threshold_burnin,
            )[2]

            null_detection_index = calculate_ews_trigger_index(
                null_exceeds_threshold;
                consecutive_thresholds = ews_params.consecutive_thresholds,
            )

            # Update counts
            if !isnothing(detection_index)
                true_positives += 1
            end
            if isnothing(null_detection_index)
                true_negatives += 1
            end
        end
    end

    # Calculate metrics
    sensitivity = calculate_sensitivity(true_positives, ensemble_nsims)
    specificity = calculate_specificity(true_negatives, ensemble_nsims)
    accuracy = calculate_accuracy(sensitivity, specificity)
    loss = 1.0 - accuracy

    # Update tracker if this is the best solution so far
    if loss < tracker.best_loss
        tracker.best_loss = loss
        tracker.best_accuracy = accuracy
        tracker.best_sensitivity = sensitivity
        tracker.best_specificity = specificity
        tracker.best_params = copy(params_vec)
    end

    return loss  # Return scalar loss for optimizer
end

function simulate_timeseries_arrs(
        scenario::NamedTuple,
        data_arrs::NamedTuple
    )
    @unpack noise_specification, test_specification, percent_tested,
        ews_enddate_type = scenario

    @unpack ensemble_single_incarr, null_single_incarr,
        ensemble_specification, ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs = data_arrs

    # Create noise array
    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )[1]

    # Create test arrays
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

    # Select appropriate thresholds
    thresholds = SumTypes.@cases ews_enddate_type begin
        [Reff_start, Reff_end] => ensemble_single_Reff_thresholds_vec
        [Outbreak_start, Outbreak_end, Outbreak_middle] => ensemble_single_periodsum_vecs
    end

    return (; ensemble_single_incarr, testarr, null_testarr, thresholds)
end

"""
    map_continuous_to_ews_parameters(params_vec)

Map continuous optimization parameters to discrete EWS parameters.
Note: burnin is now passed as part of the scenario, not optimized.
"""
function map_continuous_to_ews_parameters(params_vec::Vector{Float64})
    return (
        threshold_percentile = params_vec[1],  # Already in correct range
        consecutive_thresholds = round(Int, params_vec[2]),
    )
end

"""
    create_optimization_scenarios(specification_vecs)

Create all scenario combinations for optimization.
"""
function create_optimization_scenarios(specification_vecs)
    @unpack noise_specification_vec, test_specification_vec,
        percent_tested_vec, ews_metric_specification_vec,
        ews_enddate_type_vec, ews_threshold_window_vec,
        ews_threshold_burnin_vec, ews_metric_vec = specification_vecs

    scenarios = []

    for noise_spec in noise_specification_vec
        for test_spec in test_specification_vec
            for percent_tested in percent_tested_vec
                for ews_metric_spec in ews_metric_specification_vec
                    for ews_enddate_type in ews_enddate_type_vec
                        for ews_window in ews_threshold_window_vec
                            for ews_burnin in ews_threshold_burnin_vec
                                for ews_metric in ews_metric_vec
                                    push!(
                                        scenarios, (
                                            noise_specification = noise_spec,
                                            test_specification = test_spec,
                                            percent_tested = percent_tested,
                                            ews_metric_specification = ews_metric_spec,
                                            ews_enddate_type = ews_enddate_type,
                                            ews_threshold_window = ews_window,
                                            ews_threshold_burnin = ews_burnin,
                                            ews_metric = ews_metric,
                                        )
                                    )
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return scenarios
end

"""
    create_results_dataframe(results, scenarios)

Convert optimization results to DataFrame format.
"""
function create_results_dataframe(results::Vector, scenarios::Vector)
    @assert length(scenarios) == length(results)
    df = DataFrame(
        noise_specification = [s.noise_specification for s in scenarios],
        test_specification = [s.test_specification for s in scenarios],
        percent_tested = [s.percent_tested for s in scenarios],
        ews_metric_specification = [s.ews_metric_specification for s in scenarios],
        ews_enddate_type = [s.ews_enddate_type for s in scenarios],
        ews_threshold_window = [s.ews_threshold_window for s in scenarios],
        ews_metric = [s.ews_metric for s in scenarios],
        ews_threshold_percentile = [r.optimal_threshold_percentile for r in results],
        ews_consecutive_thresholds = [r.optimal_consecutive_thresholds for r in results],
        ews_threshold_burnin = [s.ews_threshold_burnin for s in scenarios],
        accuracy = [r.accuracy for r in results],
        sensitivity = [r.sensitivity for r in results],
        specificity = [r.specificity for r in results],
        # n_evaluations = [r.n_evaluations for r in results],
        # convergence_status = [r.convergence_status for r in results],
    )


    return df
end

"""
    load_previous_multistart_results(filedir, filename_base)

Load previous multistart optimization results if they exist.
"""
function load_previous_multistart_results(filedir::String, filename_base::String)
    if !isdir(filedir)
        return nothing
    end

    # Get most recent file
    load_filepath = get_most_recent_hyperparam_filepath(filename_base, filedir)

    if Try.iserr(load_filepath)
        return nothing
    end

    try
        data = JLD2.load(Try.unwrap(load_filepath))
        return get(data, "multistart_optimization_df", nothing)
    catch
        return nothing
    end
end

"""
    filter_optimal_multistart_results(multistart_df; kwargs...)

Filter multistart optimization results to find optimal parameters.
Similar to filter_optimal_ews_hyperparam_gridsearch but for multistart results.
"""
function filter_optimal_multistart_results(
        multistart_df;
        subset_optimal_parameters = [],
        optimal_grouping_parameters = [
            :noise_specification,
            :test_specification,
            :percent_tested,
            :ews_metric_specification,
            :ews_enddate_type,
            :ews_metric,
            :ews_threshold_window,
        ],
    )
    return map(
        collect(
            groupby(
                subset(multistart_df, subset_optimal_parameters),
                optimal_grouping_parameters,
            ),
        ),
    ) do df
        max_accuracy = maximum(df[!, :accuracy])
        return subset(df, :accuracy => ByRow(==(max_accuracy)))
    end |>
        x -> vcat(x...; cols = :union)
end
