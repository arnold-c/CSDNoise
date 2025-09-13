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
using REPL.TerminalMenus: RadioMenu, request

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

struct OptimizedValues
    threshold_percentile::Float64
    consecutive_thresholds::Int64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
end


"""
    CachedSimulationData

Pre-computed simulation data that can be reused across parameter evaluations.
This avoids expensive recomputation of noise arrays and test arrays.
"""
struct CachedSimulationData
    ensemble_single_incarr::Array{Int64, 3}
    testarr::Array{Int64, 3}
    null_testarr::Array{Int64, 3}
    thresholds::Vector{Matrix{Int64}}
    ews_metrics::Vector{EWSMetrics}
    null_ews_metrics::Vector{EWSMetrics}
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
        # Caching options
        batch_size = 10,
        save_every_n = 10,
        disable_time_check = false,
    )
    if !isdir(filedir)
        mkpath(filedir)
    end

    # Setup checkpoint directory
    checkpoint_dir = joinpath(filedir, "checkpoints")

    # Create all scenario combinations as DataFrame
    all_scenarios_df = create_scenarios_dataframe(specification_vecs)
    n_total_scenarios = nrow(all_scenarios_df)

    if verbose
        println(styled"{green:Starting Multistart Optimization with Caching}")
        println(styled"Total scenarios defined: {yellow:$n_total_scenarios}")
        println(styled"Sobol points per scenario: {blue:$n_sobol_points}")
    end

    # Load existing results (including checkpoints)
    existing_results_df = if force
        create_empty_results_dataframe()
    else
        load_previous_multistart_results_df(filedir, optimization_filename_base)
    end

    n_existing = nrow(existing_results_df)
    if verbose && n_existing > 0
        println(styled"Found {cyan:$n_existing} existing results")
    end

    # Find missing scenarios
    missing_scenarios_df = find_missing_scenarios(all_scenarios_df, existing_results_df)
    n_missing = nrow(missing_scenarios_df)

    if verbose
        println(styled"Missing scenarios to optimize: {yellow:$n_missing}")
    end

    # Check with user if needed
    if n_missing > 0
        if !confirm_optimization_run(
                missing_scenarios_df,
                n_sobol_points;
                disable_time_check = disable_time_check
            )
            @info "Optimization cancelled by user"
            return return_df ? existing_results_df : nothing
        end
    else
        @info "All scenarios already optimized"
        return return_df ? existing_results_df : nothing
    end

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

    # Run optimization for missing scenarios in batches
    new_results_df = optimize_scenarios_in_batches(
        missing_scenarios_df,
        data_arrs,
        bounds,
        optim_config;
        batch_size = batch_size,
        executor = executor,
        save_every_n = save_every_n,
        save_results = save_results,
        checkpoint_dir = checkpoint_dir,
        verbose = verbose
    )

    # Combine existing and new results
    final_results_df = if nrow(existing_results_df) > 0
        merge_results_safely(existing_results_df, new_results_df)
    else
        new_results_df
    end

    # Save final results
    if save_results
        save_results_atomically(final_results_df, optimization_output_filepath)

        # Clean up checkpoints after successful save
        cleanup_checkpoints(checkpoint_dir)
    end

    if verbose
        n_final = nrow(final_results_df)
        println(styled"{green:Optimization complete! Total results: {yellow:$n_final}}")
    end

    return return_df ? final_results_df : nothing
end

"""
    optimize_single_scenario(scenario, data_arrs, bounds, config)

Optimize EWS parameters for a single scenario using multistart optimization.
"""
function optimize_single_scenario(
        scenario::OptimizationScenario,
        data_arrs::T1,
        bounds::T2,
        config::T3
    ) where {T1 <: NamedTuple, T2 <: NamedTuple, T3 <: NamedTuple}
    @unpack ews_metric = scenario

    # Pre-compute expensive simulation data once per scenario
    cached_data = create_cached_simulation_data(scenario, data_arrs)

    # Create tracker instance for this scenario
    tracker = OptimizationTracker()

    # Create objective function closure that updates tracker
    objective = params -> ews_objective_function_with_tracking(
        params,
        scenario,
        cached_data,
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

    return OptimizedValues(
        optimal_params.threshold_percentile,
        optimal_params.consecutive_thresholds,
        tracker.best_accuracy,
        tracker.best_sensitivity,
        tracker.best_specificity,
    )
end

"""
    ews_objective_function_with_tracking(params_vec, scenario, cached_data, tracker)

Objective function for EWS parameter optimization that tracks the best solution.
Returns 1 - accuracy for minimization and updates tracker with best metrics.
Uses pre-computed cached simulation data for efficiency.
"""
function ews_objective_function_with_tracking(
        params_vec::Vector{Float64},
        scenario::OptimizationScenario,
        cached_data::CachedSimulationData,
        tracker::OptimizationTracker
    )
    @unpack ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window, ews_threshold_burnin, ews_metric = scenario

    @unpack ensemble_single_incarr,
        testarr,
        null_testarr,
        thresholds,
        ews_metrics,
        null_ews_metrics = cached_data

    # Map continuous parameters to EWS parameters
    # TODO: Confirm if it's possible to optimize with integer values
    ews_params = map_continuous_to_ews_parameters(params_vec)

    # Pre-compute values outside the loop for type stability
    ews_metric_symbol = Symbol(ews_metric)
    threshold_percentile = ews_params.threshold_percentile

    # Calculate accuracy
    valid_nsims = length(ews_metrics)
    true_positives = 0
    true_negatives = 0

    for sim in eachindex(ews_metrics)
        # Use pre-computed EWS metrics
        ews_vals = ews_metrics[sim]
        null_ews_vals = null_ews_metrics[sim]

        # Check threshold exceedances
        exceeds_threshold = expanding_ews_thresholds(
            ews_vals,
            ews_metric_symbol,
            ews_threshold_window,
            threshold_percentile,
            ews_threshold_burnin,
            # percentiles = threshold_percentile,
            # burn_in = ews_threshold_burnin,
        )

        detection_index = calculate_ews_trigger_index(
            exceeds_threshold;
            consecutive_thresholds = ews_params.consecutive_thresholds,
        )

        null_exceeds_threshold = expanding_ews_thresholds(
            null_ews_vals,
            ews_metric_symbol,
            ews_threshold_window,
            threshold_percentile,
            ews_threshold_burnin,
            # percentiles = threshold_percentile,
            # burn_in = ews_threshold_burnin,
        )

        null_detection_index = calculate_ews_trigger_index(
            null_exceeds_threshold;
            consecutive_thresholds = ews_params.consecutive_thresholds,
        )

        # Update counts
        if Try.isok(detection_index)
            true_positives += 1
        end
        if Try.iserr(null_detection_index)
            true_negatives += 1
        end
    end

    # Calculate metrics
    sensitivity = calculate_sensitivity(true_positives, valid_nsims)
    specificity = calculate_specificity(true_negatives, valid_nsims)
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

"""
    calculate_ews_metrics_for_simulation(ews_metric_specification, testarr_view, null_testarr_view, threshold, ews_enddate_type)

Calculate EWS metrics for a single simulation given the test array views and threshold.
Returns a tuple of (ews_metrics, null_ews_metrics) wrapped in Try types.
"""
function calculate_ews_metrics_for_simulation(
        ews_metric_specification,
        testarr_view,
        null_testarr_view,
        threshold,
        ews_enddate_type
    )
    enddate = calculate_ews_enddate(threshold, ews_enddate_type)

    if Try.isok(enddate)
        enddate_val = Try.unwrap(enddate)

        # Calculate EWS metrics
        ews_vals = EWSMetrics(
            ews_metric_specification,
            @view(testarr_view[1:enddate_val, 5])
        )

        null_ews_vals = EWSMetrics(
            ews_metric_specification,
            @view(null_testarr_view[1:enddate_val, 5])
        )

        return Try.Ok((ews_vals, null_ews_vals))
    else
        return enddate  # Both are Try.Err
    end
end

"""
    create_cached_simulation_data(scenario, data_arrs)

Pre-compute expensive simulation data once per scenario to avoid repeated computation
during parameter optimization. This includes noise arrays, test arrays, and EWS metrics.
"""
function create_cached_simulation_data(
        scenario::OptimizationScenario,
        data_arrs::NamedTuple
    )
    @unpack noise_specification, test_specification, percent_tested,
        ews_enddate_type, ews_metric_specification = scenario

    @unpack ensemble_single_incarr, null_single_incarr,
        ensemble_specification, ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs = data_arrs

    # Create noise array once per scenario (expensive operation)
    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr,
        ensemble_specification;
        seed = 1234,
    )[1]

    # Create test arrays once per scenario (expensive operation)
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
    thresholds = get_enddate_thresholds(ews_enddate_type, data_arrs)

    # Pre-compute EWS metrics for all simulations
    ensemble_nsims = size(ensemble_single_incarr, 3)
    ews_metrics = Vector{EWSMetrics}()
    null_ews_metrics = Vector{EWSMetrics}()

    # Keep track of valid simulation indices
    valid_sim_indices = Vector{Int}()

    for sim in 1:ensemble_nsims
        ews_metric_result = calculate_ews_metrics_for_simulation(
            ews_metric_specification,
            @view(testarr[:, :, sim]),
            @view(null_testarr[:, :, sim]),
            thresholds[sim],
            ews_enddate_type
        )
        if Try.isok(ews_metric_result)
            ews_result, null_ews_result = Try.unwrap(ews_metric_result)

            push!(ews_metrics, ews_result)
            push!(null_ews_metrics, null_ews_result)
            push!(valid_sim_indices, sim)
        end
    end

    # Ensure we have at least some valid simulations
    if isempty(ews_metrics)
        error("No valid EWS metrics could be computed for any simulation")
    end

    return CachedSimulationData(
        ensemble_single_incarr,
        testarr,
        null_testarr,
        thresholds,
        ews_metrics,
        null_ews_metrics
    )
end

"""
    map_continuous_to_ews_parameters(params_vec)

Map continuous optimization parameters to discrete EWS parameters.
"""
function map_continuous_to_ews_parameters(params_vec::Vector{Float64})
    return (
        threshold_percentile = params_vec[1],  # Already in correct range
        consecutive_thresholds = round(Int, params_vec[2]),
    )
end

"""
    create_scenarios_dataframe(specification_vecs)

Create a DataFrame containing all scenario combinations for optimization.
This replaces the nested loop approach with a DataFrame-based approach for better caching.
"""
function create_scenarios_dataframe(specification_vecs)
    scenarios = create_optimization_scenarios(specification_vecs) |> StructVector

    # Convert StructArray to DataFrame
    return DataFrame(scenarios)
end

"""
    create_optimization_scenarios(specification_vecs)

Create all scenario combinations for optimization.
Validates specification vectors before creating combinations.
"""
@unstable function create_optimization_scenarios(specification_vecs)
    # Validate specification vectors before processing
    _validate_specification_vectors(specification_vecs)

    @unpack noise_specification_vec, test_specification_vec,
        percent_tested_vec, ews_metric_specification_vec,
        ews_enddate_type_vec, ews_threshold_window_vec,
        ews_threshold_burnin_vec, ews_metric_vec = specification_vecs

    # Use Iterators.product instead of nested loops
    combinations = Iterators.product(
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec
    ) |> collect

    scenarios_vec = mapreduce(vcat, combinations; init = OptimizationScenario[]) do (
                noise_spec,
                test_spec,
                percent_tested,
                ews_metric_spec,
                ews_enddate_type,
                ews_window,
                ews_burnin,
                ews_metric,
            )
        OptimizationScenario(
            noise_spec,
            test_spec,
            percent_tested,
            ews_metric_spec,
            ews_enddate_type,
            ews_window,
            ews_burnin,
            ews_metric,
        )
    end

    # Convert to StructArray for easier manipulation
    scenarios = StructVector(scenarios_vec)

    return scenarios
end

"""
    validate_specification_vectors(specification_vecs)

Validate that specification vectors contain valid values and types.
Throws ArgumentError for invalid specifications and warns for empty vectors.
"""
function _validate_specification_vectors(specification_vecs)
    # Define required fields
    required_fields = Set(
        [
            :noise_specification_vec,
            :test_specification_vec,
            :percent_tested_vec,
            :ews_metric_specification_vec,
            :ews_enddate_type_vec,
            :ews_threshold_window_vec,
            :ews_threshold_burnin_vec,
            :ews_metric_vec,
        ]
    )

    # Get provided fields
    provided_fields = Set(keys(specification_vecs))

    # Check for missing required fields
    missing_fields = setdiff(required_fields, provided_fields)
    if !isempty(missing_fields)
        throw(ArgumentError("Missing required fields: $(join(sort(collect(missing_fields)), ", "))"))
    end

    # Check for extra fields
    extra_fields = setdiff(provided_fields, required_fields)
    if !isempty(extra_fields)
        throw(ArgumentError("Unexpected extra fields provided: $(join(sort(collect(extra_fields)), ", ")). Only these fields are allowed: $(join(sort(collect(required_fields)), ", "))"))
    end

    @unpack noise_specification_vec, test_specification_vec,
        percent_tested_vec, ews_metric_specification_vec,
        ews_enddate_type_vec, ews_threshold_window_vec,
        ews_threshold_burnin_vec, ews_metric_vec = specification_vecs

    # Validate noise specifications
    if !all(spec -> spec isa NoiseSpecification, noise_specification_vec)
        throw(ArgumentError("All noise specifications must be of type NoiseSpecification"))
    end

    # Validate test specifications
    if !all(spec -> spec isa IndividualTestSpecification, test_specification_vec)
        throw(ArgumentError("All test specifications must be of type IndividualTestSpecification"))
    end

    # Validate percent tested values
    if !all(p -> 0.0 <= p <= 1.0, percent_tested_vec)
        throw(ArgumentError("All percent_tested values must be between 0.0 and 1.0"))
    end

    # Validate EWS metric specifications
    # if !all(spec -> spec isa EWSMetricSpecification, ews_metric_specification_vec)
    #     throw(ArgumentError("All EWS metric specifications must be of type EWSMetricSpecification"))
    # end

    # Validate EWS end date types
    if !all(enddate -> enddate isa EWSEndDateType, ews_enddate_type_vec)
        throw(ArgumentError("All EWS end date types must be of type EWSEndDateType"))
    end

    # Validate threshold window types
    if !all(window -> window isa EWSThresholdWindowType, ews_threshold_window_vec)
        throw(ArgumentError("All threshold windows must be either ExpandingThresholdWindow or RollingThresholdWindow"))
    end

    # Validate burnin periods
    if !all(burnin -> burnin isa Union{Dates.Day, Dates.Year}, ews_threshold_burnin_vec)
        throw(ArgumentError("All burnin periods must be of type Day or Year"))
    end

    # Validate EWS metric names
    valid_metrics = [
        "autocorrelation", "autocovariance", "coefficient_of_variation",
        "index_of_dispersion", "kurtosis", "mean", "skewness", "variance",
    ]
    if !all(metric -> metric in valid_metrics, ews_metric_vec)
        invalid_metrics = setdiff(ews_metric_vec, valid_metrics)
        throw(ArgumentError("Invalid EWS metrics: $(invalid_metrics). Valid metrics are: $(valid_metrics)"))
    end

    # Check for empty vectors (warn but don't error)
    empty_fields = String[]
    for field in required_fields
        if isempty(getfield(specification_vecs, field))
            push!(empty_fields, string(field))
        end
    end

    if !isempty(empty_fields)
        @warn "Empty specification vectors detected: $(join(empty_fields, ", ")). This will result in zero scenarios."
    end

    return nothing
end


"""
    find_missing_scenarios(all_scenarios_df, completed_results_df)

Find scenarios that haven't been computed yet using antijoin.
Returns DataFrame of missing scenarios.
"""
function find_missing_scenarios(all_scenarios_df::DataFrame, completed_results_df::DataFrame)
    if nrow(completed_results_df) == 0
        return all_scenarios_df
    end

    # Define scenario columns for comparison
    scenario_cols = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
    ]

    # Find scenarios that are in all_scenarios but not in completed_results
    missing_scenarios_df = antijoin(
        all_scenarios_df,
        completed_results_df;
        on = scenario_cols
    )

    return missing_scenarios_df
end

"""
    estimate_optimization_time(missing_scenarios_df, n_sobol_points)

Estimate total time needed for optimization based on number of scenarios and complexity.
"""
function estimate_optimization_time(missing_scenarios_df::DataFrame, n_sobol_points::Int)
    n_missing = nrow(missing_scenarios_df)

    if n_missing == 0
        return 0.0
    end

    # Base time estimate per scenario (in seconds)
    # This is a rough estimate - adjust based on your system performance
    base_time_per_scenario = 30.0  # seconds

    # Scale by number of Sobol points (more points = more time)
    sobol_scaling = n_sobol_points / 100.0  # Normalize to 100 points baseline

    # Scale by scenario complexity (some metrics/configurations take longer)
    complexity_scaling = 1.0  # Could be refined based on specific parameters

    total_time_s = n_missing * base_time_per_scenario * sobol_scaling * complexity_scaling

    return total_time_s
end

"""
    confirm_optimization_run(missing_scenarios_df, n_sobol_points; disable_time_check=false)

Ask user for confirmation before running optimization, showing time estimate.
"""
function confirm_optimization_run(
        missing_scenarios_df::DataFrame,
        n_sobol_points::Int;
        disable_time_check::Bool = false
    )
    n_missing = nrow(missing_scenarios_df)

    if n_missing == 0
        @info "No missing scenarios to optimize"
        return false
    end

    if disable_time_check
        return true
    end

    # Estimate time
    total_time_s = estimate_optimization_time(missing_scenarios_df, n_sobol_points)

    # Format time message
    time_message = if total_time_s < 60
        "approximately $(round(total_time_s; digits = 0)) seconds"
    elseif total_time_s < 3600
        minutes = round(total_time_s / 60; digits = 1)
        "approximately $(minutes) minutes"
    else
        hours = floor(total_time_s / 3600)
        minutes = round((total_time_s % 3600) / 60; digits = 0)
        "approximately $(hours) hours and $(minutes) minutes"
    end

    println(styled"{yellow:Found $(n_missing) missing scenarios to optimize}")
    println(styled"Estimated time: {blue:$(time_message)}")
    println(styled"Sobol points per scenario: {cyan:$(n_sobol_points)}")

    # Ask for confirmation using REPL menu
    choice = request(
        "Do you want to continue with the optimization?",
        RadioMenu(["No", "Yes"]; ctrl_c_interrupt = false)
    )

    return choice == 2
end

"""
    optimize_scenarios_in_batches(
        missing_scenarios_df, data_arrs, bounds, optim_config;
        batch_size=10, executor=FLoops.SequentialEx(), save_every_n=5, checkpoint_dir=""
    )

Optimize scenarios in batches with thread safety and incremental saving.

# Thread Safety Design
This function follows a fork-join pattern with clear phase separation:

1. **Parallel Phase**: Each batch is processed with FLoops.@floop where:
   - Each task writes to a unique, pre-allocated index in batch_results[idx]
   - No shared mutable state between tasks
   - No use of threadid() or thread-local storage

2. **Serial Phase**: Between batches, all I/O and state management happens single-threaded:
   - DataFrame creation and merging
   - Progress updates
   - Checkpoint saving
   - Result accumulation

This avoids the threadid() anti-pattern described in:
https://julialang.org/blog/2023/07/PSA-dont-use-threadid/

The design is robust to task migration, thread adoption, and other Julia threading
implementation details because it reasons about tasks, not threads.
"""
function optimize_scenarios_in_batches(
        missing_scenarios_df::DataFrame,
        data_arrs::T1,
        bounds::T2,
        optim_config::T3;
        batch_size::Int = 10,
        executor = FLoops.SequentialEx(),
        save_every_n::Int = 10,
        save_results = true,
        checkpoint_dir::String = "",
        verbose::Bool = true
    ) where {T1 <: NamedTuple, T2 <: NamedTuple, T3 <: NamedTuple}
    n_missing = nrow(missing_scenarios_df)

    if n_missing == 0
        return DataFrame()
    end

    # Auto-configure batch size for threaded execution
    # Note: We use nthreads() only as a heuristic, not for correctness
    if executor isa FLoops.ThreadedEx
        if batch_size == 10  # Default value
            # Use nthreads() as a heuristic for load balancing, but cap at reasonable bounds
            suggested_batch_size = max(1, n_missing ÷ (Threads.nthreads() * 4))
            batch_size = clamp(suggested_batch_size, 1, 100)  # Cap between 1 and 100
            verbose && @info "Auto-configured batch size for threading: $batch_size (suggested: $suggested_batch_size)"
        end
    end

    # Initialize results storage
    all_results = DataFrame[]

    # Setup progress tracking
    if verbose
        prog = Progress(n_missing; desc = "Optimizing scenarios: ", showspeed = true)
    end

    # Process scenarios in batches
    scenario_batches = collect(Iterators.partition(1:n_missing, batch_size))

    for (batch_idx, batch_indices) in pairs(scenario_batches)
        batch_scenarios = missing_scenarios_df[batch_indices, :]
        batch_size_actual = length(batch_indices)

        verbose && println(styled"{green:Processing batch $batch_idx/$(length(scenario_batches)) ($(batch_size_actual) scenarios)}")

        # Process batch in parallel (thread-safe within batch)
        # THREAD SAFETY: Pre-allocate fixed-size array where each task writes to unique index
        # This avoids the threadid() anti-pattern and follows task-based parallelism principles
        batch_results = Vector{OptimizedValues}(undef, batch_size_actual)

        FLoops.@floop executor for (idx, scenario_row) in pairs(eachrow(batch_scenarios))
            # Convert DataFrame row to NamedTuple
            scenario = dataframe_row_to_scenario(scenario_row)

            # Run optimization (each task works independently)
            result = optimize_single_scenario(
                scenario,
                data_arrs,
                bounds,
                optim_config
            )

            # Store result at unique index (thread-safe: no overlapping writes)
            batch_results[idx] = result
        end

        # SERIAL PHASE: Convert batch results to DataFrame (single-threaded)
        # All parallel tasks have completed at this point (implicit barrier)
        batch_scenarios_vec = StructVector([dataframe_row_to_scenario(row) for row in eachrow(batch_scenarios)])
        batch_df = create_results_dataframe(
            batch_scenarios_vec,
            StructVector(batch_results),
        )

        # Add to results (single-threaded - no race conditions)
        push!(all_results, batch_df)

        # Update progress (single-threaded)
        verbose && update!(prog, sum(length.(scenario_batches[1:batch_idx])))

        # Save checkpoint periodically (single-threaded I/O)
        if batch_idx % save_every_n == 0 && !isempty(checkpoint_dir) && save_results
            combined_df = vcat(all_results...; cols = :union)
            save_checkpoint_atomic(combined_df, checkpoint_dir, batch_idx)
        end
    end

    # Combine all results
    final_results_df = if !isempty(all_results)
        vcat(all_results...; cols = :union)
    else
        create_empty_results_dataframe()
    end

    return final_results_df
end

"""
    dataframe_row_to_scenario(row)

Convert a DataFrame row to an OptimizationScenario struct.
"""
function dataframe_row_to_scenario(row::DataFrameRow)
    return OptimizationScenario(
        row.noise_specification,
        row.test_specification,
        row.percent_tested,
        row.ews_metric_specification,
        row.ews_enddate_type,
        row.ews_threshold_window,
        row.ews_threshold_burnin,
        row.ews_metric,
    )
end

"""
    save_checkpoint_atomic(df, checkpoint_dir, batch_idx)

Save checkpoint file atomically to prevent corruption.
"""
function save_checkpoint_atomic(df::DataFrame, checkpoint_dir::String, batch_idx::Int)
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    # Create unique temporary filename with .jld2 extension
    # Use random number instead of threadid() to avoid relying on thread implementation details
    unique_id = rand(UInt32)
    temp_file = joinpath(
        checkpoint_dir,
        "checkpoint_$(getpid())_$(unique_id)_$(time_ns()).tmp.jld2"
    )

    # Save to temporary file
    return try
        @tagsave(temp_file, Dict("partial_results" => df))

        # Atomic rename to final filename
        final_file = joinpath(
            checkpoint_dir,
            "checkpoint_batch_$(batch_idx)_$(Dates.now()).jld2"
        )

        mv(temp_file, final_file; force = false)
        @info "Saved checkpoint: $final_file"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_file) && rm(temp_file; force = true)
        @warn "Failed to save checkpoint: $e"
    end
end

"""
    save_results_atomically(results_df, filepath)

Save final results atomically to prevent corruption during write.
"""
function save_results_atomically(results_df::DataFrame, filepath::String)
    # Ensure directory exists
    dir = dirname(filepath)
    !isdir(dir) && mkpath(dir)

    # Create temporary file with .jld2 extension
    temp_filepath = filepath * ".tmp.jld2"

    try
        # Save to temporary file
        @tagsave(temp_filepath, Dict("multistart_optimization_df" => results_df))

        # Atomic rename
        mv(temp_filepath, filepath; force = true)
        @info "Saved optimization results to $filepath"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_filepath) && rm(temp_filepath; force = true)
        rethrow(e)
    end

    return nothing
end

"""
    cleanup_checkpoints(checkpoint_dir)

Remove checkpoint files after successful completion.
"""
function cleanup_checkpoints(checkpoint_dir::String)
    if !isdir(checkpoint_dir)
        return
    end

    checkpoint_files = filter(
        f -> endswith(f, ".jld2") && startswith(f, "checkpoint_"),
        readdir(checkpoint_dir)
    )

    for file in checkpoint_files
        filepath = joinpath(checkpoint_dir, file)
        try
            rm(filepath)
        catch e
            @warn "Failed to remove checkpoint file $file: $e"
        end
    end

    # Remove checkpoint directory if empty
    try
        if isempty(readdir(checkpoint_dir))
            rm(checkpoint_dir)
        end
    catch
        # Directory not empty or other issue - leave it
    end

    return nothing
end

"""
    scenarios_equal(scenario1, scenario2)

Compare two scenarios for equality, handling floating point comparisons properly.
"""
function scenarios_equal(scenario1::OptimizationScenario, scenario2::OptimizationScenario)
    # Get all field names
    fields = fieldnames(OptimizationScenario)

    # Compare each field
    for field in fields
        val1 = getfield(scenario1, field)
        val2 = getfield(scenario2, field)

        # Handle floating point comparison
        if val1 isa AbstractFloat && val2 isa AbstractFloat
            if !isapprox(val1, val2; rtol = 1.0e-10)
                return false
            end
        else
            if val1 != val2
                return false
            end
        end
    end

    return true
end

"""
    validate_results_integrity(results_df)

Validate the integrity of results DataFrame.
"""
function validate_results_integrity(results_df::DataFrame)
    required_cols = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
        :ews_threshold_percentile,
        :ews_consecutive_thresholds,
        :accuracy,
        :sensitivity,
        :specificity,
    ]

    # Check required columns exist
    missing_cols = setdiff(required_cols, names(results_df))
    if !isempty(missing_cols)
        @warn "Missing required columns: $missing_cols"
        return false
    end

    # Check for missing values in critical columns
    for col in [:accuracy, :sensitivity, :specificity]
        if any(ismissing, results_df[!, col])
            @warn "Missing values found in column: $col"
            return false
        end
    end

    # Check accuracy bounds
    if any(x -> x < 0 || x > 1, results_df.accuracy)
        @warn "Accuracy values outside [0,1] range found"
        return false
    end

    @info "Results DataFrame validation passed"
    return true
end

"""
    summarize_optimization_results(results_df)

Provide a summary of optimization results for user feedback.
"""
function summarize_optimization_results(results_df::DataFrame)
    if nrow(results_df) == 0
        println(styled"{yellow:No results to summarize}")
        return
    end

    println(styled"{green:=== Optimization Results Summary ===}")
    println(styled"Total scenarios: {yellow:$(nrow(results_df))}")

    # Accuracy statistics
    acc_stats = (
        mean = mean(results_df.accuracy),
        median = median(results_df.accuracy),
        min = minimum(results_df.accuracy),
        max = maximum(results_df.accuracy),
    )

    println(
        styled"Accuracy - Mean: {blue:$(round(acc_stats.mean, digits=3))}, " *
            "Median: {blue:$(round(acc_stats.median, digits = 3))}, " *
            "Range: {blue:$(round(acc_stats.min, digits = 3))} - {blue:$(round(acc_stats.max, digits = 3))}"
    )

    # Count by major categories
    noise_counts = combine(groupby(results_df, :noise_specification), nrow => :count)
    println(styled"Noise specifications: {cyan:$(nrow(noise_counts))} types")

    metric_counts = combine(groupby(results_df, :ews_metric), nrow => :count)
    println(styled"EWS metrics: {cyan:$(nrow(metric_counts))} types")

    # Best performing scenarios
    best_idx = argmax(results_df.accuracy)
    best_accuracy = results_df.accuracy[best_idx]
    return println(
        styled"Best accuracy: {green:$(round(best_accuracy, digits=4))} " *
            "with metric: {yellow:$(results_df.ews_metric[best_idx])}"
    )
end

"""
    create_results_dataframe(results, scenarios)

Convert optimization results to DataFrame format.
"""
function create_results_dataframe(
        scenarios::StructVector{OptimizationScenario},
        results::StructVector{OptimizedValues},
    )
    @assert length(scenarios) == length(results)
    scenarios_df = DataFrame(scenarios)
    results_df = DataFrame(results)
    return hcat(scenarios_df, results_df)
    # df = DataFrame(
    #     noise_specification = [s.noise_specification for s in scenarios],
    #     test_specification = [s.test_specification for s in scenarios],
    #     percent_tested = [s.percent_tested for s in scenarios],
    #     ews_metric_specification = [s.ews_metric_specification for s in scenarios],
    #     ews_enddate_type = [s.ews_enddate_type for s in scenarios],
    #     ews_threshold_window = [s.ews_threshold_window for s in scenarios],
    #     ews_metric = [s.ews_metric for s in scenarios],
    #     # ews_threshold_burnin = [s.ews_threshold_burnin for s in scenarios], # NOTE: type instability as burnin <:Period
    #     ews_threshold_percentile = [r.optimal_threshold_percentile for r in results],
    #     ews_consecutive_thresholds = [r.optimal_consecutive_thresholds for r in results],
    #     accuracy = [r.accuracy for r in results],
    #     sensitivity = [r.sensitivity for r in results],
    #     specificity = [r.specificity for r in results],
    #     # n_evaluations = [r.n_evaluations for r in results],
    #     # convergence_status = [r.convergence_status for r in results],
    # )
    #
    # return df
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
    load_previous_multistart_results_df(filedir, filename_base)

Load previous multistart optimization results as DataFrame, including checkpoints.
Returns empty DataFrame if no results found.
"""
function load_previous_multistart_results_df(filedir::String, filename_base::String)
    # Try to load main results file
    main_results = load_previous_multistart_results(filedir, filename_base)

    # Try to load checkpoint files
    checkpoint_results = load_checkpoint_results(filedir)

    # Merge results
    if !isnothing(main_results) && !isnothing(checkpoint_results)
        return merge_results_safely(main_results, checkpoint_results)
    elseif !isnothing(main_results)
        return main_results
    elseif !isnothing(checkpoint_results)
        return checkpoint_results
    else
        return create_empty_results_dataframe()
    end
end

"""
    load_checkpoint_results(filedir)

Load and merge all checkpoint files in the directory.
"""
function load_checkpoint_results(filedir::String)
    checkpoint_dir = joinpath(filedir, "checkpoints")

    if !isdir(checkpoint_dir)
        return nothing
    end

    checkpoint_files = filter(
        f -> endswith(f, ".jld2") && startswith(f, "checkpoint_"),
        readdir(checkpoint_dir)
    )

    if isempty(checkpoint_files)
        return nothing
    end

    # Load all checkpoint files and merge
    all_results = DataFrame[]

    for file in checkpoint_files
        filepath = joinpath(checkpoint_dir, file)
        try
            data = JLD2.load(filepath)
            if haskey(data, "partial_results")
                push!(all_results, data["partial_results"])
            elseif haskey(data, "multistart_optimization_df")
                push!(all_results, data["multistart_optimization_df"])
            end
        catch e
            @warn "Failed to load checkpoint file $file: $e"
        end
    end

    if isempty(all_results)
        return nothing
    end

    # Merge all checkpoint results
    merged_df = all_results[1]
    for i in 2:length(all_results)
        merged_df = merge_results_safely(merged_df, all_results[i])
    end

    return merged_df
end

"""
    merge_results_safely(df1, df2)

Safely merge two results DataFrames, removing duplicates based on scenario parameters.
"""
function merge_results_safely(df1::DataFrame, df2::DataFrame)
    if nrow(df1) == 0
        return df2
    elseif nrow(df2) == 0
        return df1
    end

    # Define scenario columns for duplicate detection
    scenario_cols = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
    ]

    # Combine DataFrames
    combined_df = vcat(df1, df2; cols = :union)

    # Remove duplicates, keeping first occurrence
    # (assumes df1 contains more recent/authoritative results)
    unique_df = unique(combined_df, scenario_cols)

    return unique_df
end

"""
    create_empty_results_dataframe()

Create an empty DataFrame with the correct structure for multistart optimization results.
"""
function create_empty_results_dataframe()
    return DataFrame(
        noise_specification = NoiseSpecification[],
        test_specification = eltype(IndividualTestSpecification)[],
        percent_tested = Float64[],
        ews_metric_specification = EWSMetricSpecification[],
        ews_enddate_type = EWSEndDateType[],
        ews_threshold_window = Union{Type{ExpandingThresholdWindow}, Type{RollingThresholdWindow}}[],
        ews_metric = String[],
        ews_threshold_percentile = Float64[],
        ews_consecutive_thresholds = Int[],
        ews_threshold_burnin = Union{Dates.Day, Dates.Year}[],
        accuracy = Float64[],
        sensitivity = Float64[],
        specificity = Float64[]
    )
end
