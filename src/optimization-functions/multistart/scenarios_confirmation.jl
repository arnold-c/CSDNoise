export confirm_optimization_run_structvector,
    estimate_optimization_time,
    confirm_optimization_run

"""
    confirm_optimization_run_structvector(missing_scenarios, n_sobol_points; disable_time_check)

Confirm with user before running optimization on StructVector of scenarios.
"""
function confirm_optimization_run_structvector(
        missing_scenarios::StructVector{OptimizationScenario},
        n_sobol_points::Int;
        disable_time_check::Bool = false
    )
    n_missing = length(missing_scenarios)

    if n_missing == 0
        @info "No missing scenarios to optimize"
        return false
    end

    if disable_time_check
        return true
    end

    # Estimate time (convert to DataFrame temporarily for existing function)
    missing_df = DF.DataFrame(missing_scenarios)
    estimated_time = estimate_optimization_time(missing_df, n_sobol_points)

    if estimated_time > 300  # 5 minutes
        println(StyledStrings.styled"{yellow:Warning:} Estimated optimization time: {red:$(round(estimated_time/60, digits=1))} minutes")
        print("Continue? (y/N): ")
        response = readline()
        return lowercase(strip(response)) in ["y", "yes"]
    end

    return true
end


"""
    estimate_optimization_time(missing_scenarios_df, n_sobol_points)

Estimate total time needed for optimization based on number of scenarios and complexity.
"""
function estimate_optimization_time(
        missing_scenarios_df::DF.DataFrames.DataFrame,
        n_sobol_points::Int
    )
    n_missing = DF.nrow(missing_scenarios_df)

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
        missing_scenarios_df::DF.DataFrame,
        n_sobol_points::Int;
        disable_time_check::Bool = false
    )
    n_missing = DF.nrow(missing_scenarios_df)

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

    println(StyledStrings.styled"{yellow:Found $(n_missing) missing scenarios to optimize}")
    println(StyledStrings.styled"Estimated time: {blue:$(time_message)}")
    println(StyledStrings.styled"Sobol points per scenario: {cyan:$(n_sobol_points)}")

    # Ask for confirmation using REPL menu
    choice = RTM.request(
        "Do you want to continue with the optimization?",
        RTM.RadioMenu(["No", "Yes"]; ctrl_c_interrupt = false)
    )

    return choice == 2
end
