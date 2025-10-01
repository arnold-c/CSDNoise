export confirm_gridsearch_run_structvector

"""
    confirm_gridsearch_run_structvector(missing_scenarios; disable_time_check)

Confirm with user before running grid search on StructVector of scenarios.
"""
function confirm_gridsearch_run_structvector(
        missing_scenarios::StructVector{GridSearchScenario};
        disable_time_check::Bool = false
    )
    n_missing = length(missing_scenarios)

    if n_missing == 0
        @info "No missing grid points to evaluate"
        return false
    end

    if disable_time_check
        return true
    end

    # Estimate time (simpler than multistart since no optimization)
    estimated_time = n_missing * 0.5  # ~0.5 seconds per grid point

    if estimated_time > 300  # 5 minutes
        println(StyledStrings.styled"{yellow:Warning:} Estimated grid search time: {red:$(round(estimated_time/60, digits=1))} minutes")
        print("Continue? (y/N): ")
        response = readline()
        return lowercase(strip(response)) in ["y", "yes"]
    end

    return true
end
