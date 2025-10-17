export load_previous_gridsearch_results_structvector

"""
    load_previous_gridsearch_results_structvector(filedir, filename_base)

Load previous grid search results as StructVector, including checkpoints.
Reuses checkpoint loading logic from multistart.
"""
function load_previous_gridsearch_results_structvector(filedir::String, filename_base::String)
    if !isdir(filedir)
        error("The provided file directory $filedir isn't a valid directory. Check `filedir` again.")
    end

    # Get most recent file - reuse function from multistart
    load_filepath = get_most_recent_hyperparam_filepath(filename_base, filedir)

    checkpoint_warning_message = "Failed to load previous completed results so attempting to load the most recent checkpoint file"

    if Try.iserr(load_filepath)
        # Try to load checkpoint files
        @warn checkpoint_warning_message
        return load_checkpoint_results_structvector(filedir)
    end

    try
        data = JLD2.load(Try.unwrap(load_filepath))

        # Try to load as StructVector first (new format)
        if haskey(data, "optimization_results")
            return data["optimization_results"]
        end

        @warn "Failed to load the previous results in $load_filepath. Returning an empty StructVector to force the recreation of all scenarios."
        return StructVector(OptimizationResult[])
    catch
        # Try to load checkpoint files as fallback
        @warn checkpoint_warning_message
        return load_checkpoint_results_structvector(filedir)
    end
end
