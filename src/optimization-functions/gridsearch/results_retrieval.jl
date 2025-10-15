export load_previous_gridsearch_results_structvector

"""
    load_previous_gridsearch_results_structvector(filedir, filename_base)

Load previous grid search results as StructVector, including checkpoints.
Reuses checkpoint loading logic from multistart.
"""
function load_previous_gridsearch_results_structvector(filedir::String, filename_base::String)
    if !isdir(filedir)
        return StructVector(OptimizationResult[])
    end

    # Get most recent file - reuse function from multistart
    load_filepath = get_most_recent_hyperparam_filepath(filename_base, filedir)

    if Try.iserr(load_filepath)
        # Try to load checkpoint files
        return load_checkpoint_results_structvector(filedir)
    end

    try
        data = JLD2.load(Try.unwrap(load_filepath))

        # Try to load as StructVector first (new format)
        if haskey(data, "optimization_results")
            return data["optimization_results"]
        end

        # Fall back to DataFrame format (old format) and convert
        if haskey(data, "ews_df")
            df = data["ews_df"]
            return df_to_structvector(df, OptimizationResult)
        end

        return StructVector(OptimizationResult[])
    catch
        # Try to load checkpoint files as fallback
        return load_checkpoint_results_structvector(filedir)
    end
end
