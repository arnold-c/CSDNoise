export load_previous_multistart_results_structvector,
    load_previous_multistart_results_df,
    load_previous_multistart_results

"""
    load_previous_multistart_results_structvector(filedir, filename_base)

Load previous multistart optimization results as StructVector, including checkpoints.
Returns empty StructVector if no results found.
"""
function load_previous_multistart_results_structvector(filedir::String, filename_base::String)
    if !isdir(filedir)
        return StructVector(OptimizationResult[])
    end

    # Get most recent file
    load_filepath = get_most_recent_hyperparam_filepath(filename_base, filedir)

    if Try.iserr(load_filepath)
        # Try to load checkpoint files
        return load_checkpoint_results_structvector(filedir)
    end

    try
        data = JLD2.load(Try.unwrap(load_filepath))

        # Try to load as StructVector first (new format)
        if haskey(data, "multistart_optimization_results")
            return data["multistart_optimization_results"]::StructVector{OptimizationResult}
        end

        # Fall back to DataFrame format (old format) and convert
        if haskey(data, "multistart_optimization_df")
            df = data["multistart_optimization_df"]
            return df_to_structvector(df, OptimizationResult)
        end

        return StructVector(OptimizationResult[])
    catch
        # Try to load checkpoint files as fallback
        return load_checkpoint_results_structvector(filedir)
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
