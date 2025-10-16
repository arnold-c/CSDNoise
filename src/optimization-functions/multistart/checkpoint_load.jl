export load_checkpoint_results_structvector,
    load_checkpoint_results

"""
    load_checkpoint_results_structvector(filedir)

Load and merge all checkpoint files in the directory as StructVector.
"""
function load_checkpoint_results_structvector(filedir::String)
    checkpoint_dir = joinpath(filedir, "checkpoints")

    if !isdir(checkpoint_dir)
        return StructVector(OptimizationResult[])
    end

    checkpoint_files = filter(
        f -> endswith(f, ".jld2") && startswith(f, "checkpoint_"),
        readdir(checkpoint_dir)
    )

    if isempty(checkpoint_files)
        return StructVector(OptimizationResult[])
    end

    # Sort checkpoint files by modification time (most recent first)
    sorted_files = sort(
        checkpoint_files,
        by = f -> mtime(joinpath(checkpoint_dir, f)),
        rev = true
    )

    all_results = OptimizationResult[]
    loaded_file = nothing

    for file in sorted_files
        filepath = joinpath(checkpoint_dir, file)
        try
            data = JLD2.load(filepath)
            if haskey(data, "results")
                results = data["results"]
                if results isa StructVector{OptimizationResult}
                    append!(all_results, results)
                elseif results isa DataFrame
                    # Convert DataFrame to StructVector
                    sv = df_to_structvector(results, OptimizationResult)
                    append!(all_results, sv)
                end
                if isnothing(loaded_file)
                    loaded_file = file
                end
            end
        catch e
            if isnothing(loaded_file)
                @warn "Failed to load checkpoint file $file: $e. Trying next checkpoint..."
            else
                @warn "Failed to load checkpoint file $file: $e"
            end
        end
    end

    if isempty(all_results)
        error("All checkpoint files failed to load")
    end

    if loaded_file != sorted_files[1]
        @info "Loaded checkpoint from $loaded_file (most recent checkpoint failed)"
    end

    return StructVector(all_results)
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

    # Sort checkpoint files by modification time (most recent first)
    sorted_files = sort(
        checkpoint_files,
        by = f -> mtime(joinpath(checkpoint_dir, f)),
        rev = true
    )

    # Load all checkpoint files and merge
    all_results = DF.DataFrame[]
    loaded_file = nothing

    for file in sorted_files
        filepath = joinpath(checkpoint_dir, file)
        try
            data = JLD2.load(filepath)
            if haskey(data, "partial_results")
                push!(all_results, data["partial_results"])
                if isnothing(loaded_file)
                    loaded_file = file
                end
            elseif haskey(data, "multistart_optimization_df")
                push!(all_results, data["multistart_optimization_df"])
                if isnothing(loaded_file)
                    loaded_file = file
                end
            end
        catch e
            if isnothing(loaded_file)
                @warn "Failed to load checkpoint file $file: $e. Trying next checkpoint..."
            else
                @warn "Failed to load checkpoint file $file: $e"
            end
        end
    end

    if isempty(all_results)
        error("All checkpoint files failed to load")
    end

    if loaded_file != sorted_files[1]
        @info "Loaded checkpoint from $loaded_file (most recent checkpoint failed)"
    end

    # Merge all checkpoint results
    merged_df = all_results[1]
    for i in 2:length(all_results)
        merged_df = merge_results_safely(merged_df, all_results[i])
    end

    return merged_df
end
