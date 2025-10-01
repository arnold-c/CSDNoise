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

    all_results = OptimizationResult[]

    for file in checkpoint_files
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
            end
        catch e
            @warn "Failed to load checkpoint file $file: $e"
        end
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

    # Load all checkpoint files and merge
    all_results = DF.DataFrame[]

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
