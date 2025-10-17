export load_checkpoint_results_structvector

"""
    load_checkpoint_results_structvector(filedir)

Load the most recent valid checkpoint file from the directory as StructVector.
Falls back to older checkpoints if the most recent one is corrupted.
"""
function load_checkpoint_results_structvector(filedir::String)
    checkpoint_dir = joinpath(filedir, "checkpoints")

    if !isdir(checkpoint_dir)
        return StructVector(OptimizationResult[])
    end

    checkpoint_files = filter(
        f -> endswith(f, ".jld2") && contains(f, "checkpoint_"),
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

    for file in sorted_files
        filepath = joinpath(checkpoint_dir, file)
        try
            data = JLD2.load(filepath)
            if haskey(data, "optimization_results")
                results = data["optimization_results"]
                if results isa StructVector{OptimizationResult}
                    if file != sorted_files[1]
                        @info "Loaded checkpoint from $file (most recent checkpoint failed)"
                    end
                    return results
                end
            end
        catch e
            @warn "Failed to load checkpoint file $file: $e. Trying next checkpoint..."
        end
    end

    error("All checkpoint files failed to load")
end
