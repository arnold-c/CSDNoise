export cleanup_checkpoints

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
