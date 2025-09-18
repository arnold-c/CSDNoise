"""
    save_results_structvector(results, filepath)

Save StructVector results directly to JLD2 file.
"""
function save_results_structvector(results::StructVector{OptimizationResult}, filepath::String)
    # Ensure directory exists
    dir = dirname(filepath)
    !isdir(dir) && mkpath(dir)

    # Create temporary file with .jld2 extension
    temp_filepath = filepath * ".tmp.jld2"

    return try
        # Save StructVector directly to JLD2
        jldsave(temp_filepath; multistart_optimization_results = results)

        # Atomic rename
        mv(temp_filepath, filepath; force = true)
        @info "Saved optimization results to $filepath"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_filepath) && rm(temp_filepath; force = true)
        rethrow(e)
    end
end

"""
    save_checkpoint_structvector(results, checkpoint_dir, batch_idx)

Save checkpoint file atomically for StructVector results.
"""
function save_checkpoint_structvector(results::StructVector{OptimizationResult}, checkpoint_dir::String, batch_idx::Int)
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    checkpoint_file = joinpath(checkpoint_dir, "checkpoint_batch_$(batch_idx).jld2")
    temp_file = checkpoint_file * ".tmp"

    jldsave(temp_file; results = results, batch_idx = batch_idx)
    return mv(temp_file, checkpoint_file; force = true)
end
