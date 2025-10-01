"""
    save_checkpoint_structvector(results, checkpoint_dir, batch_idx)

Save checkpoint file atomically for StructVector results.
"""
function save_checkpoint_structvector(
        results::StructVector{OptimizationResult},
        checkpoint_dir::String,
        batch_idx::Int
    )
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    checkpoint_file = joinpath(checkpoint_dir, "checkpoint_batch_$(batch_idx).jld2")
    temp_file = checkpoint_file * ".tmp"

    JLD2.jldsave(temp_file; results = results, batch_idx = batch_idx)
    return mv(temp_file, checkpoint_file; force = true)
end
