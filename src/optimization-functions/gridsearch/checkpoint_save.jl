"""
    save_checkpoint_structvector(results, checkpoint_dir, batch_idx)

Save checkpoint file atomically for StructVector results.
"""
function save_checkpoint_structvector(
        results::StructVector{OptimizationResult},
        checkpoint_dir::String,
        checkpoint_output_filename_base::String,
        batch_idx::Int64
    )
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    checkpoint_file = joinpath(checkpoint_dir, "$(checkpoint_output_filename_base)$(batch_idx).jld2")
    temp_file = checkpoint_file * ".tmp"

    # save to temp and then move to force an overwrite if back with the same name
    # already exists (though it shouldn't as it uses the datetime in the file name)
    JLD2.jldsave(temp_file; optimization_results = results, batch_idx = batch_idx)
    return mv(temp_file, checkpoint_file; force = true)
end
