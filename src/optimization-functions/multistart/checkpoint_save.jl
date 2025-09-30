using DataFrames: DataFrame
using DrWatson: @tagsave

export save_checkpoint_atomic

"""
    save_checkpoint_atomic(df, checkpoint_dir, batch_idx)

Save checkpoint file atomically to prevent corruption.
"""
function save_checkpoint_atomic(df::DataFrame, checkpoint_dir::String, batch_idx::Int)
    if !isdir(checkpoint_dir)
        mkpath(checkpoint_dir)
    end

    # Create unique temporary filename with .jld2 extension
    # Use random number instead of threadid() to avoid relying on thread implementation details
    unique_id = rand(UInt32)
    temp_file = joinpath(
        checkpoint_dir,
        "checkpoint_$(getpid())_$(unique_id)_$(time_ns()).tmp.jld2"
    )

    # Save to temporary file
    return try
        @tagsave(temp_file, Dict("partial_results" => df))

        # Atomic rename to final filename
        final_file = joinpath(
            checkpoint_dir,
            "checkpoint_batch_$(batch_idx)_$(Dates.now()).jld2"
        )

        mv(temp_file, final_file; force = false)
        @info "Saved checkpoint: $final_file"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_file) && rm(temp_file; force = true)
        @warn "Failed to save checkpoint: $e"
    end
end
