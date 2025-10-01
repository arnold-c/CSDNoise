export save_results_atomically

"""
    save_results_atomically(results_df, filepath)

Save final results atomically to prevent corruption during write.
"""
function save_results_atomically(results_df::DF.DataFrame, filepath::String)
    # Ensure directory exists
    dir = dirname(filepath)
    !isdir(dir) && mkpath(dir)

    # Create temporary file with .jld2 extension
    temp_filepath = filepath * ".tmp.jld2"

    try
        # Save to temporary file
        DrWatson.@tagsave(temp_filepath, Dict("multistart_optimization_df" => results_df))

        # Atomic rename
        mv(temp_filepath, filepath; force = true)
        @info "Saved optimization results to $filepath"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_filepath) && rm(temp_filepath; force = true)
        rethrow(e)
    end

    return nothing
end
