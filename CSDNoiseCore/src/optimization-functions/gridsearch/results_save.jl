export save_results_structvector

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
        JLD2.jldsave(temp_filepath; optimization_results = results)

        # Atomic rename
        mv(temp_filepath, filepath; force = true)
        @info "Saved optimization results to $filepath"

    catch e
        # Clean up temp file if something went wrong
        isfile(temp_filepath) && rm(temp_filepath; force = true)
        rethrow(e)
    end
end
