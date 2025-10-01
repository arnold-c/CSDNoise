# Modern logging utilities for benchmarking scripts
# Provides functions that output to both console (with styling) and markdown files

using Dates
using Logging
using LoggingExtras
using StyledStrings

export styled_to_markdown, format_markdown_table, setup_dual_logging,
    cleanup_logging, log_both

# Global variables for logging
const ORIGINAL_LOGGER = Ref{AbstractLogger}()
const MARKDOWN_FILE = Ref{IOStream}()

"""
Helper function for dual output - console with styling, markdown clean
"""
function log_both(msg::AbstractString)
    # Print styled version to console
    println(msg)

    # Convert to markdown and log (will only go to file now)
    markdown_msg = styled_to_markdown(msg)
    return @info markdown_msg
end

"""
Convert styled string to markdown format by analyzing annotations and applying transformation rules
"""
function styled_to_markdown(styled_str::AbstractString)
    text = string(styled_str)

    # Simple pattern matching rules
    rules = [
        # Add spacing before lines ending with ":"
        (r":$", m -> "\n" * text),

        # Convert "  label: value" to markdown list items with bold labels
        (r"^  ([^:]+): (.+)$", m -> "- **$(m.captures[1]):** $(m.captures[2])"),
    ]

    # Apply transformation rules
    for (pattern, transform) in rules
        m = match(pattern, text)
        if m !== nothing
            return transform(m)
        end
    end

    # Default: return text as-is
    return text
end

"""
Format a table as markdown
"""
function format_markdown_table(headers::Vector{String}, rows::Vector{Vector{String}})
    if isempty(headers) || isempty(rows)
        return ""
    end

    # Create header row
    header_row = "| " * join(headers, " | ") * " |"

    # Create separator row
    separator = "| " * join(["---" for _ in headers], " | ") * " |"

    # Create data rows
    data_rows = [
        "| " * join(row, " | ") * " |"
            for row in rows
    ]

    return join([header_row, separator, data_rows...], "\n")
end

"""
Setup dual logging to console and markdown file
"""
function setup_dual_logging(script_name::String = "benchmark"; title::String = "Benchmark Report", output_dir::String = ".", filename_base::String = "")
    # Store original logger
    ORIGINAL_LOGGER[] = global_logger()

    # Create markdown filename
    if isempty(filename_base)
        # Use default timestamped naming
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        filename = "$(script_name)_$(timestamp).md"
    else
        # Use provided base name with _log.md suffix
        filename = "$(filename_base)_log.md"
    end

    # Create full path with output directory
    full_path = joinpath(output_dir, filename)

    # Open markdown file
    MARKDOWN_FILE[] = open(full_path, "w")

    # Create FormatLogger for clean markdown output
    markdown_logger = FormatLogger(MARKDOWN_FILE[]) do io, args
        # Only output the message content, no metadata
        println(io, args.message)
    end

    # Set ONLY the markdown logger as global (not a TeeLogger)
    global_logger(markdown_logger)

    # Write markdown header
    @info "# $title"
    @info ""
    @info "**Generated:** $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))"
    @info "**Script:** $script_name"
    @info ""

    return full_path
end

"""
Cleanup logging and close files
"""
function cleanup_logging()
    return try
        # Restore original logger
        global_logger(ORIGINAL_LOGGER[])

        # Close markdown file
        if isassigned(MARKDOWN_FILE) && isopen(MARKDOWN_FILE[])
            close(MARKDOWN_FILE[])
        end
    catch e
        @warn "Error during logging cleanup: $e"
    end
end
