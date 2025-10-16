export get_ews_metric_specification_description,
    method_string

"""
    get_ews_metric_specification_description(ews_metric_specification)

Generate a human-readable description string for an EWS metric specification.

This function creates a formatted string that summarizes the key parameters of an
`EWSMetricSpecification`, providing a concise overview of the calculation settings
for display, logging, or debugging purposes.

# Arguments
- `ews_metric_specification::EWSMetricSpecification`: The EWS specification to describe

# Returns
- `String`: Formatted description containing method, aggregation, bandwidth, and lag information

# Example
```julia
spec = EWSMetricSpecification(
    EWSMethod(Backward()),
    Dates.Day(7),
    Dates.Day(30),
    0
)
description = get_ews_metric_specification_description(spec)
# Returns: "Method: backward, Aggregation: 7 days, Bandwidth: 30 days, Lag: 0"
```

# See Also
- [`method_string`](@ref): Converts EWSMethod to lowercase string representation
- [`EWSMetricSpecification`](@ref): The specification type being described
"""
function get_ews_metric_specification_description(ews_metric_specification)
    return "Method: $(method_string(ews_metric_specification.method)), Aggregation: $(ews_metric_specification.aggregation), Bandwidth: $(ews_metric_specification.bandwidth), Lag: $(ews_metric_specification.lag)"
end

"""
    method_string(method::EWSMethod)

Convert an EWSMethod to its lowercase string representation.

This function extracts the method name from an `EWSMethod` sum type and returns
it as a lowercase string, suitable for use in file paths, descriptions, or
other string-based identifiers.

# Arguments
- `method::EWSMethod`: The EWS method to convert (either Backward or Centered)

# Returns
- `String`: Lowercase string representation of the method ("backward" or "centered")

# Example
```julia
backward_method = EWSMethod(Backward())
method_str = method_string(backward_method)
# Returns: "backward"

centered_method = EWSMethod(Centered())
method_str = method_string(centered_method)
# Returns: "centered"
```

# Implementation Notes
This function parses the string representation of the sum type to extract the
variant name, removing the parentheses and converting to lowercase.

# See Also
- [`EWSMethod`](@ref): The sum type representing EWS calculation methods
- [`get_ews_metric_specification_description`](@ref): Uses this function for method descriptions
"""
function method_string(method::EWSMethod)
    return lowercase(split(string(method), "(")[2])
end
