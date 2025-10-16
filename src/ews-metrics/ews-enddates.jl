export calculate_all_ews_enddates,
    calculate_ews_enddate

"""
    calculate_all_ews_enddates(
        vec_of_thresholds::Vector{T},
        ews_enddate_type::EWSEndDateType
    ) where {T<:AbstractThresholds}

Calculate EWS endpoints for all simulations based on enddate type.

# Arguments
- `thresholds`: Vector of threshold objects (one per simulation)
- `ews_enddate_type`: Type of endpoint to calculate (ReffStart, ReffEnd, OutbreakStart, etc.)

# Returns
- `Vector{Int64}`: Vector of endpoints, one per simulation

# Example
```julia
endpoints = calculate_ews_endpoints(thresholds, EWSEndDateType(ReffStart()))
```
"""
function calculate_all_ews_enddates(
        vec_of_thresholds::StructVector{T},
        ews_enddate_type::EWSEndDateType
    )::Union{Try.Ok, Try.Err} where {T <: AbstractThresholds}

    nsims = length(vec_of_thresholds)
    enddates = Vector{Int64}(undef, nsims)

    for (sim, thresholds) in pairs(vec_of_thresholds)
        ews_enddate_result = calculate_ews_enddate(thresholds, ews_enddate_type)
        Try.iserr(ews_enddate_result) && return ews_enddate_result
        enddates[sim] = Try.unwrap(ews_enddate_result)
    end

    return Try.Ok(enddates)
end


"""
    calculate_ews_enddate(thresholds, enddate_type)

Calculate a single EWS endpoint for a given threshold object and enddate type.

This function computes the endpoint for early warning signal analysis based on
the specified enddate type and threshold data. It serves as a wrapper around
the internal `_calculate_ews_enddate` function, providing error handling and
meaningful error messages when endpoint calculation fails.

# Arguments
- `thresholds::AbstractThresholds`: Threshold object containing boundary information
- `enddate_type::EWSEndDateType`: Type of endpoint to calculate

# Supported Enddate Types
- `ReffStart`: Start of effective reproduction number threshold period
- `ReffEnd`: End of effective reproduction number threshold period
- `OutbreakStart`: Start of outbreak threshold period
- `OutbreakEnd`: End of outbreak threshold period
- `OutbreakMiddle`: Middle point of outbreak threshold period

# Returns
- `Try.Ok{Int64}`: Successfully calculated endpoint index
- `Try.Err{String}`: Error with descriptive message if calculation fails

# Error Conditions
Returns an error if:
- The threshold object doesn't contain the required boundary data
- The enddate type is incompatible with the threshold type
- Vaccination rate post burn-in is too high
- Simulation length is insufficient

# Example
```julia
# Calculate outbreak start endpoint
result = calculate_ews_enddate(outbreak_thresholds, EWSEndDateType(OutbreakStart()))
if Try.isok(result)
    endpoint = Try.unwrap(result)
    println("Outbreak starts at index: $endpoint")
else
    println("Error: $(Try.unwrap_err(result))")
end
```

# See Also
- [`calculate_all_ews_enddates`](@ref): Batch calculation for multiple simulations
- [`_calculate_ews_enddate`](@ref): Internal implementation functions
- [`AbstractThresholds`](@ref): Base type for threshold objects
"""
function calculate_ews_enddate(
        thresholds::T,
        enddate_type::EWSEndDateType,
    ) where {T <: AbstractThresholds}

    enddate = _calculate_ews_enddate(thresholds, enddate_type)

    if Try.isok(enddate)
        return enddate
    end

    return Try.Err(
        "Failed to calculate ews_enddate for $enddate_type. The vaccination rate post burnin may need to be lowered, or increase the simulation length"
    )
end

_calculate_ews_enddate(
    thresholds::T,
    enddate_type::EWSEndDateType
) where {T <: AbstractThresholds} = _calculate_ews_enddate(thresholds, LightSumTypes.variant(enddate_type))

_calculate_ews_enddate(
    thresholds::Thresholds,
    enddate_type::ReffStart
) = TryExperimental.trygetindex(thresholds.lower_bounds, 1)

_calculate_ews_enddate(
    thresholds::Thresholds,
    enddate_type::ReffEnd
) = TryExperimental.trygetindex(thresholds.upper_bounds, 1)

_calculate_ews_enddate(
    thresholds::OutbreakThresholds,
    enddate_type::OutbreakStart
) = TryExperimental.trygetindex(thresholds.lower_bounds, 1)

_calculate_ews_enddate(
    thresholds::OutbreakThresholds,
    enddate_type::OutbreakEnd
) = TryExperimental.trygetindex(thresholds.upper_bounds, 1)

function _calculate_ews_enddate(thresholds::OutbreakThresholds, enddate_type::OutbreakMiddle)
    return Try.Ok(
        TryExperimental.trygetindex(thresholds.lower_bounds, 1) +
            div((TryExperimental.trygetindex(thresholds.duration, 1), 2))
    )
end
