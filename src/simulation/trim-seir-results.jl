export trim_seir_results

"""
    trim_seir_results(seir_results, enddates) -> StructVector{SEIRRun}

Create a filtered version of SEIRRun results, keeping only data up to specified endpoints
for incidence and states properties. The Reff property is preserved as-is.

# Arguments
- `seir_results`: StructVector of SEIRRun containing simulation results
- `enddates`: Vector of endpoints, one per simulation

# Returns
- `StructVector{SEIRRun}`: Filtered results with truncated incidence and states

# Example
```julia
filtered_results = trim_seir_results(seir_results, endpoints)
```
"""
function trim_seir_results(
        seir_results::StructVector{SEIRRun},
        enddates::Vector{Int64}
    )

    nsims = length(seir_results)
    @assert nsims == length(enddates) "Number of simulations must match number of endpoints"

    # Pre-allocate vectors for filtered data
    trimmed_incidence = Vector{Vector{Int64}}(undef, nsims)
    trimmed_states = Vector{Vector{StaticArrays.SVector{5, Int64}}}(undef, nsims)
    trimmed_Reff = Vector{Vector{Float64}}(undef, nsims)

    for (sim, enddate) in pairs(enddates)
        if enddate > 0 && enddate <= length(seir_results[sim].incidence)
            # Filter incidence up to endpoint
            trimmed_incidence[sim] = Vector{Int64}(
                seir_results[sim].incidence[1:enddate]
            )

            # Filter states up to endpoint
            trimmed_states[sim] = Vector{StaticArrays.SVector{5, Int64}}(
                seir_results[sim].states[1:enddate]
            )

            trimmed_Reff[sim] = Vector{Float64}(
                seir_results[sim].Reff[1:enddate]
            )
        else
            error("Sim $sim has enddate ($enddate) outside of the incidence vectors size ($(length(seir_results[sim].incidence)))")
        end
    end

    # Create new StructVector with filtered data, preserving Reff as-is
    return StructVector{SEIRRun}(
        states = trimmed_states,
        incidence = trimmed_incidence,
        Reff = trimmed_Reff
    )
end
