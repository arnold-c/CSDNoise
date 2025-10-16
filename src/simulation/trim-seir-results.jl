export trim_ensemble_simulations,
    trim_seir_results

"""
    trim_ensemble_simulations(ensemble_run, enddates)

Trim ensemble simulation results to specified end dates.

This function takes an existing ensemble of SEIR simulations and trims both
the emergent and null dynamics scenarios to the specified end dates. This is
useful for post-processing ensemble results to focus analysis on specific
time periods without re-running the full simulation, and required before
creating the noise simulations.

# Arguments
- `ensemble_run::EnsembleSEIRRun`: An ensemble simulation result containing
  both emergent and null SEIR runs as StructVectors of individual simulation
  results.
- `enddates::Vector{Int64}`: Vector of end dates (in simulation time units)
  specifying where to trim each simulation in the ensemble. The length should
  match the number of simulations in the ensemble.

# Returns
- `EnsembleSEIRRun`: A new ensemble structure with both emergent and null
  SEIR runs trimmed to the specified end dates. The structure maintains the
  same format as the input but with shortened time series.
"""
function trim_ensemble_simulations(
        ensemble_run::EnsembleSEIRRun,
        enddates::Vector{Int64}
    )

    trimmed_emergent_run = trim_seir_results(
        ensemble_run.emergent_seir_run,
        enddates
    )

    trimmed_null_run = trim_seir_results(
        ensemble_run.null_seir_run,
        enddates
    )

    return EnsembleSEIRRun(
        emergent_seir_run = trimmed_emergent_run,
        null_seir_run = trimmed_null_run
    )

end


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
