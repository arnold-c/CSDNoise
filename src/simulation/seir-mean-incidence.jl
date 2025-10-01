using StructArrays: StructVector
using StatsBase: StatsBase
using Bumper: @no_escape, @alloc

export calculate_mean_incidence

"""
    calculate_mean_incidence(seir_results::StructVector{SEIRRun})

Calculate mean incidence for each simulation up to its enddate.

# Arguments
- `seir_results`: StructVector of SEIR simulation results (either the full-length simulations, or pre-filtered to the enddate)

# Returns
- `overall_mean`: The overall mean

# Example
```julia
overall_mean = calculate_mean_incidence(seir_results)
```
"""
function calculate_mean_incidence(seir_results::StructVector{SEIRRun})

    nsims = length(seir_results)

    @no_escape begin
        incidence_means = @alloc(Float64, nsims)

        for sim in eachindex(seir_results.incidence)
            # Calculate mean up to the endpoint
            incidence_means[sim] = StatsBase.mean(seir_results[sim].incidence)
        end

        # Calculate overall mean across all simulations
        overall_mean = StatsBase.mean(incidence_means)
    end

    return overall_mean
end

# end
