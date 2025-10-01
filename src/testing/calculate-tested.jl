export calculate_tested_vec!

"""
    calculate_tested_vec!(
        outvec::AbstractVector{Int64},
        invec::Vector{Int64},
        perc_tested::Float64
    )

Calculate the number of individuals tested from incidence data.

Applies the testing percentage to each day's incidence, rounding to nearest integer.

# Arguments
- `outvec`: Output vector to store number tested (modified in-place)
- `invec`: Input incidence vector
- `perc_tested`: Percentage of individuals tested (0.0 to 1.0)
"""
function calculate_tested_vec!(
        outvec::AbstractVector{Int64},
        invec::Vector{Int64},
        perc_tested::Float64
    )
    @inbounds for i in eachindex(invec)
        outvec[i] = round(Int64, invec[i] * perc_tested)
    end

    return nothing
end

function calculate_tested!(outvec, invec, perc_tested)
    return @. outvec = round(invec * perc_tested)
end
