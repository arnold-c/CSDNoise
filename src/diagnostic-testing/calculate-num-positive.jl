export calculate_positives_vec!

"""
    calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64
    )

Calculate test-positive individuals accounting for test result lag.

Applies the test performance multiplier (sensitivity or 1-specificity) to the number
tested, with results appearing after the specified lag period.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of individuals tested
- `sim_length`: Length of simulation
- `lag`: Test result lag in days
- `tested_multiplier`: Test performance multiplier (sensitivity or 1-specificity)
"""
function calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64
    )
    # Initialize to zero
    fill!(npos_vec, 0)

    @inbounds for day in eachindex(npos_vec)
        result_day = day + lag
        if result_day <= sim_length
            npos_vec[result_day] = rand(Distributions.Binomial(tested_vec[day], tested_multiplier))
        end
    end

    return nothing
end

function calculate_noise_positives!(outvec, tested_vec, tlength, lag, spec)
    tested_multiplier = 1.0 - spec
    calculate_positives!(outvec, tested_vec, tlength, lag, tested_multiplier)
    return nothing
end

function calculate_true_positives!(outvec, tested_vec, tlength, lag, sens)
    calculate_positives!(outvec, tested_vec, tlength, lag, sens)
    return nothing
end

function calculate_positives!(
        npos_vec, tested_vec, tlength, lag, tested_multiplier
    )
    @inbounds for day in eachindex(tested_vec)
        if day + lag <= tlength
            npos_vec[day + lag] = Int64(
                round(tested_vec[day] * tested_multiplier)
            )
        end
    end
    return nothing
end
