export sample_vaccination_coverage

"""
    sample_vaccination_coverage(min_coverage, max_coverage, digits = 4)

Sample a vaccination coverage value from a uniform distribution between specified bounds.

This function generates a random vaccination coverage value uniformly distributed
between the minimum and maximum coverage values, rounded to the specified number
of decimal places.

# Arguments
- `min_coverage`: Minimum vaccination coverage (lower bound)
- `max_coverage`: Maximum vaccination coverage (upper bound)
- `digits`: Number of decimal places for rounding (default: 4)

# Returns
- `Float64`: Randomly sampled vaccination coverage value
"""
function sample_vaccination_coverage(
        min_coverage,
        max_coverage,
        digits = 4
    )
    return round(
        rand(
            Distributions.Uniform(
                min_coverage,
                max_coverage
            )
        );
        digits = digits
    )
end
