export calculate_min_max_vaccination_range,
    calculate_vaccination_rate_to_achieve_Reff

function calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range = 0.2,
    )
    @assert mean_vaccination_coverage <= 1.0

    min_vaccination_range = minimum(
        [
            max_vaccination_range,
            1.0 - mean_vaccination_coverage,
            mean_vaccination_coverage,
        ]
    )

    min_vaccination_coverage = round(
        mean_vaccination_coverage - min_vaccination_range; digits = 4
    )
    max_vaccination_coverage = round(
        mean_vaccination_coverage + min_vaccination_range; digits = 4
    )
    return min_vaccination_coverage, max_vaccination_coverage
end

export calculate_vaccination_rate_to_achieve_Reff

"""
    calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, initial_states, dynamics_specification
    )

Assuming that in the burnin period the rate of infections is negligible, the time to reach a target Reff can be calculated by the difference in the rates in and out of the Susceptible group.

μN(1-ρ) -> S -> μS
dS/dt = μ(N - Nρ - S)
"""
function calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, initial_states, R_0, mu
    )
    @unpack S, N = initial_states

    return calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, S, N, R_0, mu
    )
end

function calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, S, N, R_0, mu
    )
    @assert target_Reff > 0
    @assert target_Reff < 1.2

    vaccination_coverage =
        1 -
        (((target_Reff * N) / R_0 - S) / (365 * target_years * mu) + S) / N

    if vaccination_coverage > 1 || vaccination_coverage < 0
        return error(
            "Target Reff cannot be reached in the burn-in period. Initial Reff = $(R_0 * S / N). Try a longer burn-in period or a smaller target Reff."
        )
    end
    return round(vaccination_coverage; digits = 4)
end

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
