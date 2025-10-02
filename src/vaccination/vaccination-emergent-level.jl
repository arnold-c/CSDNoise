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
        time_parameters::SimTimeParameters,
        state_parameters::StateParameters,
        target_Reff,
        R_0,
        mu
    )
    @unpack S, N = state_parameters.init_states
    @unpack burnin = time_parameters

    return calculate_vaccination_rate_to_achieve_Reff(
        burnin,
        S,
        N,
        target_Reff,
        R_0,
        mu
    )
end

function calculate_vaccination_rate_to_achieve_Reff(
        target_days,
        S,
        N,
        target_Reff,
        R_0,
        mu
    )
    @assert target_Reff > 0
    @assert target_Reff < 1.2

    vaccination_coverage =
        1 -
        (((target_Reff * N) / R_0 - S) / (target_days * mu) + S) / N

    if vaccination_coverage > 1 || vaccination_coverage < 0
        return error(
            "Target Reff cannot be reached in the burn-in period. Initial Reff = $(R_0 * S / N). Try a longer burn-in period or a smaller target Reff."
        )
    end
    return round(vaccination_coverage; digits = 4)
end
