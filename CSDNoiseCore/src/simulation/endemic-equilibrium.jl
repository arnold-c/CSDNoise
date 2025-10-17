export calculate_endemic_equilibrium_proportions

"""
    calculate_endemic_equilibrium_proportions(dynamics_params, vaccination_coverage)

Calculate endemic equilibrium proportions for an SEIR model with vaccination.

# Arguments
- `dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification}`: Disease dynamics parameters
- `vaccination_coverage::Float64`: Proportion of population vaccinated (must be in [0, 1))

# Returns
A named tuple containing:
- `s_prop`: Proportion susceptible at equilibrium
- `e_prop`: Proportion exposed at equilibrium
- `i_prop`: Proportion infectious at equilibrium
- `r_prop`: Proportion recovered at equilibrium
"""
function calculate_endemic_equilibrium_proportions(
        dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification},
        vaccination_coverage::Float64
    )
    return calculate_endemic_equilibrium_proportions(
        dynamics_params.R_0,
        dynamics_params.beta_mean,
        dynamics_params.sigma,
        dynamics_params.gamma,
        dynamics_params.mu,
        vaccination_coverage
    )
end

"""
    calculate_endemic_equilibrium_proportions(R_0, beta_mean, sigma, gamma, mu, vaccination_coverage)

Calculate endemic equilibrium proportions for an SEIR model with vaccination using explicit parameters.

# Arguments
- `R_0::Float64`: Basic reproduction number
- `beta_mean::Float64`: Mean transmission rate
- `sigma::Float64`: Rate of progression from exposed to infectious (1/latent period)
- `gamma::Float64`: Recovery rate (1/infectious period)
- `mu::Float64`: Birth/death rate
- `vaccination_coverage::Float64`: Proportion of population vaccinated (must be in [0, 1])

# Returns
A named tuple containing:
- `s_prop`: Proportion susceptible at equilibrium
- `e_prop`: Proportion exposed at equilibrium
- `i_prop`: Proportion infectious at equilibrium
- `r_prop`: Proportion recovered at equilibrium

# Errors
- Throws error if `vaccination_coverage` is not in [0, 1]
- Throws error if effective R₀ ≤ 1 (no endemic equilibrium exists)
- Throws error if calculated recovered proportion is negative (parameter inconsistency)
"""
function calculate_endemic_equilibrium_proportions(
        R_0::Float64,
        beta_mean::Float64,
        sigma::Float64,
        gamma::Float64,
        mu::Float64,
        vaccination_coverage::Float64
    )
    if vaccination_coverage < 0.0 || vaccination_coverage > 1.0
        error("Vaccination coverage rho must be in [0, 1]. Got rho = $vaccination_coverage")
    end

    R_eff = R_0 * (1.0 - vaccination_coverage)

    if R_eff <= 1.0
        return Try.Err(
            "Effective R₀ must be > 1 for endemic equilibrium. Got R_eff = R₀(1 - ρ) = $(round(R_eff; digits = 2)). "
        )
    end

    s_prop = 1.0 / R_0

    i_prop = mu * (R_eff - 1.0) / beta_mean

    e_prop = i_prop * (gamma + mu) / sigma

    r_prop = 1.0 - s_prop - e_prop - i_prop

    if r_prop < 0.0
        error("Calculated negative recovered proportion. Check parameter consistency.")
    end

    return Try.Ok(
        (
            s_prop = s_prop,
            e_prop = e_prop,
            i_prop = i_prop,
            r_prop = r_prop,
        )
    )
end
