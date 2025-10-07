export calculate_beta,
    calculate_gamma,
    calculateReffective,
    calculateR0,
    calculate_mu,
    calculate_import_rate

"""
    calculate_beta(R_0, sigma, gamma, mu)

Calculate the transmission rate beta for an SEIR model given R_0 and other parameters.

For an SEIR model with a single population:
R_0 = (β * σ) / ((σ + μ) * (γ + μ))

Solving for β:
β = R_0 * (σ + μ) * (γ + μ) / σ

# Arguments
- `R_0`: Basic reproduction number
- `sigma`: Rate of progression from E to I (1/latent_period)
- `gamma`: Recovery rate (1/infectious_period)
- `mu`: Death rate
"""
function calculate_beta(
        R_0::Float64,
        sigma::Float64,
        gamma::Float64,
        mu::Float64
    )::Float64
    return R_0 * (sigma + mu) * (gamma + mu) / sigma
end

"""
    calculate_gamma(R_0, beta, sigma, mu)

Calculate the recovery rate gamma for an SEIR model given R_0, beta, and other parameters.

For an SEIR model with a single population:
R_0 = (β * σ) / ((σ + μ) * (γ + μ))

Solving for γ:
γ = (β * σ) / (R_0 * (σ + μ)) - μ

# Arguments
- `R_0`: Basic reproduction number
- `beta`: Transmission rate
- `sigma`: Rate of progression from E to I (1/latent_period)
- `mu`: Death rate
"""
function calculate_gamma(
        R_0::Float64,
        beta::Float64,
        sigma::Float64,
        mu::Float64
    )::Float64
    gamma = (beta * sigma) / (R_0 * (sigma + mu)) - mu

    if gamma <= 0
        error("Calculated gamma is non-positive. Check parameter consistency.")
    end

    return gamma
end

"""
    calculate_beta_amp!(beta_vec, dynamics_parameters, time_parameters)

Calculate the amplitude of the transmission rate beta as a function of time for each time point,
storing results in-place in `beta_vec`.

This function modifies `beta_vec` in-place, calculating the time-varying transmission rate
using the mean transmission rate, seasonal forcing amplitude, and seasonality function
specified in `dynamics_parameters`.

# Arguments
- `beta_vec`: Vector to store the calculated beta amplitudes (modified in-place)
- `dynamics_parameters`: Parameters containing beta_mean, beta_force, and seasonality
- `time_parameters`: Time parameters containing the time range (trange)

# Returns
- `nothing` (modifies `beta_vec` in-place)
"""
function calculate_beta_amp!(
        beta_vec::V,
        dynamics_parameters::Union{DynamicsParameterSpecification, DynamicsParameters},
        time_parameters::SimTimeParameters
    ) where {V <: AbstractVector{<:AbstractFloat}}
    beta_mean = dynamics_parameters.beta_mean
    beta_force = dynamics_parameters.beta_force
    seasonality = dynamics_parameters.seasonality
    trange = time_parameters.trange

    for i in eachindex(beta_vec)
        beta_vec[i] = _calculate_beta_amp(
            beta_mean,
            beta_force,
            trange[i],
            LightSumTypes.variant(seasonality)
        )
    end
    return nothing
end

"""
    calculate_beta_amp(beta_mean, beta_force, t; seasonality)

Calculate the amplitude of the transmission rate beta as a function of time.
`beta_mean` is the mean transmission rate, `beta_force` is the amplitude of the `seasonality` function.
`seasonality` should be a SeasonalityFunction sum type.
"""
function calculate_beta_amp(beta_mean, beta_force, t; seasonality::SeasonalityFunction)
    return _calculate_beta_amp(beta_mean, beta_force, t, LightSumTypes.variant(seasonality))
end

# Internal dispatch functions for seasonality variants

"""
    _calculate_beta_amp(beta_mean, beta_force, t, ::CosineSeasonality)

Internal function to calculate beta amplitude using cosine seasonality.
Returns `beta_mean * (1 + beta_force * cos(2π * t / 365))`.
"""
_calculate_beta_amp(beta_mean, beta_force, t, ::CosineSeasonality) =
    beta_mean * (1 + beta_force * cos(2π * t / 365))

"""
    _calculate_beta_amp(beta_mean, beta_force, t, ::SineSeasonality)

Internal function to calculate beta amplitude using sine seasonality.
Returns `beta_mean * (1 + beta_force * sin(2π * t / 365))`.
"""
_calculate_beta_amp(beta_mean, beta_force, t, ::SineSeasonality) =
    beta_mean * (1 + beta_force * sin(2π * t / 365))

"""
    calculateReffective_t!(Reff_vec, beta_vec, dynamics_params, seir_arr)

Calculate the effective reproduction number, R_eff, at each time step for a given set of parameters.
"""
function calculateReffective_t!(
        Reff_vec::AbstractVector{Float64},
        beta_vec::Vector{Float64},
        dynamics_params::DynamicsParameters,
        seir_arr
    )::Nothing
    for i in eachindex(Reff_vec)
        Reff_vec[i] = calculateReffective(
            beta_vec[i],
            dynamics_params,
            seir_arr[i][1],
        )
    end

    return nothing
end

"""
    calculateReffective(beta_t, dynamics_params, S, N)

Calculate the effective reproduction number, R_eff, for a given set of parameters.
R_eff = R_0 * (S / N), where R_0 is calculated from beta and other parameters.
"""
function calculateReffective(
        beta_t::Float64,
        dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification},
        S::Int64,
        N::Int64
    )::Float64
    R_0 = calculateR0(beta_t, dynamics_params)
    Reff = R_0 * (S / N)

    return Reff
end

"""
    calculateReffective(beta_t, dynamics_params, vaccination_coverage)

Calculate the effective reproduction number, R_eff, for a given set of parameters.
R_eff = R_0 * (1 - vaccination_coverage), where R_0 is calculated from beta and other parameters.
"""
function calculateReffective(
        beta_t::Float64,
        dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification},
        vaccination_coverage::Float64,
    )::Float64
    R_0 = calculateR0(beta_t, dynamics_params)
    Reff = R_0 * (1 - vaccination_coverage)

    return Reff
end


"""
    calculateR0(beta, dynamics_params)

Calculate the basic reproduction number R_0 for an SEIR model using a DynamicsParameters struct.

This is a convenience wrapper that extracts sigma, gamma, and mu from the dynamics_params
and calls the main calculateR0 function.

# Arguments
- `beta`: Transmission rate
- `dynamics_params`: DynamicsParameters struct containing sigma, gamma, and mu

# Returns
- `Float64`: The basic reproduction number R_0
"""
function calculateR0(
        beta::Float64,
        dynamics_params::Union{DynamicsParameters, DynamicsParameterSpecification}
    )
    sigma = dynamics_params.sigma
    gamma = dynamics_params.gamma
    mu = dynamics_params.mu
    R_0 = calculateR0(beta, sigma, gamma, mu)
    return R_0
end

"""
    calculateR0(beta, sigma, gamma, mu)

Calculate the basic reproduction number R_0 for an SEIR model.

For an SEIR model with a single population:
R_0 = (β * σ) / ((σ + μ) * (γ + μ))

# Arguments
- `beta`: Transmission rate
- `sigma`: Rate of progression from E to I (1/latent_period)
- `gamma`: Recovery rate (1/infectious_period)
- `mu`: Death rate

# Returns
- `Float64`: The basic reproduction number R_0
"""
@inline function calculateR0(
        beta::Float64,
        sigma::Float64,
        gamma::Float64,
        mu::Float64
    )::Float64
    return (beta * sigma) / ((sigma + mu) * (gamma + mu))
end

"""
    calculate_mu(annual_births_per_k)

Calculate the death rate (mu) from annual births per thousand population.

This function converts annual births per thousand to a daily death rate by:
1. Converting births per thousand to life expectancy in years (1000 / annual_births_per_k)
2. Converting life expectancy to daily death rate (1 / (life_expectancy_years * 365))

# Arguments
- `annual_births_per_k`: Annual births per thousand population

# Returns
- `Float64`: Daily death rate (mu)
"""
function calculate_mu(annual_births_per_k)
    life_expectancy_years = 1000 / annual_births_per_k
    return 1 / (life_expectancy_years * 365)
end

"""
    calculate_import_rate(mu, R_0, N)

Calulate the rate of new infectious individuals imported into the simulation using the commuter import formula defined in p210 of Keeling & Rohani
"""
function calculate_import_rate(mu, R_0, N)
    return (1.06 * mu * (R_0 - 1)) / sqrt(N)
end
