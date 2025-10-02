export SeasonalityFunction,
    CosineSeasonality,
    SineSeasonality,
    TargetDiseaseDynamicsParameters,
    CommonDiseaseDynamicsParameters,
    DynamicsParameterSpecification,
    DynamicsParameters

abstract type AbstractSeasonalityFunction end
struct CosineSeasonality end
struct SineSeasonality end
LightSumTypes.@sumtype SeasonalityFunction(
    CosineSeasonality,
    SineSeasonality
) <: AbstractSeasonalityFunction


Base.@kwdef struct TargetDiseaseDynamicsParameters
    R_0::Float64
    latent_period_days::Float64
    infectious_duration_days::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    max_burnin_vaccination_coverage::Float64
    min_post_burnin_vaccination_coverage::Float64
    max_post_burnin_vaccination_coverage::Float64
end

Base.@kwdef struct CommonDiseaseDynamicsParameters
    births_per_k_pop::Float64
    nsims::Int64
    burnin_target_Reff::Float64
end

"""
    DynamicsParameterSpecification

Setup struct containing the majority of parameters required for constructing a `DynamicsParameters` instance.

The constructor for `DynamicsParameters` samples from the burnin and post-burnin vaccination coverage ranges to generate specific vaccination coverage values for each simulation.

# Fields
- `beta_mean::Float64`: Mean transmission rate
- `beta_force::Float64`: Forcing amplitude for seasonal transmission
- `seasonality::SeasonalityFunction`: Type of seasonal forcing function
- `sigma::Float64`: Rate of progression from latent to infectious (1/latent_period)
- `gamma::Float64`: Recovery rate (1/infectious_duration)
- `mu::Float64`: Birth/death rate
- `annual_births_per_k::Float64`: Annual births per 1000 population
- `epsilon::Float64`: Vaccine efficacy
- `R_0::Float64`: Basic reproduction number
- `min_burnin_vaccination_coverage::Float64`: Minimum vaccination coverage during burnin
- `max_burnin_vaccination_coverage::Float64`: Maximum vaccination coverage during burnin
- `min_post_burnin_vaccination_coverage::Float64`: Minimum vaccination coverage post-burnin
- `max_post_burnin_vaccination_coverage::Float64`: Maximum vaccination coverage post-burnin
"""
Base.@kwdef struct DynamicsParameterSpecification
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    min_burnin_vaccination_coverage::Float64
    max_burnin_vaccination_coverage::Float64
    min_post_burnin_vaccination_coverage::Float64
    max_post_burnin_vaccination_coverage::Float64
end

function DynamicsParameterSpecification(
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_post_burnin_vaccination_coverage::Nothing,
        max_post_burnin_vaccination_coverage::Nothing,
    )
    return DynamicsParameterSpecification(
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
    )
end

"""
    DynamicsParameters

Complete set of disease dynamics parameters for a single simulation, including sampled vaccination coverage values.

This struct contains all parameters from `DynamicsParameterSpecification` plus the specific `burnin_vaccination_coverage` and `post_burnin_vaccination_coverage` values sampled from their respective ranges for this particular simulation instance.

# Fields
- `beta_mean::Float64`: Mean transmission rate
- `beta_force::Float64`: Forcing amplitude for seasonal transmission
- `seasonality::SeasonalityFunction`: Type of seasonal forcing function
- `sigma::Float64`: Rate of progression from latent to infectious (1/latent_period)
- `gamma::Float64`: Recovery rate (1/infectious_duration)
- `mu::Float64`: Birth/death rate
- `annual_births_per_k::Float64`: Annual births per 1000 population
- `epsilon::Float64`: Vaccine efficacy
- `R_0::Float64`: Basic reproduction number
- `min_burnin_vaccination_coverage::Float64`: Minimum vaccination coverage during burnin
- `max_burnin_vaccination_coverage::Float64`: Maximum vaccination coverage during burnin
- `min_post_burnin_vaccination_coverage::Float64`: Minimum vaccination coverage post-burnin
- `max_post_burnin_vaccination_coverage::Float64`: Maximum vaccination coverage post-burnin
- `burnin_vaccination_coverage::Float64`: Sampled vaccination coverage for burnin period
- `post_burnin_vaccination_coverage::Float64`: Sampled vaccination coverage for post-burnin period
"""
Base.@kwdef struct DynamicsParameters
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    min_burnin_vaccination_coverage::Float64
    max_burnin_vaccination_coverage::Float64
    min_post_burnin_vaccination_coverage::Float64
    max_post_burnin_vaccination_coverage::Float64
    burnin_vaccination_coverage::Float64
    post_burnin_vaccination_coverage::Float64
end

function DynamicsParameters(
        dynamic_parameter_specification::DynamicsParameterSpecification;
        seed = 1234
    )
    Random.seed!(seed)

    burnin_vaccination_coverage =
    if dynamic_parameter_specification.min_burnin_vaccination_coverage ==
            dynamic_parameter_specification.max_burnin_vaccination_coverage
        dynamic_parameter_specification.min_burnin_vaccination_coverage
    else
        sample_vaccination_coverage(
            dynamic_parameter_specification.min_burnin_vaccination_coverage,
            dynamic_parameter_specification.max_burnin_vaccination_coverage,
        )
    end

    post_burnin_vaccination_coverage =
    if dynamic_parameter_specification.min_burnin_vaccination_coverage ==
            dynamic_parameter_specification.min_post_burnin_vaccination_coverage &&
            dynamic_parameter_specification.max_burnin_vaccination_coverage ==
            dynamic_parameter_specification.max_post_burnin_vaccination_coverage
        burnin_vaccination_coverage
    else
        sample_vaccination_coverage(
            dynamic_parameter_specification.min_post_burnin_vaccination_coverage,
            dynamic_parameter_specification.max_post_burnin_vaccination_coverage,
        )
    end

    dynamics_parameters = DynamicsParameters(
        dynamic_parameter_specification.beta_mean,
        dynamic_parameter_specification.beta_force,
        dynamic_parameter_specification.seasonality,
        dynamic_parameter_specification.sigma,
        dynamic_parameter_specification.gamma,
        dynamic_parameter_specification.mu,
        dynamic_parameter_specification.annual_births_per_k,
        dynamic_parameter_specification.epsilon,
        dynamic_parameter_specification.R_0,
        dynamic_parameter_specification.min_burnin_vaccination_coverage,
        dynamic_parameter_specification.max_burnin_vaccination_coverage,
        dynamic_parameter_specification.min_post_burnin_vaccination_coverage,
        dynamic_parameter_specification.max_post_burnin_vaccination_coverage,
        burnin_vaccination_coverage,
        post_burnin_vaccination_coverage,
    )
    return dynamics_parameters
end
