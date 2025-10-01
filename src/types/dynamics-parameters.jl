export SeasonalityFunction,
    CosineSeasonality,
    SineSeasonality,
    DynamicsParameterSpecification,
    DynamicsParameters

abstract type AbstractSeasonalityFunction end
struct CosineSeasonality end
struct SineSeasonality end
LightSumTypes.@sumtype SeasonalityFunction(
    CosineSeasonality,
    SineSeasonality
) <: AbstractSeasonalityFunction


struct DynamicsParameterSpecification
    contact_matrix::Matrix{Int64}
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
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
end

function DynamicsParameterSpecification(
        contact_matrix,
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
        min_vaccination_coverage::Nothing,
        max_vaccination_coverage::Nothing,
    )
    return DynamicsParameterSpecification(
        contact_matrix,
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

struct DynamicsParameters
    contact_matrix::Matrix{Int64}
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
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
    burnin_vaccination_coverage::Float64
    vaccination_coverage::Float64
end

function DynamicsParameters(
        dynamic_parameter_specification::DynamicsParameterSpecification; seed = 1234
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

    vaccination_coverage =
    if dynamic_parameter_specification.min_burnin_vaccination_coverage ==
            dynamic_parameter_specification.min_vaccination_coverage &&
            dynamic_parameter_specification.max_burnin_vaccination_coverage ==
            dynamic_parameter_specification.max_vaccination_coverage
        burnin_vaccination_coverage
    else
        sample_vaccination_coverage(
            dynamic_parameter_specification.min_vaccination_coverage,
            dynamic_parameter_specification.max_vaccination_coverage,
        )
    end

    dynamics_parameters = DynamicsParameters(
        dynamic_parameter_specification.contact_matrix,
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
        dynamic_parameter_specification.min_vaccination_coverage,
        dynamic_parameter_specification.max_vaccination_coverage,
        burnin_vaccination_coverage,
        vaccination_coverage,
    )
    return dynamics_parameters
end
