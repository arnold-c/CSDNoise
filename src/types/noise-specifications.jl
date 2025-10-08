export NoiseSpecification,
    PoissonNoise,
    DynamicalNoise,
    DynamicalNoiseSpecification,
    NoiseVaccinationOptimizationParameters

struct DynamicalNoiseSpecification
    R_0::Float64
    latent_period::Int64
    duration_infection::Int64
    correlation::String
    poisson_component::Float64
    vaccination_bounds::Vector{Float64}
    function DynamicalNoiseSpecification(
            R_0::Float64,
            latent_period::Int64,
            duration_infection::Int64,
            correlation::String,
            poisson_component::Float64,
            vaccination_bounds::Vector{Float64},
        )

        @assert length(vaccination_bounds) == 2
        @assert vaccination_bounds[1] < vaccination_bounds[2]
        return new(
            R_0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            vaccination_bounds,
        )
    end
end

function DynamicalNoiseSpecification(;
        R_0::Float64,
        latent_period::Int64,
        duration_infection::Int64,
        correlation::String,
        poisson_component::Float64,
        vaccination_bounds::Vector{Float64} = [0.0, 1.0],
    )
    return DynamicalNoiseSpecification(
        R_0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        vaccination_bounds,
    )
end

abstract type AbstractNoiseSpecification end

Base.@kwdef struct PoissonNoise
    noise_mean_scaling::Float64
end

Base.@kwdef struct DynamicalNoise
    R_0::Float64
    latent_period::Int64
    duration_infection::Int64
    correlation::String
    noise_mean_scaling::Float64
    vaccination_coverage::Float64
end

function DynamicalNoise(
        dynamical_noise_specification::DynamicalNoiseSpecification,
        vaccination_coverage
    )
    @assert vaccination_coverage >= 0.0 && vaccination_coverage <= 1.0
    return DynamicalNoise(
        R_0 = dynamical_noise_specification.R_0,
        latent_period = dynamical_noise_specification.latent_period,
        duration_infection = dynamical_noise_specification.duration_infection,
        correlation = dynamical_noise_specification.correlation,
        noise_mean_scaling = dynamical_noise_specification.poisson_component,
        vaccination_coverage = vaccination_coverage
    )
end

LightSumTypes.@sumtype NoiseSpecification(PoissonNoise, DynamicalNoise) <: AbstractNoiseSpecification


struct NoiseVaccinationOptimizationParameters
    n_sobol_points::Int64
    local_algorithm
    maxeval::Int64
    xtol_rel::Float64
    xtol_abs::Float64
    atol::Float64
end

"""
    NoiseVaccinationOptimizationParameters(;
        n_sobol_points::Int64 = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval::Int64 = 1000,
        xtol_rel::Float64 = 1.0e-3,
        xtol_abs::Float64 = 1.0e-3,
        atol::Float64 = 1.0e-4
    )

Parameters for optimizing noise vaccination levels.

# Arguments
- `n_sobol_points::Int64`: Number of Sobol sequence points for global search
- `local_algorithm`: NLopt algorithm for local optimization
- `maxeval::Int64`: Maximum number of function evaluations for local optimization
- `xtol_rel::Float64`: Relative tolerance on parameter changes for local optimization
- `xtol_abs::Float64`: Absolute tolerance on parameter changes for local optimization
- `atol::Float64`: Absolute difference tolerance threshold; program errors if not met
"""
function NoiseVaccinationOptimizationParameters(;
        n_sobol_points::Int64 = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval::Int64 = 1000,
        xtol_rel::Float64 = 1.0e-3,
        xtol_abs::Float64 = 1.0e-3,
        atol::Float64 = 1.0e-3
    )
    return NoiseVaccinationOptimizationParameters(
        n_sobol_points,
        local_algorithm,
        maxeval,
        xtol_rel,
        xtol_abs,
        atol
    )
end
