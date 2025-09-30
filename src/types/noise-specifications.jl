using LightSumTypes: @sumtype
using NLopt: NLopt

export NoiseSpecification,
    PoissonNoise,
    DynamicalNoise,
    DynamicalNoiseSpecification,
    NoiseVaccinationOptimizationParameters

abstract type AbstractNoiseSpecification end

struct PoissonNoise
    noise_mean_scaling::Float64
end

struct DynamicalNoise
    R_0::Float64
    latent_period::Int64
    duration_infection::Int64
    correlation::String
    noise_mean_scaling::Float64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
end

@sumtype NoiseSpecification(PoissonNoise, DynamicalNoise) <: AbstractNoiseSpecification


struct DynamicalNoiseSpecification
    R0::Float64
    latent_period::Int64
    duration_infection::Int64
    correlation::String
    poisson_component::Float64
    vaccination_bounds::Vector{Float64}
    susceptible_bounds::Vector{Float64}
    max_vaccination_range::Float64
    function DynamicalNoiseSpecification(
            R0::Float64,
            latent_period::Int64,
            duration_infection::Int64,
            correlation::String,
            poisson_component::Float64,
            vaccination_bounds::Vector{Float64},
            susceptible_bounds::Vector{Float64},
            max_vaccination_range::Float64,
        )

        @assert length(vaccination_bounds) == 2
        @assert vaccination_bounds[1] < vaccination_bounds[2]
        @assert length(susceptible_bounds) == 2
        @assert susceptible_bounds[1] < susceptible_bounds[2]
        return new(
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            vaccination_bounds,
            susceptible_bounds,
            max_vaccination_range

        )
    end
end

function DynamicalNoiseSpecification(;
        R0::Float64,
        latent_period::Int64,
        duration_infection::Int64,
        correlation::String,
        poisson_component::Float64,
        vaccination_bounds::Vector{Float64} = [0.0, 1.0],
        susceptible_bounds::Vector{Float64} = [0.01, 0.99],
        max_vaccination_range::Float64 = 0.2
    )
    return DynamicalNoiseSpecification(
        R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        vaccination_bounds,
        susceptible_bounds,
        max_vaccination_range
    )
end

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
