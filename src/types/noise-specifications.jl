export NoiseType,
    PoissonNoiseType,
    DynamicalNoiseType,
    NoiseSpecification,
    PoissonNoise,
    DynamicalNoise,
    DynamicalNoiseSpecification,
    NoiseVaccinationOptimizationParameters

"""
    DynamicalNoiseSpecification

Configuration for dynamical noise parameters in epidemiological simulations.

This struct defines the parameters needed to generate dynamical noise that incorporates
epidemiological dynamics, including disease transmission parameters, correlation structure,
and vaccination bounds for optimization.

# Fields
- `R_0::Float64`: Basic reproduction number for the disease
- `latent_period::Int64`: Duration of the latent period in days
- `duration_infection::Int64`: Duration of the infectious period in days
- `correlation::String`: Type of correlation structure for noise generation
- `poisson_component::Float64`: Scaling factor for the Poisson noise component
- `vaccination_bounds::Vector{Float64}`: Lower and upper bounds for vaccination coverage optimization

# Constructor
    DynamicalNoiseSpecification(R_0, latent_period, duration_infection, correlation, poisson_component, vaccination_bounds)

The constructor validates that `vaccination_bounds` has exactly 2 elements and that the lower bound
is less than the upper bound.

# Example
```julia
spec = DynamicalNoiseSpecification(
    R_0 = 2.5,
    latent_period = 4,
    duration_infection = 6,
    correlation = "in-phase",
    poisson_component = 0.1,
    vaccination_bounds = [0.0, 0.8]
)
```
"""
abstract type AbstractNoiseType end
struct PoissonNoiseType end
struct DynamicalNoiseType end

LightSumTypes.@sumtype NoiseType(PoissonNoiseType, DynamicalNoiseType) <: AbstractNoiseType

struct DynamicalNoiseSpecification
    R_0::Float64
    latent_period::Float64
    duration_infection::Float64
    correlation::String
    poisson_component::Float64
    vaccination_bounds::Vector{Float64}
    function DynamicalNoiseSpecification(
            R_0::Float64,
            latent_period::Float64,
            duration_infection::Float64,
            correlation::String,
            poisson_component::Float64,
            vaccination_bounds::Vector{Float64},
        )

        @assert length(vaccination_bounds) == 2
        @assert vaccination_bounds[1] >= 0
        @assert vaccination_bounds[2] <= 1
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

"""
    DynamicalNoiseSpecification(; R_0, latent_period, duration_infection, correlation, poisson_component, vaccination_bounds = [0.0, 1.0])

Keyword constructor for `DynamicalNoiseSpecification` with default vaccination bounds.

# Arguments
- `R_0::Float64`: Basic reproduction number for the disease
- `latent_period::Int64`: Duration of the latent period in days
- `duration_infection::Int64`: Duration of the infectious period in days
- `correlation::String`: Type of correlation structure for noise generation
- `poisson_component::Float64`: Scaling factor for the Poisson noise component
- `vaccination_bounds::Vector{Float64}`: Lower and upper bounds for vaccination coverage (default: [0.0, 1.0])

# Returns
- `DynamicalNoiseSpecification`: Configured dynamical noise specification

# Example
```julia
spec = DynamicalNoiseSpecification(
    R_0 = 2.5,
    latent_period = 4,
    duration_infection = 6,
    correlation = "exponential",
    poisson_component = 0.1
)
```
"""
function DynamicalNoiseSpecification(;
        R_0::Float64,
        latent_period::Float64,
        duration_infection::Float64,
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

"""
    AbstractNoiseSpecification

Abstract base type for all noise specification types in the CSDNoise package.

This type serves as the parent type for different noise specification implementations,
enabling polymorphic behavior and type constraints in function signatures.
"""
abstract type AbstractNoiseSpecification end

"""
    PoissonNoise

Simple Poisson noise specification for epidemiological simulations.

This struct defines parameters for generating Poisson-distributed noise that is added
to the mean incidence values in simulations. The noise is independent across time points
and follows a Poisson distribution scaled by the mean incidence.

# Fields
- `noise_mean_scaling::Float64`: Scaling factor applied to the mean incidence before Poisson sampling

# Constructor
    PoissonNoise(; noise_mean_scaling)

# Example
```julia
poisson_noise = PoissonNoise(noise_mean_scaling = 0.1)
```
"""
Base.@kwdef struct PoissonNoise
    noise_mean_scaling::Float64
end

"""
    DynamicalNoise

Dynamical noise specification with specific vaccination coverage.

This struct represents a fully specified dynamical noise configuration that includes
both the epidemiological parameters and a specific vaccination coverage level. It is
typically created from a `DynamicalNoiseSpecification` by specifying the vaccination
coverage within the allowed bounds.

# Fields
- `R_0::Float64`: Basic reproduction number for the disease
- `latent_period::Int64`: Duration of the latent period in days
- `duration_infection::Int64`: Duration of the infectious period in days
- `correlation::String`: Type of correlation structure for noise generation
- `noise_mean_scaling::Float64`: Scaling factor for the Poisson noise component
- `vaccination_coverage::Float64`: Specific vaccination coverage level (0.0 to 1.0)

# Constructor
    DynamicalNoise(; R_0, latent_period, duration_infection, correlation, noise_mean_scaling, vaccination_coverage)

# Example
```julia
dyn_noise = DynamicalNoise(
    R_0 = 2.5,
    latent_period = 4,
    duration_infection = 6,
    correlation = "in-phase",
    noise_mean_scaling = 0.1,
    vaccination_coverage = 0.6
)
```
"""
Base.@kwdef struct DynamicalNoise
    R_0::Float64
    latent_period::Float64
    duration_infection::Float64
    correlation::String
    poisson_component::Float64
    vaccination_coverage::Float64
end

"""
    DynamicalNoise(dynamical_noise_specification::DynamicalNoiseSpecification, vaccination_coverage)

Create a `DynamicalNoise` instance from a specification and vaccination coverage.

This constructor takes a `DynamicalNoiseSpecification` and a specific vaccination coverage
value to create a fully specified `DynamicalNoise` instance. The vaccination coverage must
be between 0.0 and 1.0.

# Arguments
- `dynamical_noise_specification::DynamicalNoiseSpecification`: The base specification containing epidemiological parameters
- `vaccination_coverage`: Specific vaccination coverage level (must be between 0.0 and 1.0)

# Returns
- `DynamicalNoise`: Fully specified dynamical noise configuration

# Example
```julia
spec = DynamicalNoiseSpecification(
    R_0 = 2.5,
    latent_period = 4,
    duration_infection = 6,
    correlation = "exponential",
    poisson_component = 0.1,
    vaccination_bounds = [0.0, 0.8]
)
dyn_noise = DynamicalNoise(spec, 0.6)
```
"""
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
        poisson_component = dynamical_noise_specification.poisson_component,
        vaccination_coverage = vaccination_coverage
    )
end

"""
    NoiseSpecification

Sum type representing either Poisson or dynamical noise specifications.

This sum type allows functions to accept either `PoissonNoise` or `DynamicalNoise`
specifications in a type-safe manner. It uses `LightSumTypes.jl` to create an
efficient tagged union that can be pattern matched.

# Variants
- `PoissonNoise`: Simple Poisson noise specification
- `DynamicalNoise`: Dynamical noise with epidemiological parameters

# Example
```julia
# Can hold either type
poisson_spec = NoiseSpecification(PoissonNoise(noise_mean_scaling = 0.1))
dyn_spec = NoiseSpecification(DynamicalNoise(R_0 = 2.5, ...))

# Pattern matching
function process_noise(spec::NoiseSpecification)
	process_noise(LightSumTypes.variant(spec))
end
process_noise(spec::PoissonNoise) = ...
process_noise(spec::DynamicalNoise) = ...
```
"""
LightSumTypes.@sumtype NoiseSpecification(PoissonNoise, DynamicalNoise) <: AbstractNoiseSpecification

"""
    NoiseVaccinationOptimizationParameters

Parameters for optimizing vaccination coverage in dynamical noise scenarios.

This struct contains all the optimization parameters needed for finding optimal
vaccination coverage levels to achieve a target mean noise level when using
dynamical noise. It combines global search using Sobol sequences with local
optimization using NLopt algorithms.

# Fields
- `n_sobol_points::Int64`: Number of Sobol sequence points for global search initialization
- `local_algorithm`: NLopt algorithm used for local optimization (e.g., NLopt.LN_BOBYQA)
- `maxeval::Int64`: Maximum number of function evaluations allowed for local optimization
- `xtol_rel::Float64`: Relative tolerance on parameter changes for convergence
- `xtol_abs::Float64`: Absolute tolerance on parameter changes for convergence
- `atol::Float64`: Absolute difference tolerance threshold for optimization success

# Constructor
    NoiseVaccinationOptimizationParameters(n_sobol_points, local_algorithm, maxeval, xtol_rel, xtol_abs, atol)

# Example
```julia
opt_params = NoiseVaccinationOptimizationParameters(
    n_sobol_points = 200,
    local_algorithm = NLopt.LN_BOBYQA,
    maxeval = 2000,
    xtol_rel = 1e-4,
    xtol_abs = 1e-4,
    atol = 1e-5
)
```
"""
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

Parameters for optimizing noise vaccination levels for a given mean noise level.

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
