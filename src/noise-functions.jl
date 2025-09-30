# module NoiseFunctions
#
# export create_noise_arr, create_noise_arr!

# include("ensemble-functions.jl")
# using .EnsembleFunctions
using UnPack
using LightSumTypes: variant
using StaticArrays: SVector
using MultistartOptimization
using NLopt
using StatsBase: mean
using StructArrays: StructVector
using Bumper: @no_escape, @alloc
using Random: Random
using Distributions: Poisson

"""
    create_noise_vecs(
        noise_specification::NoiseSpecification,
        ensemble_specification::EnsembleSpecification,
        endpoints::Vector{Int64};
        kwargs...
    )

Create a NoiseRun struct that contains a vector of noise vectors with variable lengths matching the provided endpoints, and the summary statistics.

This function generates noise for each simulation up to its specific endpoint,
allowing for variable-length noise generation based on EWS endpoint calculations.

# Arguments
- `noise_specification`: Specification for noise generation
- `ensemble_specification`: Ensemble simulation parameters
- `endpoints`: Vector of endpoints, one per simulation
- `kwargs...`: Additional keyword arguments

# Returns
- `NoiseRun`: Noise run with variable-length incidence vectors

# Example
```julia
noise_result = create_noise_vecs(
    noise_spec,
    ensemble_spec,
	enddates
)
```
"""
function create_noise_vecs(
        noise_specification::NoiseSpecification,
        varargs...;
        kwargs...
    )
    return create_noise_vecs(
        variant(noise_specification),
        varargs...;
        kwargs...
    )
end

function create_noise_vecs(
        noise_specification::DynamicalNoise,
        ensemble_specification::EnsembleSpecification,
        enddates::Vector{Int64};
        seed = 1234,
        kwargs...,
    )
    seed *= 10

    @unpack state_parameters,
        dynamics_parameter_specification, time_parameters,
        nsims =
        ensemble_specification
    @unpack tlength = time_parameters

    @assert nsims == length(enddates) "Number of simulations must match number of endpoints"

    @unpack init_states = state_parameters
    @unpack N = init_states
    @unpack noise_mean_scaling = noise_specification

    init_states_sv = SVector(init_states)

    noise_beta_force = if noise_specification.correlation == "none"
        0.0
    else
        dynamics_parameter_specification.beta_force
    end

    noise_seasonality = if noise_specification.correlation == "out-of-phase"
        if variant(dynamics_parameter_specification.seasonality) isa CosineSeasonality
            SeasonalityFunction(SineSeasonality())
        elseif variant(dynamics_parameter_specification.seasonality) isa SineSeasonality
            SeasonalityFunction(CosineSeasonality())
        else
            dynamics_parameter_specification.seasonality
        end
    else
        dynamics_parameter_specification.seasonality
    end

    noise_gamma = 1 / noise_specification.duration_infection
    noise_sigma = 1 / noise_specification.latent_period

    noise_beta_mean = calculate_beta(
        noise_specification.R_0,
        noise_gamma,
        dynamics_parameter_specification.mu,
        1,
        N
    )

    noise_dynamics_parameter_specification = DynamicsParameterSpecification(
        dynamics_parameter_specification.contact_matrix,
        noise_beta_mean,
        noise_beta_force,
        noise_seasonality,
        noise_sigma,
        noise_gamma,
        dynamics_parameter_specification.mu,
        dynamics_parameter_specification.annual_births_per_k,
        calculate_import_rate(
            dynamics_parameter_specification.mu,
            noise_specification.R_0,
            N
        ),
        noise_specification.R_0,
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
    )

    # Generate a single long vec of beta values that will be trimmed
    # for each simulation as beta doesn't have any stochasticity itself
    beta_vec = Vector{Float64}(undef, tlength)
    calculate_beta_amp!(
        beta_vec,
        dynamics_parameter_specification,
        time_parameters
    )

    # Initialize vectors for variable lengths
    mean_poisson_noise_vec = Vector{Float64}(undef, nsims)
    mean_dynamical_noise_vec = Vector{Float64}(undef, nsims)

    # Create vector of variable-sized vectors for incidence data
    incidence_vecs = Vector{Vector{Int64}}(undef, nsims)

    for sim in eachindex(incidence_vecs)
        run_seed = seed + (sim - 1)
        enddate = enddates[sim]
        Random.seed!(run_seed)

        local noise_dynamics_parameters = DynamicsParameters(
            noise_dynamics_parameter_specification; seed = run_seed
        )

        local beta_worker_vec = @view(beta_vec[1:enddate])

        _calculate_dynamic_noise_values!(
            incidence_vecs,
            mean_dynamical_noise_vec,
            mean_poisson_noise_vec,
            beta_worker_vec,
            init_states_sv,
            noise_dynamics_parameters,
            noise_mean_scaling,
            time_parameters,
            enddate,
            sim,
        )

    end

    mean_dynamical_noise = mean(mean_dynamical_noise_vec)
    mean_poisson_noise = mean(mean_poisson_noise_vec)
    mean_noise = mean_dynamical_noise + mean_poisson_noise

    return NoiseRun(
        incidence_vecs,
        mean_noise,
        mean_poisson_noise,
        mean_dynamical_noise
    )
end

function create_noise_vecs(
        noise_specification::PoissonNoise,
        ensemble_specification::EnsembleSpecification,
        enddates::Vector{Int64},
        seir_results::StructVector{SEIRRun};
        seed = 1234,
        kwargs...,
    )
    @unpack nsims = ensemble_specification
    @assert nsims == length(enddates) "Number of simulations must match number of endpoints"
    @assert nsims == length(seir_results) "Number of simulations must match SEIR results"

    @unpack noise_mean_scaling = noise_specification

    # Initialize vectors for variable lengths
    incidence_vecs = Vector{Vector{Int64}}(undef, nsims)
    mean_poisson_noise_vec = Vector{Float64}(undef, nsims)

    for sim in eachindex(mean_poisson_noise_vec)
        run_seed = seed + (sim - 1)
        enddate = enddates[sim]
        Random.seed!(run_seed)

        mean_dynamical_noise_incidence = mean(seir_results.incidence[sim])

        _calculate_poisson_noise_values!(
            incidence_vecs,
            mean_poisson_noise_vec,
            mean_dynamical_noise_incidence,
            noise_mean_scaling,
            enddate,
            sim,
        )

    end

    mean_noise = mean(mean_poisson_noise_vec)

    return NoiseRun(
        incidence_vecs,
        mean_noise,
        mean_noise,
        0.0
    )
end

function _calculate_dynamic_noise_values!(
        incidence_vecs,
        mean_dynamical_noise_vec,
        mean_poisson_noise_vec,
        beta_worker_vec,
        init_states_sv,
        noise_dynamics_parameters,
        noise_mean_scaling,
        time_parameters,
        enddate,
        sim,
    )
    @no_escape begin
        seir_worker_vec = @alloc(SVector{5, Int64}, enddate)
        incidence_worker_vec = @alloc(Int64, enddate)
        Reff_worker_vec = @alloc(Float64, enddate)

        seir_mod!(
            seir_worker_vec,
            incidence_worker_vec,
            Reff_worker_vec,
            beta_worker_vec,
            init_states_sv,
            noise_dynamics_parameters,
            time_parameters,
        )

        mean_dynamical_noise_incidence = mean(incidence_worker_vec)

        poisson_noise_worker_vec = @alloc(Int64, enddate)

        _add_poisson_noise!(
            poisson_noise_worker_vec,
            mean_dynamical_noise_incidence,
            noise_mean_scaling,
        )

        mean_poisson_noise = mean(poisson_noise_worker_vec)

        combined_noise_vec = Vector{Int64}(undef, enddate)
        @inbounds for i in eachindex(combined_noise_vec)
            combined_noise_vec[i] = incidence_worker_vec[i] + poisson_noise_worker_vec[i]
        end

        incidence_vecs[sim] = combined_noise_vec
        mean_dynamical_noise_vec[sim] = mean_dynamical_noise_incidence
        mean_poisson_noise_vec[sim] = mean_poisson_noise
    end
    return nothing
end

function _calculate_poisson_noise_values!(
        incidence_vecs,
        mean_poisson_noise_vec,
        mean_dynamical_noise_incidence,
        noise_mean_scaling,
        enddate,
        sim,
    )
    poisson_noise_vec = Vector{Int64}(undef, enddate)

    _add_poisson_noise!(
        poisson_noise_vec,
        mean_dynamical_noise_incidence,
        noise_mean_scaling,
    )

    mean_poisson_noise = mean(poisson_noise_vec)

    incidence_vecs[sim] = poisson_noise_vec
    mean_poisson_noise_vec[sim] = mean_poisson_noise
    return nothing
end

function _add_poisson_noise!(
        noise_vec::AIV,
        mean_dynamical_noise_incidence::Float64,
        noise_mean_scaling::Float64,
    ) where {AIV <: AbstractVector{<:Integer}}
    poisson_rate = noise_mean_scaling * mean_dynamical_noise_incidence
    @inbounds @simd for i in eachindex(noise_vec)
        noise_vec[i] = rand(Poisson(poisson_rate))
    end

    return nothing
end
