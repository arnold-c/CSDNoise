export create_noise_vecs

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
        LightSumTypes.variant(noise_specification),
        varargs...;
        kwargs...
    )
end

function create_noise_vecs(
        noise_specification::DynamicalNoise,
        ensemble_specification::EnsembleSpecification,
        noise_dynamics_parameter_specification::DynamicsParameterSpecification,
        endemic_props_result::Union{Try.Ok, Try.Err},
        enddates_vec::Vector{Int64};
        seed = 1234,
        kwargs...,
    )
    seed *= 10

    @unpack state_parameters,
        time_parameters,
        nsims = ensemble_specification
    @unpack tlength = time_parameters
    @unpack init_states = state_parameters
    @unpack N = init_states

    @assert nsims == length(enddates_vec) "Number of simulations must match number of endpoints"

    @unpack noise_mean_scaling = noise_specification

    init_states_sv = StaticArrays.SVector(init_states)


    # Generate a single long vec of beta values that will be trimmed
    # for each simulation as beta doesn't have any stochasticity itself
    beta_vec = Vector{Float64}(undef, tlength)
    calculate_beta_amp!(
        beta_vec,
        noise_dynamics_parameter_specification,
        time_parameters
    )

    # Initialize vectors for variable lengths
    mean_poisson_noise_vec = Vector{Float64}(undef, nsims)
    mean_dynamical_noise_vec = Vector{Float64}(undef, nsims)

    # Create vector of variable-sized vectors for incidence data
    incidence_vecs = Vector{Vector{Int64}}(undef, nsims)

    for sim in eachindex(incidence_vecs)
        run_seed = seed + (sim - 1)
        enddate = enddates_vec[sim]
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

    mean_dynamical_noise = StatsBase.mean(mean_dynamical_noise_vec)
    mean_poisson_noise = StatsBase.mean(mean_poisson_noise_vec)
    mean_noise = mean_dynamical_noise + mean_poisson_noise

    return NoiseRun(
        incidence = incidence_vecs,
        mean_noise = mean_noise,
        mean_poisson_noise = mean_poisson_noise,
        mean_dynamic_noise = mean_dynamical_noise
    )
end

function create_noise_vecs(
        noise_specification::PoissonNoise,
        ensemble_specification::EnsembleSpecification,
        enddates_vec::Vector{Int64},
        seir_results::StructVector{SEIRRun};
        seed = 1234,
        kwargs...,
    )
    @unpack nsims = ensemble_specification
    @assert nsims == length(enddates_vec) "Number of simulations must match number of endpoints"
    @assert nsims == length(seir_results) "Number of simulations must match SEIR results"

    @unpack noise_mean_scaling = noise_specification

    # Initialize vectors for variable lengths
    incidence_vecs = Vector{Vector{Int64}}(undef, nsims)
    mean_poisson_noise_vec = Vector{Float64}(undef, nsims)

    mean_incidence = calculate_mean_incidence(seir_results)

    for sim in eachindex(incidence_vecs)
        run_seed = seed + (sim - 1)
        Random.seed!(run_seed)

        enddate = enddates_vec[sim]

        @assert length(seir_results.incidence[sim]) == enddate

        _calculate_poisson_noise_values!(
            incidence_vecs,
            mean_poisson_noise_vec,
            mean_incidence,
            noise_mean_scaling,
            enddate,
            sim,
        )

    end

    mean_noise = StatsBase.mean(mean_poisson_noise_vec)
    @assert isapprox(mean_incidence, mean_noise; atol = 1.0e-3)

    return NoiseRun(
        incidence = incidence_vecs,
        mean_noise = mean_noise,
        mean_poisson_noise = mean_noise,
        mean_dynamic_noise = 0.0
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
        seir_worker_vec = @alloc(StaticArrays.SVector{5, Int64}, enddate)
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

        mean_dynamical_noise_incidence = StatsBase.mean(incidence_worker_vec)

        poisson_noise_worker_vec = @alloc(Int64, enddate)

        _add_poisson_noise!(
            poisson_noise_worker_vec,
            mean_dynamical_noise_incidence,
            noise_mean_scaling,
        )

        mean_poisson_noise = StatsBase.mean(poisson_noise_worker_vec)

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
        mean_incidence,
        noise_mean_scaling,
        enddate,
        sim,
    )
    poisson_noise_vec = Vector{Int64}(undef, enddate)

    _add_poisson_noise!(
        poisson_noise_vec,
        mean_incidence,
        noise_mean_scaling,
    )

    incidence_vecs[sim] = poisson_noise_vec
    mean_poisson_noise_vec[sim] = StatsBase.mean(poisson_noise_vec)
    return nothing
end

function _add_poisson_noise!(
        noise_vec::AIV,
        mean_incidence::Float64,
        noise_mean_scaling::Float64,
    ) where {AIV <: AbstractVector{<:Integer}}
    poisson_rate = noise_mean_scaling * mean_incidence
    @inbounds @simd for i in eachindex(noise_vec)
        noise_vec[i] = rand(Distributions.Poisson(poisson_rate))
    end

    return nothing
end
