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

function create_noise_arr(
        noise_specification::NoiseSpecification,
        args...;
        kwargs...
    )
    return create_noise_arr(variant(noise_specification), args...; kwargs...)
end

function create_noise_arr(
        noise_specification::DynamicalNoise,
        ensemble_specification::EnsembleSpecification,
        args...;
        seed = 1234,
        kwargs...,
    )
    seed *= 10

    @unpack state_parameters,
        dynamics_parameter_specification, time_parameters,
        nsims =
        ensemble_specification
    @unpack tlength = time_parameters

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
        state_parameters.init_states.N,
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
            state_parameters.init_states.N,
        ),
        noise_specification.R_0,
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
    )

    seir_vec = FixedSizeVector{SVector{5, Int64}}(undef, tlength)
    Reff_vec = FixedSizeVector{Float64}(undef, tlength)
    beta_vec = FixedSizeVector{Float64}(undef, tlength)
    calculate_beta_amp!(
        beta_vec,
        dynamics_parameter_specification,
        time_parameters
    )

    # Initialize vectors
    poisson_noise_worker = FixedSizedVector{Int64}(undef, tlength)
    mean_poisson_noise_vec = FixedSizedVector{Float64}(undef, nsims)
    mean_dynamical_noise_vec = FixedSizedVector{Float64}(undef, nsims)

    # Create vector of sized vectors for incidence data and apply Poisson noise in-place
    incidence_vecs = FixedSizeVector{FixedSizedVector{Int64}}(undef, nsims)
    incidence_worker_vec = similar(poisson_noise_worker)

    for sim in eachindex(incidence_vecs)
        run_seed = seed + (sim - 1)

        local noise_dynamics_parameters = DynamicsParameters(
            noise_dynamics_parameter_specification; seed = run_seed
        )

        seir_mod!(
            seir_vec,
            incidence_worker_vec,
            Reff_vec,
            beta_vec,
            SVector(state_parameters.init_states),
            noise_dynamics_parameters,
            time_parameters,
            run_seed,
        )

        local mean_dynamical_noise_incidence = mean(incidence_worker_vec)

        mean_dynamical_noise_vec[sim] = mean_dynamical_noise_incidence

        add_poisson_noise!(
            poisson_noise_worker,
            mean_dynamical_noise_incidence,
            noise_specification.noise_mean_scaling,
            run_seed
        )

        mean_poisson_noise_vec[sim] = mean(poisson_noise_worker)
        incidence_vecs[sim] = incidence_worker_vec .+ poisson_noise_worker
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

function create_noise_arr(
        noise_specification::PoissonNoise,
        ensemble_specification::EnsembleSpecification,
        seir_results::Vector{SEIRRun};
        seed = 1234,
        kwargs...,
    )
    # Create vector of sized vectors for incidence data and apply Poisson noise in-place
    tlength = ensemble_specification.time_parameters.tlength
    nsims = ensemble_specification.nsims
    noise_vecs = FixedSizeVector{FixedSizedVector{Int64}}(undef, nsims)
    noise_worker_vec = FixedSizeVector{Int64}(undef, tlength)
    mean_noise_vec = FixedSizeVector{Float64}(undef, nsims)
    for sim in eachindex(noise_vecs)
        run_seed = seed + (sim - 1)
        mean_dynamical_noise_incidence = mean(seir_results[sim].incidence)
        mean_noise_vec[sim] = mean_dynamical_noise_incidence

        add_poisson_noise!(
            noise_worker_vec,
            mean_dynamical_noise_incidence,
            noise_specification.noise_mean_scaling,
            run_seed
        )

        noise_vecs[sim] = noise_worker_vec
    end
    mean_noise = mean(mean_noise_vec)

    return NoiseRun(
        noise_vecs,
        mean_noise,
        mean_noise,
        0.0
    )
end

function add_poisson_noise!(
        noise_vec::AIV,
        mean_dynamical_noise_incidence::Float64,
        noise_mean_scaling::Float64,
        seed = 1234
    ) where {AIV <: AbstractVector{<:Integer}}

    noise_vec .= rand(
        Poisson(noise_mean_scaling * mean_dynamical_noise_incidence),
        length(noise_vec)
    )

    return nothing
end
