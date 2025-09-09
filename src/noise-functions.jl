# module NoiseFunctions
#
# export create_noise_arr, create_noise_arr!

# include("ensemble-functions.jl")
# using .EnsembleFunctions
using UnPack

function create_noise_arr(
        noise_specification::DynamicalNoiseSpecification,
        incarr;
        ensemble_specification::EnsembleSpecification,
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
        if dynamics_parameter_specification.seasonality == cos
            sin
        elseif dynamics_parameter_specification.seasonality == sin
            cos
        end
    else
        dynamics_parameter_specification.seasonality
    end

    noise_dynamics_parameters = Vector{DynamicsParameters}(undef, nsims)

    noise_dynamics_parameter_specification = DynamicsParameterSpecification(
        dynamics_parameter_specification.beta_mean,
        noise_beta_force,
        noise_seasonality,
        1 / noise_specification.latent_period,
        1 / noise_specification.duration_infection,
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

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states), 2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_inc_vecs = Array{typeof(SVector(0)), 2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_beta_arr = zeros(Float64, tlength)

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        noise_dynamics_parameters[sim] = DynamicsParameters(
            noise_dynamics_parameter_specification; seed = run_seed
        )

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            noise_dynamics_parameters[sim],
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_inc_arr = zeros(
        Int64, size(ensemble_inc_vecs, 1), size(ensemble_inc_vecs, 2)
    )

    for sim in axes(ensemble_inc_vecs, 2)
        convert_svec_to_matrix!(
            @view(ensemble_inc_arr[:, sim]),
            @view(ensemble_inc_vecs[:, sim])
        )
    end

    poisson_noise = zeros(
        Float64, size(ensemble_inc_arr, 1), size(ensemble_inc_arr, 3)
    )

    add_poisson_noise_arr!(
        poisson_noise, ensemble_inc_arr, noise_specification.noise_mean_scaling;
        seed = seed,
    )

    mean_poisson_noise = NaNMath.mean(poisson_noise)
    mean_rubella_noise = StatsBase.mean(ensemble_inc_arr)

    poisson_noise_prop = mean_poisson_noise / mean_rubella_noise

    ensemble_inc_arr .+= poisson_noise

    return ensemble_inc_arr,
        (;
            mean_noise = mean_rubella_noise + mean_poisson_noise,
            mean_poisson_noise = mean_poisson_noise,
            mean_rubella_noise = mean_rubella_noise,
            poisson_noise_prop = poisson_noise_prop,
        )
end

function create_noise_arr(
        noise_specification::PoissonNoiseSpecification,
        incarr;
        seed = 1234,
        kwargs...,
    )
    noise_arr = zeros(Int64, size(incarr, 1), size(incarr, 3))

    add_poisson_noise_arr!(
        noise_arr, incarr, noise_specification.noise_mean_scaling; seed = seed
    )

    noise_rubella_prop = ones(Float64, size(incarr, 1), size(incarr, 3))
    return noise_arr, noise_rubella_prop
end

function add_poisson_noise_arr!(
        noise_arr, incarr, noise_mean_scaling; seed = 1234
    )
    Random.seed!(seed)

    @assert size(incarr, 3) == size(noise_arr, 2)
    return @inbounds for sim in axes(incarr, 3)
        @views noise_arr[:, sim] += rand(
            Poisson(noise_mean_scaling * mean(incarr[:, 1, sim])),
            size(incarr, 1),
        )
    end
end

# end
