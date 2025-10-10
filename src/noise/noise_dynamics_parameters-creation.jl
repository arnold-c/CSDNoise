function recreate_noise_dynamics_spec(
        noise_specification::DynamicalNoise,
        ensemble_specification::EnsembleSpecification
    )
    @unpack state_parameters,
        emergent_dynamics_parameter_specification = ensemble_specification

    N = state_parameters.init_states.N

    noise_beta_force = if noise_specification.correlation == "none"
        0.0
    else
        emergent_dynamics_parameter_specification.beta_force
    end

    noise_seasonality = if noise_specification.correlation == "out-of-phase"
        if LightSumTypes.variant(emergent_dynamics_parameter_specification.seasonality) isa CosineSeasonality
            SeasonalityFunction(SineSeasonality())
        elseif LightSumTypes.variant(emergent_dynamics_parameter_specification.seasonality) isa SineSeasonality
            SeasonalityFunction(CosineSeasonality())
        else
            emergent_dynamics_parameter_specification.seasonality
        end
    else
        emergent_dynamics_parameter_specification.seasonality
    end

    noise_gamma = 1 / noise_specification.duration_infection
    noise_sigma = 1 / noise_specification.latent_period

    noise_beta_mean = calculate_beta(
        noise_specification.R_0,
        noise_sigma,
        noise_gamma,
        emergent_dynamics_parameter_specification.mu
    )

    noise_dynamics_parameter_specification = DynamicsParameterSpecification(
        beta_mean = noise_beta_mean,
        beta_force = noise_beta_force,
        seasonality = noise_seasonality,
        sigma = noise_sigma,
        gamma = noise_gamma,
        mu = emergent_dynamics_parameter_specification.mu,
        annual_births_per_k = emergent_dynamics_parameter_specification.annual_births_per_k,
        epsilon = calculate_import_rate(
            emergent_dynamics_parameter_specification.mu,
            noise_specification.R_0,
            N
        ),
        R_0 = noise_specification.R_0,
        min_burnin_vaccination_coverage = noise_specification.vaccination_coverage,
        max_burnin_vaccination_coverage = noise_specification.vaccination_coverage,
        min_post_burnin_vaccination_coverage = noise_specification.vaccination_coverage,
        max_post_burnin_vaccination_coverage = noise_specification.vaccination_coverage,
    )

    return noise_dynamics_parameter_specification
end
