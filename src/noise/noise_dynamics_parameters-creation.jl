function recreate_noise_dynamics_spec(
        noise_specification::DynamicalNoise,
        ensemble_specification::EnsembleSpecification
    )
    @unpack state_parameters,
        dynamics_parameter_specification = ensemble_specification

    N = state_parameters.init_states.N

    noise_beta_force = if noise_specification.correlation == "none"
        0.0
    else
        dynamics_parameter_specification.beta_force
    end

    noise_seasonality = if noise_specification.correlation == "out-of-phase"
        if LightSumTypes.variant(dynamics_parameter_specification.seasonality) isa CosineSeasonality
            SeasonalityFunction(SineSeasonality())
        elseif LightSumTypes.variant(dynamics_parameter_specification.seasonality) isa SineSeasonality
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
        noise_sigma,
        noise_gamma,
        dynamics_parameter_specification.mu
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
        noise_specification.vaccination_coverage,
        noise_specification.vaccination_coverage,
        noise_specification.vaccination_coverage,
        noise_specification.vaccination_coverage,
    )

    return noise_dynamics_parameter_specification
end
