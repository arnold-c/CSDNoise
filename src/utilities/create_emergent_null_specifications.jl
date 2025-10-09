export create_ensemble_specifications

function create_ensemble_specifications(
        time_specification,
        state_specification,
        target_disease_dynamics_params::TargetDiseaseDynamicsParameters,
        common_disease_dynamics_params::CommonDiseaseDynamicsParameters,
        dynamical_noise_specification::DynamicalNoiseSpecification
    )

    mu = calculate_mu(common_disease_dynamics_params.births_per_k_pop)
    sigma = 1 / target_disease_dynamics_params.latent_period_days
    gamma = 1 / target_disease_dynamics_params.infectious_duration_days

    beta_mean = calculate_beta(
        target_disease_dynamics_params.R_0,
        sigma,
        gamma,
        mu,
    )

    epsilon = calculate_import_rate(
        mu,
        target_disease_dynamics_params.R_0,
        state_specification.init_states.N
    )

    min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
        time_specification,
        state_specification,
        common_disease_dynamics_params.burnin_target_Reff,
        target_disease_dynamics_params.R_0,
        mu,
    )


    emergent_dynamics_specification = DynamicsParameterSpecification(
        beta_mean,
        target_disease_dynamics_params.beta_force,
        target_disease_dynamics_params.seasonality,
        sigma,
        gamma,
        mu,
        common_disease_dynamics_params.births_per_k_pop,
        epsilon,
        target_disease_dynamics_params.R_0,
        min_burnin_vaccination_coverage,
        target_disease_dynamics_params.max_burnin_vaccination_coverage,
        target_disease_dynamics_params.min_post_burnin_vaccination_coverage,
        target_disease_dynamics_params.max_post_burnin_vaccination_coverage,
    )

    null_dynamics_specification = DynamicsParameterSpecification(
        map(
            pn -> getproperty(emergent_dynamics_specification, pn),
            filter(
                name ->
                name != :min_post_burnin_vaccination_coverage &&
                    name != :max_post_burnin_vaccination_coverage,
                propertynames(emergent_dynamics_specification),
            ),
        )...,
        nothing,
        nothing,
    )

    ensemble_specification = EnsembleSpecification(
        state_specification,
        time_specification,
        emergent_dynamics_specification,
        null_dynamics_specification,
        dynamical_noise_specification,
        common_disease_dynamics_params.nsims,
    )

    return ensemble_specification
end
