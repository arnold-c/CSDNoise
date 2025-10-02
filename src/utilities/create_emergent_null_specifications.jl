export create_ensemble_specifications

# struct EnsembleSpecsParameters
#     burnin_years::Int64
#     nyears::Int64
#     annual_births_per_k::Int64
#     ensemble_state_specification::StateParameters
#     R_0::Float64
#     gamma::Float64
#     sigma::Float64
#     target_Reff::Float64
#     target_years::Int64
#     min_vaccination_coverage::Float64
#     max_vaccination_coverage::Float64
#     nsims::Int64
# end
# create_ensemble_specs(
#     EnsembleSpecsParameters(
#         burnin_years = 5,
#         tmax_years = 20,
#         annual_births_per_k = ANNUAL_BIRTHS_PER_K,
#         ensemble_state_specification = StateParameters(
#             500_000,
#             Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
#         ),
#         R_0 = R0,
#         gamma = GAMMA,
#         sigma = SIGMA,
#         target_Reff = 0.9,
#         target_years = 10,
#         min_vaccination_coverage = 0.6,
#         max_vaccination_coverage = 0.8,
#         nsims = 100
#     )
# )
function create_ensemble_specifications(
        time_specification,
        state_specification,
        target_disease_dynamics_params::TargetDiseaseDynamicsParameters,
        common_disease_dynamics_params::CommonDiseaseDynamicsParameters,
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

    emergent_specification = EnsembleSpecification(
        state_specification,
        emergent_dynamics_specification,
        time_specification,
        common_disease_dynamics_params.nsims,
    )

    null_specification = EnsembleSpecification(
        state_specification,
        null_dynamics_specification,
        time_specification,
        common_disease_dynamics_params.nsims,
    )

    return (;
        emergent_specification,
        null_specification,
    )
end
