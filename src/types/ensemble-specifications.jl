export EnsembleSpecification,
    EnsembleSpecsParameters

Base.@kwdef struct EnsembleSpecification
    state_parameters::StateParameters
    time_parameters::SimTimeParameters
    emergent_dynamics_parameter_specification::DynamicsParameterSpecification
    null_dynamics_parameter_specification::DynamicsParameterSpecification
    dynamical_noise_specification::DynamicalNoiseSpecification
    nsims::Int64
    dirpath::String
end

function EnsembleSpecification(
        state_parameters::StateParameters,
        time_parameters::SimTimeParameters,
        emergent_dynamics_parameter_specification::DynamicsParameterSpecification,
        null_dynamics_parameter_specification::DynamicsParameterSpecification,
        dynamical_noise_specification::DynamicalNoiseSpecification,
        nsims::Int64,
    )
    dirpath = outdir(
        "ensemble",
        "seasonal-infectivity-import",
        "tau-leaping",
        "N_$(state_parameters.init_states.N)",
        "r_$(state_parameters.init_state_props.r_prop)",
        "nsims_$(nsims)",
        "R0_$(emergent_dynamics_parameter_specification.R_0)",
        "latent_period_$(round(1 / emergent_dynamics_parameter_specification.sigma; digits = 2))",
        "infectious_period_$(round(1 / emergent_dynamics_parameter_specification.gamma; digits = 2))",
        "noise_R0_$(dynamical_noise_specification.R_0)",
        "noise_latent_period_$(round(dynamical_noise_specification.latent_period; digits = 2))",
        "noise_infectious_period_$(round(dynamical_noise_specification.duration_infection; digits = 2))",
        "min_burnin_vaccination_coverage_$(emergent_dynamics_parameter_specification.min_burnin_vaccination_coverage)",
        "max_burnin_vaccination_coverage_$(emergent_dynamics_parameter_specification.max_burnin_vaccination_coverage)",
        "min_post_burnin_vaccination_coverage_$(emergent_dynamics_parameter_specification.min_post_burnin_vaccination_coverage)",
        "max_post_burnin_vaccination_coverage_$(emergent_dynamics_parameter_specification.max_post_burnin_vaccination_coverage)",
        "births_per_k_$(emergent_dynamics_parameter_specification.annual_births_per_k)",
        "beta_force_$(emergent_dynamics_parameter_specification.beta_force)",
        "burnin_$(time_parameters.burnin)",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)",
    )

    return EnsembleSpecification(
        state_parameters,
        time_parameters,
        emergent_dynamics_parameter_specification,
        null_dynamics_parameter_specification,
        dynamical_noise_specification,
        nsims,
        dirpath,
    )
end

# TODO: Figure out if this is actually needed
struct EnsembleSpecsParameters
    burnin_years::Int64
    tmax_years::Int64
    annual_births_per_k::Int64
    ensemble_state_specification::StateParameters
    R_0::Float64
    gamma::Float64
    sigma::Float64
    target_Reff::Float64
    target_years::Int64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
    nsims::Int64
    function EnsembleSpecsParameters(
            burnin_years,
            tmax_years,
            annual_births_per_k,
            ensemble_state_specification,
            R_0,
            gamma,
            sigma,
            target_Reff,
            target_years,
            min_vaccination_coverage,
            max_vaccination_coverage,
            nsims,
        )
        @assert tmax_years >= target_years
        @assert min_vaccination_coverage < max_vaccination_coverage

        return new(
            burnin_years,
            tmax_years,
            annual_births_per_k,
            ensemble_state_specification,
            R_0,
            gamma,
            sigma,
            target_Reff,
            target_years,
            min_vaccination_coverage,
            max_vaccination_coverage,
            nsims
        )
    end

end

function EnsembleSpecsParameters(;
        burnin_years::Int,
        tmax_years::Int,
        annual_births_per_k::Int64 = ANNUAL_BIRTHS_PER_K,
        ensemble_state_specification::StateParameters = StateParameters(
            500_000,
            Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
        ),
        R_0::Float64 = R0,
        gamma::Float64 = GAMMA,
        sigma::Float64 = SIGMA,
        target_Reff::Float64 = 0.9,
        target_years::Int = 2 * burnin_years,
        min_vaccination_coverage::Float64 = 0.6,
        max_vaccination_coverage::Float64 = 0.8,
        nsims::Int = 1000
    )

    return EnsembleSpecsParameters(
        burnin_years,
        tmax_years,
        annual_births_per_k,
        ensemble_state_specification,
        R_0,
        gamma,
        sigma,
        target_Reff,
        target_years,
        min_vaccination_coverage,
        max_vaccination_coverage,
        nsims
    )
end
