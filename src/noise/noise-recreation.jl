export recreate_noise_vecs

"""
    recreate_noise_vecs(
        dynamical_noise_spec::DynamicalNoiseSpecification,
        mean_vaccination_coverage,
        susceptible_proportion,
        ensemble_specification::EnsembleSpecification,
        enddates_vec,
        dynamics_parameter_specification
    )

Recreate noise vectors with updated vaccination coverage and susceptible proportion parameters.

Takes updated susceptible proportion and vaccination coverage values (either during or after
optimization), recreates the state parameters, ensemble specification, and dynamics
specifications, then generates new noise vectors using the updated parameters.

If `dynamics_parameter_specification` is provided, calculates endemic equilibrium proportions
for E, I, and R compartments based on the vaccination coverage. Otherwise, uses the original
ensemble specification's initial state proportions.

# Arguments
- `dynamical_noise_spec`: DynamicalNoiseSpecification containing fixed noise parameters
- `mean_vaccination_coverage`: Mean vaccination coverage level (0-1) - optimization variable
- `susceptible_proportion`: Initial proportion of population that is susceptible (0-1) - optimization variable
- `ensemble_specification`: Parameters for ensemble simulation
- `enddates_vec`: Vector of simulation end dates
- `dynamics_parameter_specification`: Optional DynamicsParameterSpecification for endemic equilibrium calculation

# Returns
- `NoiseRun`: StructVector containing the generated noise simulation results

# Example
```julia
noise_spec = DynamicalNoiseSpecification(
    R0 = 12.0,
    latent_period = 8,
    duration_infection = 7,
    correlation = "in-phase",
    poisson_component = 1.0
)
noise_result = recreate_noise_vecs(
    noise_spec,
    0.85,     # mean_vaccination_coverage
    0.7,      # susceptible_proportion
    ensemble_spec,
    enddates_vec,
    dynamics_param_spec
)
```
"""
function recreate_noise_vecs(
        dynamical_noise_spec::DynamicalNoiseSpecification,
        mean_vaccination_coverage,
        ensemble_specification::EnsembleSpecification,
        enddates_vec;
        verbose = false
    )
    @unpack R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        max_vaccination_range = dynamical_noise_spec

    # Create final EnsembleSpecification with optimal parameters for verification
    @unpack state_parameters,
        dynamics_parameter_specification,
        time_parameters,
        nsims,
        dirpath = ensemble_specification
    @unpack init_states, init_state_props = state_parameters
    @unpack N = init_states

    min_vaccination_coverage, max_vaccination_coverage = calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range,
    )

    updated_dynamical_noise_spec = DynamicalNoise(
        R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        min_vaccination_coverage,
        max_vaccination_coverage,
    )

    updated_dynamics_parameter_specification = recreate_noise_dynamics_spec(
        updated_dynamical_noise_spec,
        ensemble_specification,
    )

    # Calculate endemic equilibrium proportions if dynamics_parameter_specification is provided
    endemic_props_result = calculate_endemic_equilibrium_proportions(
        updated_dynamics_parameter_specification,
        mean_vaccination_coverage
    )

    updated_state_parameters = if Try.isok(endemic_props_result)
        endemic_props = Try.unwrap(endemic_props_result)
        StateParameters(
            ; N = N,
            s_prop = endemic_props.s_prop,
            e_prop = endemic_props.e_prop,
            i_prop = endemic_props.i_prop,
        )
    else
        if verbose
            @warn Try.unwrap_err(endemic_props_result) *
                "\nDefaulting to no initial infections and proportion susceptible = 1 - vaccination coverage"
        end

        StateParameters(
            ; N = N,
            s_prop = 1 - mean_vaccination_coverage,
            e_prop = 0.0,
            i_prop = 0.0
        )

    end

    updated_ensemble_specification = EnsembleSpecification(
        updated_state_parameters,
        updated_dynamics_parameter_specification,
        time_parameters,
        nsims,
        dirpath
    )

    noise_result = create_noise_vecs(
        updated_dynamical_noise_spec,
        updated_ensemble_specification,
        updated_dynamics_parameter_specification,
        endemic_props_result,
        enddates_vec,
    )

    return noise_result
end
