using UnPack: @unpack

export recreate_noise_vecs

"""
    recreate_noise_vecs(
        dynamical_noise_spec::DynamicalNoiseSpecification,
        mean_vaccination_coverage,
        susceptible_proportion,
        ensemble_specification::EnsembleSpecification,
        enddates_vec
    )

Recreate noise vectors with updated vaccination coverage and susceptible proportion parameters.

Takes updated susceptible proportion and vaccination coverage values (either during or after
optimization), recreates the state parameters, ensemble specification, and dynamics
specifications, then generates new noise vectors using the updated parameters.

# Arguments
- `dynamical_noise_spec`: DynamicalNoiseSpecification containing fixed noise parameters
- `mean_vaccination_coverage`: Mean vaccination coverage level (0-1) - optimization variable
- `susceptible_proportion`: Initial proportion of population that is susceptible (0-1) - optimization variable
- `ensemble_specification`: Parameters for ensemble simulation
- `enddates_vec`: Vector of simulation end dates

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
    enddates_vec
)
```
"""
function recreate_noise_vecs(
        dynamical_noise_spec::DynamicalNoiseSpecification,
        mean_vaccination_coverage,
        susceptible_proportion,
        ensemble_specification::EnsembleSpecification,
        enddates_vec
    )
    @unpack R0, latent_period, duration_infection, correlation, poisson_component, max_vaccination_range = dynamical_noise_spec

    # Create final EnsembleSpecification with optimal parameters for verification
    @unpack state_parameters, dynamics_parameter_specification, time_parameters, nsims, dirpath = ensemble_specification
    @unpack init_states = state_parameters
    @unpack N = init_states

    final_state_parameters = StateParameters(;
        N = N,
        s_prop = susceptible_proportion,
        e_prop = 0.0,
        i_prop = 0.0,
    )

    final_ensemble_specification = EnsembleSpecification(
        final_state_parameters,
        dynamics_parameter_specification,
        time_parameters,
        nsims,
        dirpath
    )

    min_vaccination_coverage, max_vaccination_coverage = calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range,
    )

    dynamical_noise_spec_sumtype = NoiseSpecification(
        DynamicalNoise(
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            min_vaccination_coverage,
            max_vaccination_coverage,
        )
    )

    noise_result = create_noise_vecs(
        dynamical_noise_spec_sumtype,
        final_ensemble_specification,
        enddates_vec,
    )

    return noise_result
end
