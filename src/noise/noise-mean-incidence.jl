export calculate_mean_dynamical_noise

"""
    calculate_mean_dynamical_noise(
        dynamical_noise_spec::DynamicalNoiseSpecification,
        mean_vaccination_coverage,
        susceptible_proportion,
        ensemble_specification::EnsembleSpecification,
        enddates_vec
    )

Calculate the mean dynamical noise level for given vaccination coverage parameters.

Creates a dynamical noise specification with vaccination coverage bounds and runs
ensemble simulations to compute the resulting mean noise level.

# Arguments
- `dynamical_noise_spec`: DynamicalNoiseSpecification containing fixed noise parameters
- `mean_vaccination_coverage`: Mean vaccination coverage level (0-1) - optimization variable
- `susceptible_proportion`: Initial susceptible proportion (0-1) - optimization variable
- `ensemble_specification`: Parameters for ensemble simulation
- `enddates_vec`: Vector of simulation end dates

# Returns
- `Float64`: Mean noise level from the ensemble simulation

# Example
```julia
noise_spec = DynamicalNoiseSpecification(
    R0 = 12.0,
    latent_period = 8,
    duration_infection = 7,
    correlation = "in-phase",
    poisson_component = 1.0
)
mean_noise = calculate_mean_dynamical_noise(
    noise_spec,
    0.85,     # mean_vaccination_coverage
    0.7,      # susceptible_proportion
    ensemble_spec,
    enddates_vec
)
```
"""
function calculate_mean_dynamical_noise(
        dynamical_noise_spec::DynamicalNoiseSpecification,
        mean_vaccination_coverage,
        ensemble_specification::EnsembleSpecification,
        enddates_vec;
        verbose = false
    )

    noise_result = recreate_noise_vecs(
        dynamical_noise_spec,
        mean_vaccination_coverage,
        ensemble_specification,
        enddates_vec;
        verbose = verbose
    )
    )[1]

    return noise_result.mean_noise
end
