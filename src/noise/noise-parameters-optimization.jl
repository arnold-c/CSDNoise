export optimize_dynamic_noise_params_wrapper,
    optimize_dynamic_noise_params


"""
    optimize_dynamic_noise_params_wrapper(
        target_scaling,
        seir_results::StructVector{SEIRRun},
        thresholds::StructVector{T},
        ews_enddate_type::EWSEndDateType,
        dynamical_noise_spec::DynamicalNoiseSpecification,
        ensemble_specification::EnsembleSpecification,
        optimization_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters()
    ) where {T<:AbstractThresholds}

Wrapper function that calculates vaccination coverage using SEIR results and EWS endpoints.

This function prepares data by calculating endpoints and mean incidence, then calls the
original optimization function with the computed mean incidence value.

# Arguments
- `target_scaling`: Multiplicative factor for target noise level
- `seir_results`: StructVector of SEIR simulation results
- `thresholds`: StructVector of threshold objects for endpoint calculation
- `ews_enddate_type`: Type of endpoint to calculate
- `dynamical_noise_spec`: DynamicalNoiseSpecification containing noise model parameters
- `ensemble_specification`: Ensemble simulation parameters
- `optimization_params`: NoiseVaccinationOptimizationParameters containing algorithm settings (optional)

# Returns
- Same as `optimize_dynamic_noise_params`

# Example
```julia
noise_spec = DynamicalNoiseSpecification(
    R0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 1.0
)
result = optimize_dynamic_noise_params(
    7.0,  # target_scaling
    seir_results,
    thresholds,
    EWSEndDateType(ReffStart()),
    noise_spec,
    ensemble_spec
)
```
"""
function optimize_dynamic_noise_params_wrapper(
        target_scaling,
        seir_results::StructVector{SEIRRun},
        thresholds::StructVector{T},
        ews_enddate_type::EWSEndDateType,
        dynamical_noise_spec::DynamicalNoiseSpecification,
        ensemble_specification::EnsembleSpecification,
        optimization_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters(),
        verbose = false
    ) where {T <: AbstractThresholds}

    # Calculate endpoints for all simulations
    enddates_vec = calculate_all_ews_enddates(thresholds, ews_enddate_type)

    filtered_seir_results = trim_seir_results(seir_results, enddates_vec)

    # Calculate filtered mean incidence up to endpoints
    overall_mean = calculate_mean_incidence(filtered_seir_results)

    # Call the original optimization function with the computed mean
    return optimize_dynamic_noise_params(
        target_scaling,
        overall_mean,
        dynamical_noise_spec,
        ensemble_specification,
        enddates_vec,
        optimization_params;
        verbose = verbose
    )
end

"""
    optimize_dynamic_noise_params(
        target_scaling,
        measles_daily_incidence,
        dynamical_noise_spec::DynamicalNoiseSpecification,
        ensemble_specification::EnsembleSpecification,
        enddates_vec,
        optimization_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters()
    )

Optimize vaccination coverage and initial susceptible proportion using multistart optimization to achieve target noise scaling.

Uses Sobol sequences for initial points and local optimization to find the vaccination level
and initial susceptible proportion that produces mean noise equal to `target_scaling * measles_daily_incidence`.

# Arguments
- `target_scaling`: Multiplicative factor for target noise level
- `measles_daily_incidence`: Daily measles incidence to scale
- `dynamical_noise_spec`: DynamicalNoiseSpecification containing noise model parameters and search bounds
- `ensemble_specification`: Ensemble simulation parameters
- `enddates_vec`: A vector of all simulation enddates
- `optimization_params`: NoiseVaccinationOptimizationParameters containing algorithm settings (optional, uses defaults if not provided)
- `verbose`: Describe steps being conducted

# Returns
- Named tuple with `optimal_vaccination`, `optimal_susceptible_proportion`, `mean_noise`, `target_noise`, and `difference`

# Example
```julia
noise_spec = DynamicalNoiseSpecification(
    R0 = 5.0,
    latent_period = 7,
    duration_infection = 14,
    correlation = "in-phase",
    poisson_component = 1.0
)
result = calculate_dynamic_vaccination_coverage(
    7.0,  # target_scaling
    3.0,  # measles_daily_incidence (target noise = 21.0)
    noise_spec,
    ensemble_spec,
    enddates_vec
)
```
"""
function optimize_dynamic_noise_params(
        target_scaling,
        measles_daily_incidence,
        dynamical_noise_spec::DynamicalNoiseSpecification,
        ensemble_specification::EnsembleSpecification,
        enddates_vec,
        optimization_params::NoiseVaccinationOptimizationParameters = NoiseVaccinationOptimizationParameters();
        verbose = false,
    )

    @unpack vaccination_bounds,
        susceptible_bounds = dynamical_noise_spec

    @unpack n_sobol_points, local_algorithm,
        xtol_rel,
        xtol_abs,
        maxeval,
        atol = optimization_params

    target_noise = target_scaling * measles_daily_incidence

    if verbose
        println("Target noise level: $target_noise")
        println("Vaccination bounds: $vaccination_bounds")
        println("Susceptible proportion bounds: $susceptible_bounds")
        println("Starting multistart optimization with $n_sobol_points points...")
    end

    # Pre-extract ensemble specification components to avoid repeated unpacking
    N = ensemble_specification.state_parameters.init_states.N

    # Define objective function: minimize squared difference from target
    objective = let target_noise = target_noise,
            dynamical_noise_spec = dynamical_noise_spec,
            ensemble_specification = ensemble_specification,
            enddates_vec = enddates_vec,
            N = N
        function (params)
            vaccination_coverage = params[1]
            susceptible_proportion = 1 - vaccination_coverage

            # Ensure susceptible proportion is valid and will result in positive compartments
            # Use stricter bounds to prevent numerical issues with very small populations
            min_safe_prop = max(1.0 / N, 0.001)  # At least 1 person or 0.1%, whichever is larger
            max_safe_prop = min(1.0 - 1.0 / N, 0.999)  # At most N-1 people or 99.9%, whichever is smaller

            if susceptible_proportion <= min_safe_prop || susceptible_proportion >= max_safe_prop
                return 1.0e10  # Large penalty for unsafe proportions
            end

            noise_level = calculate_mean_dynamical_noise(
                dynamical_noise_spec,
                vaccination_coverage,
                susceptible_proportion,
                ensemble_specification,
                enddates_vec
            )

            # Return squared error from target
            return (noise_level - target_noise)^2
        end
    end

    # Setup multistart optimization problem
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [vaccination_bounds[1]], # lower bounds
        [vaccination_bounds[2]]  # upper bounds
    )

    # Configure local optimization method
    local_method = MultistartOptimization.NLopt_local_method(
        local_algorithm;
        xtol_rel = xtol_rel,
        xtol_abs = xtol_abs,
        maxeval = maxeval,
    )

    # Configure multistart method using Sobol sequences
    multistart_method = MultistartOptimization.TikTak(n_sobol_points)

    if verbose
        println("Running multistart optimization...")
    end

    # Run optimization with threading disabled to avoid potential race conditions
    result = MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    if !in(result.ret, [:SUCCESS, :XTOL_REACHED, :FTOL_REACHED, :STOPVAL_REACHED]) ||
            sqrt(result.value) > atol

        optimal_vaccination_coverage = result.location[1]
        optimal_susceptible_proportion = 1 - optimal_vaccination_coverage

        optimized_noise = calculate_mean_dynamical_noise(
            dynamical_noise_spec,
            optimal_vaccination_coverage,
            optimal_susceptible_proportion,
            ensemble_specification,
            enddates_vec
        )
        error_msg = "Unsuccessful optimization." *
            "\nReturn code: $(result.ret)" *
            "\nAbsolute difference: $(sqrt(result.value))" *
            "\nTarget value: $(target_noise)" *
            "\nOptimized value: $(optimized_noise)" *
            "\nOptimization parameters: Vax. prop = $(optimal_vaccination_coverage), Sus. prop = $(optimal_susceptible_proportion)"


        error(error_msg)
    end

    return result
end
