using UnPack: @unpack
using FixedSizeArrays: FixedSizeVector
using Try: Try

"""
    calculate_dynamic_vaccination_coverage_multistart(
        target_scaling,
        measles_daily_incidence,
        dynamical_noise_specification_parameters,
        ensemble_specification;
        vaccination_bounds = [0.0, 1.0],
        max_vaccination_range = 0.2,
        n_sobol_points = 50,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
        verbose = false
    )

Optimize vaccination coverage using multistart optimization to achieve target noise scaling.

Uses Sobol sequences for initial points and local optimization to find the vaccination level
that produces mean noise equal to `target_scaling * measles_daily_incidence`.

# Arguments
- `target_scaling`: Multiplicative factor for target noise level
- `measles_daily_incidence`: Daily measles incidence to scale
- `dynamical_noise_specification_parameters`: Parameters for dynamical noise model
- `ensemble_specification`: Ensemble simulation parameters
- `vaccination_bounds`: [min, max] bounds for vaccination coverage search
- `max_vaccination_range`: Maximum range around mean vaccination coverage
- `n_sobol_points`: Number of Sobol sequence starting points for multistart
- `local_algorithm`: NLopt algorithm for local optimization
- `maxeval`: Maximum function evaluations per local optimization
- `xtol_rel`, `xtol_abs`, `ftol_rel`: Convergence tolerances
- `verbose`: Print optimization progress

# Returns
- `(optimal_vaccination, achieved_noise)`: Tuple of optimal vaccination coverage and resulting noise level

# Example
```julia
optimal_vaccination, achieved_noise = calculate_dynamic_vaccination_coverage_multistart(
    7.0,  # target_scaling
    3.0,  # measles_daily_incidence (target noise = 21.0)
    dynamical_noise_params,
    ensemble_spec
)
```
"""
function calculate_dynamic_vaccination_coverage_multistart(
        target_scaling,
        measles_daily_incidence,
        dynamical_noise_specification_parameters,
        ensemble_specification;
        vaccination_bounds = [0.0, 1.0],
        max_vaccination_range = 0.2,
        n_sobol_points = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval = 1000,
        xtol_rel = 1.0e-3,
        xtol_abs = 1.0e-3,
        ftol_rel = 1.0e-4,
        verbose = false
    )
    @unpack R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component = dynamical_noise_specification_parameters

    @assert length(vaccination_bounds) == 2
    @assert vaccination_bounds[1] < vaccination_bounds[2]

    target_noise = target_scaling * measles_daily_incidence

    if verbose
        println("Target noise level: $target_noise")
        println("Vaccination bounds: $vaccination_bounds")
        println("Starting multistart optimization with $n_sobol_points points...")
    end

    # Define objective function: minimize squared difference from target
    function objective(params)
        vaccination_coverage = params[1]

        # Ensure vaccination coverage is within bounds
        noise_level = calculate_mean_dynamical_noise(
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            vaccination_coverage,
            max_vaccination_range,
            ensemble_specification,
        )

        # Return squared error from target
        return (noise_level - target_noise)^2
    end

    # Setup multistart optimization problem
    problem = MultistartOptimization.MinimizationProblem(
        objective,
        [vaccination_bounds[1]],  # lower bounds
        [vaccination_bounds[2]]   # upper bounds
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

    # Run optimization
    result = MultistartOptimization.multistart_minimization(
        multistart_method,
        local_method,
        problem
    )

    optimal_vaccination = result.location[1]
    optimal_objective = result.value

    # Calculate the actual noise level achieved
    achieved_noise = calculate_mean_dynamical_noise(
        R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        optimal_vaccination,
        max_vaccination_range,
        ensemble_specification,
    )

    if verbose
        println("Optimization complete!")
        println("Optimal vaccination coverage: $optimal_vaccination")
        println("Achieved noise level: $achieved_noise")
        println("Target noise level: $target_noise")
        println("Final objective value (squared error): $optimal_objective")
        println("Absolute error: $(abs(achieved_noise - target_noise))")
    end

    return (;
        optimal_vaccination = round(optimal_vaccination; digits = 4),
        mean_noise = round(achieved_noise; digits = 4),
        target_noise = round(target_noise; digits = 4),
        difference = round(achieved_noise - target_noise; digits = 4),
    )
end

"""
    calculate_mean_dynamical_noise(
        R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        mean_vaccination_coverage,
        max_vaccination_range,
        ensemble_specification
    )

Calculate the mean dynamical noise level for given vaccination coverage parameters.

Creates a dynamical noise specification with vaccination coverage bounds and runs
ensemble simulations to compute the resulting mean noise level.

# Arguments
- `R0`: Basic reproduction number of the noise disease
- `latent_period`: Duration of latent period of the noise in the SEIR model
- `duration_infection`: Duration of infectious period of the noise
- `correlation`: Temporal correlation parameter for noise relative to the target disease ("in-phase", "out-of-phase")
- `poisson_component`: The amount of additional Poisson noise to add (relative to the mean dynamical noise generated in a simulation)
- `mean_vaccination_coverage`: Mean vaccination coverage level (0-1)
- `max_vaccination_range`: Maximum range around mean coverage for bounds to (uniformly) sample within
- `ensemble_specification`: Parameters for ensemble simulation

# Returns
- `Float64`: Mean noise level from the ensemble simulation

# Example
```julia
mean_noise = calculate_mean_dynamical_noise(
    12.0,     # R0
    8.0,      # latent_period
    7.0,      # duration_infection
    0.9,      # correlation
    1.0,      # poisson_component
    0.85,     # mean_vaccination_coverage
    0.2,      # max_vaccination_range
    ensemble_spec
)
```
"""
function calculate_mean_dynamical_noise(
        R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        mean_vaccination_coverage,
        max_vaccination_range,
        ensemble_specification,
    )
    min_vaccination_coverage,
        max_vaccination_coverage = calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range,
    )

    dynamical_noise_spec = NoiseSpecification(
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

    noise_result = create_noise_arr(
        dynamical_noise_spec,
        ensemble_specification
    )

    return noise_result.mean_noise
end


# end
