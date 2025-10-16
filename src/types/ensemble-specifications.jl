export EnsembleSpecification,
    EnsembleSpecsParameters

"""
    EnsembleSpecification

Complete specification for ensemble disease simulation parameters and configuration.

This struct defines all parameters needed to run an ensemble of SEIR disease simulations,
including initial population states, time parameters, disease dynamics for both emergent
and null scenarios, noise characteristics, and output configuration. It serves as the
primary configuration object for ensemble-based epidemiological modeling and analysis.

The `dirpath` field is automatically constructed from the component specifications to
create a unique filesystem path for storing ensemble simulation results, incorporating
all relevant parameter values in a hierarchical directory structure.

# Fields
- `label::String`: The disease name
- `state_parameters::StateParameters`: Initial population state configuration (S, E, I, R compartments)
- `time_parameters::SimTimeParameters`: Simulation time configuration including burnin, duration, and timestep
- `emergent_dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for outbreak scenarios
- `null_dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for null/baseline scenarios
- `dynamical_noise_specification::DynamicalNoiseSpecification`: Specification for dynamical noise in disease parameters
- `nsims::Int64`: Number of simulation replicates in the ensemble
- `dirpath::String`: Automatically generated directory path for storing ensemble results

# Constructor
    EnsembleSpecification(; label, state_parameters, time_parameters, emergent_dynamics_parameter_specification,
                         null_dynamics_parameter_specification, dynamical_noise_specification, nsims, dirpath)

# Example
```julia
# Create ensemble specification with all components
ensemble_spec = EnsembleSpecification(
    state_parameters = StateParameters(N=500_000, s_prop=0.9, e_prop=0.05, i_prop=0.05),
    time_parameters = SimTimeParameters(burnin=Year(10), tmax=20.0, tstep=0.1),
    emergent_dynamics_parameter_specification = emergent_dynamics,
    null_dynamics_parameter_specification = null_dynamics,
    dynamical_noise_specification = noise_spec,
    nsims = 1000,
    dirpath = "/path/to/results"
)

# Access ensemble configuration
println("Running \$(ensemble_spec.nsims) simulations")
println("Results stored in: \$(ensemble_spec.dirpath)")
```

# See Also
- [`StateParameters`](@ref): Initial population state configuration
- [`SimTimeParameters`](@ref): Simulation time parameters
- [`DynamicsParameterSpecification`](@ref): Disease dynamics specifications
- [`DynamicalNoiseSpecification`](@ref): Dynamical noise parameters
- [`simulate_ensemble_seir_results`](@ref): Function that uses EnsembleSpecification for simulation
"""
AutoHashEquals.@auto_hash_equals Base.@kwdef struct EnsembleSpecification
    label::String
    state_parameters::StateParameters
    time_parameters::SimTimeParameters
    emergent_dynamics_parameter_specification::DynamicsParameterSpecification
    null_dynamics_parameter_specification::DynamicsParameterSpecification
    dynamical_noise_specification::DynamicalNoiseSpecification
    nsims::Int64
    dirpath::String
end

"""
    EnsembleSpecification(label, state_parameters, time_parameters, emergent_dynamics_parameter_specification,
                         null_dynamics_parameter_specification, dynamical_noise_specification, nsims)

Construct an EnsembleSpecification with automatic directory path generation.

This constructor creates a complete ensemble specification by combining all the individual
component specifications and automatically generating a hierarchical directory path that
uniquely identifies this ensemble configuration. The directory path incorporates all
relevant parameter values to ensure unique storage locations for different ensemble runs.

The generated directory structure follows this hierarchy:
- Base: "ensemble/seasonal-infectivity-import/tau-leaping"
- Population parameters: N, initial recovered proportion
- Simulation parameters: number of simulations, time parameters
- Disease dynamics: R₀, latent period, infectious period for both emergent and null scenarios
- Noise parameters: R₀ noise, latent period noise, infectious period noise
- Vaccination parameters: burnin and post-burnin coverage ranges
- Birth rate and transmission parameters

# Arguments
- `label::String`: The disease name
- `state_parameters::StateParameters`: Initial population state configuration
- `time_parameters::SimTimeParameters`: Simulation time configuration
- `emergent_dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for outbreak scenarios
- `null_dynamics_parameter_specification::DynamicsParameterSpecification`: Disease dynamics for baseline scenarios
- `dynamical_noise_specification::DynamicalNoiseSpecification`: Dynamical noise specification
- `nsims::Int64`: Number of simulation replicates

# Returns
- `EnsembleSpecification`: Complete ensemble specification with auto-generated directory path

# Example
```julia
# Create ensemble specification with automatic path generation
ensemble_spec = EnsembleSpecification(
	label,
    state_params,
    time_params,
    emergent_dynamics,
    null_dynamics,
    noise_spec,
    1000  # number of simulations
)

# The dirpath is automatically generated based on all parameters
println(ensemble_spec.dirpath)
# Output: "ensemble/seasonal-infectivity-import/tau-leaping/N_500000/r_0.95/nsims_1000/..."
```
"""
function EnsembleSpecification(
        label::String,
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
        "$label",
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
        "burnin_$(Float64(Dates.days(time_parameters.burnin)))",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)",
    )

    return EnsembleSpecification(
        label,
        state_parameters,
        time_parameters,
        emergent_dynamics_parameter_specification,
        null_dynamics_parameter_specification,
        dynamical_noise_specification,
        nsims,
        dirpath,
    )
end

"""
    EnsembleSpecsParameters

Legacy ensemble specification parameters for disease simulation configuration.

This struct provides an alternative parameterization for ensemble simulations using
a more direct specification of disease and vaccination parameters. It includes built-in
validation to ensure parameter consistency and logical constraints.

**Note**: This struct appears to be legacy code (marked with TODO) and may not be
actively used in the current codebase. Consider using [`EnsembleSpecification`](@ref)
for new implementations.

# Fields
- `burnin_years::Int64`: Number of years for simulation burnin period
- `tmax_years::Int64`: Total simulation duration in years
- `annual_births_per_k::Int64`: Annual births per 1000 population
- `ensemble_state_specification::StateParameters`: Initial population state configuration
- `R_0::Float64`: Basic reproduction number
- `gamma::Float64`: Recovery rate (1/infectious period)
- `sigma::Float64`: Progression rate from exposed to infectious (1/latent period)
- `target_Reff::Float64`: Target effective reproduction number for vaccination scenarios
- `target_years::Int64`: Years over which to achieve target Reff
- `min_vaccination_coverage::Float64`: Minimum vaccination coverage proportion
- `max_vaccination_coverage::Float64`: Maximum vaccination coverage proportion
- `nsims::Int64`: Number of simulation replicates

# Validation
The constructor enforces these constraints:
- `tmax_years >= target_years`: Total simulation time must accommodate target period
- `min_vaccination_coverage < max_vaccination_coverage`: Valid coverage range

# Constructors
    EnsembleSpecsParameters(burnin_years, tmax_years, annual_births_per_k, ensemble_state_specification,
                           R_0, gamma, sigma, target_Reff, target_years, min_vaccination_coverage,
                           max_vaccination_coverage, nsims)

    EnsembleSpecsParameters(; burnin_years, tmax_years, ...)

# Example
```julia
# Create ensemble parameters with validation
ensemble_params = EnsembleSpecsParameters(
    burnin_years = 50,
    tmax_years = 100,
    target_Reff = 0.9,
    target_years = 75,  # Must be <= tmax_years
    min_vaccination_coverage = 0.6,
    max_vaccination_coverage = 0.8,  # Must be > min_vaccination_coverage
    nsims = 1000
)
```

# See Also
- [`EnsembleSpecification`](@ref): Modern ensemble specification interface
- [`StateParameters`](@ref): Population state configuration
- [`DynamicsParameterSpecification`](@ref): Disease dynamics parameters
"""
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

"""
    EnsembleSpecsParameters(; burnin_years, tmax_years, annual_births_per_k=ANNUAL_BIRTHS_PER_K,
                           ensemble_state_specification=StateParameters(...), R_0=R0, gamma=GAMMA,
                           sigma=SIGMA, target_Reff=0.9, target_years=2*burnin_years,
                           min_vaccination_coverage=0.6, max_vaccination_coverage=0.8, nsims=1000)

Construct EnsembleSpecsParameters using keyword arguments with sensible defaults.

This constructor provides a convenient interface for creating ensemble parameters with
default values for most parameters. It automatically sets reasonable defaults based on
common epidemiological scenarios and project constants.

# Keyword Arguments
- `burnin_years::Int`: Number of years for simulation burnin (required)
- `tmax_years::Int`: Total simulation duration in years (required)
- `annual_births_per_k::Int64`: Annual births per 1000 population (default: `ANNUAL_BIRTHS_PER_K`)
- `ensemble_state_specification::StateParameters`: Initial state (default: 500k population, 5% susceptible, 95% recovered)
- `R_0::Float64`: Basic reproduction number (default: `R0` constant)
- `gamma::Float64`: Recovery rate (default: `GAMMA` constant)
- `sigma::Float64`: Progression rate (default: `SIGMA` constant)
- `target_Reff::Float64`: Target effective reproduction number (default: 0.9)
- `target_years::Int`: Years to achieve target Reff (default: 2 × burnin_years)
- `min_vaccination_coverage::Float64`: Minimum vaccination coverage (default: 0.6)
- `max_vaccination_coverage::Float64`: Maximum vaccination coverage (default: 0.8)
- `nsims::Int`: Number of simulation replicates (default: 1000)

# Returns
- `EnsembleSpecsParameters`: Configured ensemble parameters with validation applied

# Example
```julia
# Minimal specification with defaults
params = EnsembleSpecsParameters(
    burnin_years = 50,
    tmax_years = 100
)

# Custom specification overriding defaults
params = EnsembleSpecsParameters(
    burnin_years = 30,
    tmax_years = 80,
    target_Reff = 0.85,
    min_vaccination_coverage = 0.7,
    max_vaccination_coverage = 0.9,
    nsims = 2000
)
```
"""
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
