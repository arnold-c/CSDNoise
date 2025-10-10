export generate_single_ensemble,
    trim_ensemble_simulations

"""
    trim_ensemble_simulations(ensemble_run, enddates)

Trim ensemble simulation results to specified end dates.

This function takes an existing ensemble of SEIR simulations and trims both
the emergent and null dynamics scenarios to the specified end dates. This is
useful for post-processing ensemble results to focus analysis on specific
time periods without re-running the full simulation, and required before
creating the noise simulations.

# Arguments
- `ensemble_run::EnsembleSEIRRun`: An ensemble simulation result containing
  both emergent and null SEIR runs as StructVectors of individual simulation
  results.
- `enddates::Vector{Int64}`: Vector of end dates (in simulation time units)
  specifying where to trim each simulation in the ensemble. The length should
  match the number of simulations in the ensemble.

# Returns
- `EnsembleSEIRRun`: A new ensemble structure with both emergent and null
  SEIR runs trimmed to the specified end dates. The structure maintains the
  same format as the input but with shortened time series.
"""
function trim_ensemble_simulations(
        ensemble_run::EnsembleSEIRRun,
        enddates::Vector{Int64}
    )

    trimmed_emergent_run = trim_seir_results(
        ensemble_run.emergent_seir_run,
        enddates
    )

    trimmed_null_run = trim_seir_results(
        ensemble_run.null_seir_run,
        enddates
    )

    return EnsembleSEIRRun(
        emergent_seir_run = trimmed_emergent_run,
        null_seir_run = trimmed_null_run
    )

end

"""
    generate_single_ensemble(ensemble_spec; seed = 1234)

Generate a single ensemble of SEIR simulations with both emergent and null
dynamics scenarios.

This function creates an ensemble of stochastic SEIR model runs, simulating
both an emergent dynamics scenario (with specified vaccination coverage) and
a null dynamics scenario (counterfactual with different vaccination coverage).
Each simulation in the ensemble uses the same seasonal transmission pattern
but different random seeds to capture stochastic variation.

# Arguments
- `ensemble_spec::EnsembleSpecification`: Specification containing all
  parameters for the ensemble simulation, including state parameters, time
  parameters, dynamics specifications for both emergent and null scenarios,
  and the number of simulations to run.
- `seed::Int64 = 1234`: Base random seed for reproducibility. Each simulation
  in the ensemble uses `seed + (sim - 1)` to ensure different but reproducible
  random trajectories.

# Returns
- `EnsembleSEIRRun`: A structure containing two `StructVector{SEIRRun}`
  objects, one for the emergent dynamics scenario and one for the null
  dynamics scenario. Each vector contains `nsims` individual SEIR simulation
  results.

# Notes
Both emergent and null simulations use the same seasonal transmission pattern,
initial states, and random seeds. The only difference is in the vaccination
coverage parameters specified in their respective dynamics specifications.
"""
function generate_single_ensemble(
        ensemble_spec::EnsembleSpecification;
        seed::Int64 = 1234
    )
    @unpack state_parameters,
        emergent_dynamics_parameter_specification,
        null_dynamics_parameter_specification,
        time_parameters,
        nsims = ensemble_spec

    @unpack tlength = time_parameters

    init_states_sv = StaticArrays.SVector(state_parameters.init_states)

    # Get concrete type to avoid abstract element types
    emergent_seir_results = Vector{SEIRRun}(undef, nsims)
    null_seir_results = Vector{SEIRRun}(undef, nsims)

    beta_vec = zeros(Float64, tlength)
    calculate_beta_amp!(
        beta_vec,
        emergent_dynamics_parameter_specification,
        time_parameters
    )

    for sim in eachindex(emergent_seir_results)
        run_seed = seed + (sim - 1)

        local emergent_seir_res = _create_run_seir_results(
            emergent_dynamics_parameter_specification,
            init_states_sv,
            beta_vec,
            time_parameters,
            run_seed
        )
        emergent_seir_results[sim] = emergent_seir_res

        null_seir_results[sim] = _create_run_seir_results(
            null_dynamics_parameter_specification,
            init_states_sv,
            beta_vec,
            time_parameters,
            run_seed
        )
    end

    return EnsembleSEIRRun(
        emergent_seir_run = StructVector(emergent_seir_results),
        null_seir_run = StructVector(null_seir_results),
    )
end

"""
    _create_run_seir_results(dynamics_parameter_specification, init_states,
                             beta_vec, time_parameters, run_seed = 1234)

Internal helper function to create a single SEIR simulation run.

This function constructs a `DynamicsParameters` object from the specification
and runs the SEIR model with the provided initial conditions, transmission
pattern, and time parameters.

# Arguments
- `dynamics_parameter_specification::DynamicsParameterSpecification`:
  Specification containing disease dynamics parameters including transmission
  rates, recovery rates, and vaccination coverage ranges.
- `init_states::StaticArrays.SVector{5, Int64}`: Initial state vector with
  five compartments [S, E, I, R, V] representing susceptible, exposed,
  infectious, recovered, and vaccinated populations.
- `beta_vec::Vector{Float64}`: Pre-calculated seasonal transmission rate
  pattern over the simulation time period.
- `time_parameters::SimTimeParameters`: Time-related parameters including
  simulation length, time step, and burnin period.
- `run_seed::Int64 = 1234`: Random seed for this specific simulation run.

# Returns
- `SEIRRun`: A structure containing the simulation results including state
  trajectories, incidence time series, and effective reproduction number
  (Reff) time series.

# Notes
This is an internal function (indicated by the leading underscore) used by
`generate_single_ensemble` to avoid code duplication when creating both
emergent and null scenario simulations.
"""
function _create_run_seir_results(
        dynamics_parameter_specification::DynamicsParameterSpecification,
        init_states::StaticArrays.SVector{5, Int64},
        beta_vec::Vector{Float64},
        time_parameters::SimTimeParameters,
        run_seed = 1234
    )
    dynp = DynamicsParameters(
        dynamics_parameter_specification;
        seed = run_seed
    )

    return seir_mod(
        init_states,
        dynp,
        beta_vec,
        time_parameters;
        seed = run_seed,
    )
end
