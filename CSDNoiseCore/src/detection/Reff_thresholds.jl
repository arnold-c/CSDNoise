export calculate_all_Reff_thresholds,
    Reff_ge_than_one

"""
    calculate_all_Reff_thresholds(ensemble_run::EnsembleSEIRRun)

Calculate Reff >= 1 threshold periods for all simulations in an ensemble run.

This function processes an ensemble of SEIR simulation results and identifies all periods
where the effective reproduction number (Reff) is greater than or equal to 1 for each
simulation in the ensemble. It applies the threshold detection to the emergent (non-null)
simulations only.

The function iterates through each simulation in the ensemble and applies the
[`Reff_ge_than_one`](@ref) function to extract threshold crossing periods, returning
the results in a structured format suitable for further analysis.

# Arguments
- `ensemble_run::EnsembleSEIRRun`: Ensemble simulation results containing:
  - `emergent_seir_run`: StructVector of SEIRRun results for emergent scenarios
  - `null_seir_run`: StructVector of SEIRRun results for null scenarios (not used)

# Returns
- `StructVector{Thresholds}`: A structured vector where each element contains threshold
  information for one simulation, with fields:
  - `lower_bounds`: Starting indices of Reff >= 1 periods
  - `upper_bounds`: Ending indices of Reff >= 1 periods
  - `duration`: Duration of each Reff >= 1 period

# Examples
```julia
# Process ensemble results to find Reff >= 1 periods
ensemble_results = simulate_ensemble_seir_results(...)
reff_thresholds = calculate_all_Reff_thresholds(ensemble_results)

# Access threshold data for first simulation
first_sim_thresholds = reff_thresholds[1]
println("Number of outbreak periods: ", length(first_sim_thresholds.lower_bounds))
```

# See Also
- [`Reff_ge_than_one`](@ref): Function applied to individual Reff timeseries
- [`EnsembleSEIRRun`](@ref): Input type containing ensemble simulation results
- [`Thresholds`](@ref): Return type for threshold period information
"""
function calculate_all_Reff_thresholds(ensemble_run::EnsembleSEIRRun)
    nsims = length(ensemble_run.emergent_seir_run)

    Reff_vecs = Vector{Thresholds}(undef, nsims)
    for i in eachindex(Reff_vecs)
        Reff_vecs[i] = Reff_ge_than_one(ensemble_run.emergent_seir_run[i].Reff)
    end

    return StructVector(Reff_vecs)
end

"""
    Reff_ge_than_one(Reff_vec)

Identify periods where the effective reproduction number (Reff) is greater than or equal to 1.

This function takes a timeseries of Reff values and identifies all consecutive periods
where Reff >= 1, which typically indicates potential outbreak or epidemic growth phases.
The function uses run-length encoding to efficiently detect threshold crossings and
calculate the bounds and durations of these critical periods.

The detection process involves:
1. Converting the Reff timeseries to a boolean series (Reff >= 1)
2. Applying run-length encoding to identify consecutive periods
3. Extracting bounds and durations using [`calculate_above_threshold_bounds`](@ref)

# Arguments
- `Reff_vec`: Vector of effective reproduction numbers (typically Float64 values)

# Returns
- `Thresholds`: A struct containing:
  - `lower_bounds`: Vector of starting indices for each Reff >= 1 period
  - `upper_bounds`: Vector of ending indices for each Reff >= 1 period
  - `duration`: Vector of durations (in time steps) for each Reff >= 1 period

# Examples
```julia
# Detect outbreak periods in an Reff timeseries
Reff_timeseries = [0.8, 0.9, 1.1, 1.2, 1.3, 0.9, 0.8, 1.1, 1.0, 0.9]
outbreak_periods = Reff_ge_than_one(Reff_timeseries)

# Check results
println("Number of outbreak periods: ", length(outbreak_periods.lower_bounds))
println("First outbreak: indices ", outbreak_periods.lower_bounds[1],
        " to ", outbreak_periods.upper_bounds[1])
println("Duration: ", outbreak_periods.duration[1], " time steps")
```

# See Also
- [`calculate_above_threshold_bounds`](@ref): Core function for threshold period detection
- [`calculate_all_Reff_thresholds`](@ref): Function that applies this to ensemble results
- [`Thresholds`](@ref): Return type containing threshold period information
- [`StatsBase.rle`](@ref): Run-length encoding function used internally
"""
function Reff_ge_than_one(Reff_vec)
    Reff_rle = StatsBase.rle(Reff_vec .>= 1)
    return calculate_above_threshold_bounds(Reff_rle)
end
