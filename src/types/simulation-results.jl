export SEIRRun,
    EnsembleSEIRRun,
    NoiseRun

"""
    SEIRRun

Results from a single SEIR epidemic simulation run.

This struct stores the complete time series output from a single realization of the SEIR 
epidemic model, including population states, incidence, and effective reproduction number.
The data represents the trajectory of an epidemic over time and is used for analysis of
epidemic dynamics and early warning signal detection.

# Fields
- `states::Vector{StaticArrays.SVector{5, Int64}}`: Time series of population states [S, E, I, R, N] for each time step
- `incidence::Vector{Int64}`: Time series of new infections (incidence) at each time step
- `Reff::Vector{Float64}`: Time series of effective reproduction number at each time step

# Example
```julia
# Create a SEIRRun from simulation results
seir_run = SEIRRun(
    states = [SVector(9000, 100, 50, 850, 10000), SVector(8950, 120, 60, 870, 10000)],
    incidence = [50, 60],
    Reff = [1.2, 1.1]
)

# Access time series data
seir_run.states[1].S  # Susceptible population at time 1
seir_run.incidence[end]  # Final incidence value
```

# See Also
- [`EnsembleSEIRRun`](@ref): Container for multiple SEIR runs
- [`simulate_ensemble_seir_results`](@ref): Function that generates SEIRRun objects
"""
Base.@kwdef struct SEIRRun
    states::Vector{StaticArrays.SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

"""
    EnsembleSEIRRun

Container for ensemble simulation results comparing emergent and null scenarios.

This struct holds the results from ensemble simulations that compare two scenarios:
emergent scenarios (where an outbreak occurs) and null scenarios (baseline/control).
Each scenario contains multiple SEIR simulation runs stored efficiently using StructVector
for better memory layout and performance when processing large ensembles.

# Fields
- `emergent_seir_run::StructVector{SEIRRun}`: Collection of SEIR runs for emergent outbreak scenarios
- `null_seir_run::StructVector{SEIRRun}`: Collection of SEIR runs for null/baseline scenarios

# Example
```julia
# Create ensemble runs
emergent_runs = StructVector([seir_run1, seir_run2, seir_run3])
null_runs = StructVector([null_run1, null_run2, null_run3])

ensemble = EnsembleSEIRRun(
    emergent_seir_run = emergent_runs,
    null_seir_run = null_runs
)

# Access ensemble data
ensemble.emergent_seir_run.incidence  # All incidence time series for emergent runs
ensemble.null_seir_run[1]  # First null scenario run
```

# See Also
- [`SEIRRun`](@ref): Individual simulation run structure
- [`simulate_ensemble_seir_results`](@ref): Function that generates ensemble results
"""
Base.@kwdef struct EnsembleSEIRRun
    emergent_seir_run::StructVector{SEIRRun}
    null_seir_run::StructVector{SEIRRun}
end

"""
    AbstractNoiseRun

Abstract base type for noise simulation results.

This abstract type serves as the parent for all noise run result types, providing
a common interface for different noise modeling approaches in epidemic simulations.

# See Also
- [`DynamicalNoiseRun`](@ref): Concrete implementation for dynamical noise
- [`PoissonNoiseRun`](@ref): Concrete implementation for Poisson noise
- [`NoiseRun`](@ref): Sum type encompassing all noise run types
"""
abstract type AbstractNoiseRun end

"""
    DynamicalNoiseRun

Results from noise simulations incorporating both Poisson and dynamical noise components.

This struct stores the results of noise simulations that model both observational noise
(Poisson) and dynamical noise (process noise in the epidemic dynamics). The dynamical
noise captures uncertainty in the underlying epidemic process, while Poisson noise
represents additional uncertainty that could be the result of smaller overlapping
incidence from other noise diseases.

# Fields
- `incidence::Vector{Vector{Int64}}`: Collection of incidence time series with noise applied
- `Reff::Vector{Vector{Float64}}`: Collection of effective reproduction number time series with noise
- `mean_noise::Float64`: Overall mean noise level across all components
- `mean_poisson_noise::Float64`: Mean level of Poisson (observational) noise component
- `mean_dynamic_noise::Float64`: Mean level of dynamical (process) noise component

# Example
```julia
# Create dynamical noise run results
noise_run = DynamicalNoiseRun(
    incidence = [[45, 52, 48], [50, 55, 49]],  # Multiple noisy incidence series
    Reff = [[1.15, 1.08, 1.12], [1.18, 1.09, 1.11]],  # Multiple noisy Reff series
    mean_noise = 0.15,
    mean_poisson_noise = 0.08,
    mean_dynamic_noise = 0.07
)

# Access noise components
noise_run.incidence[1]  # First noisy incidence time series
noise_run.mean_dynamic_noise  # Dynamical noise level
```

# See Also
- [`PoissonNoiseRun`](@ref): Simpler noise model with only Poisson noise
- [`NoiseRun`](@ref): Sum type that can hold either noise run type
"""
Base.@kwdef struct DynamicalNoiseRun
    incidence::Vector{Vector{Int64}}
    Reff::Vector{Vector{Float64}}
    mean_noise::Float64
    mean_poisson_noise::Float64
    mean_dynamic_noise::Float64
end

"""
    PoissonNoiseRun

Results from noise simulations using only Poisson noise.

This struct stores results from simpler noise simulations that only incorporate
Poisson noise to model static uncertainty in incidence data.

# Fields
- `incidence::Vector{Vector{Int64}}`: Collection of incidence time series with Poisson noise applied
- `mean_noise::Float64`: Mean level of Poisson noise applied to the data

# Example
```julia
# Create Poisson noise run results
poisson_run = PoissonNoiseRun(
    incidence = [[48, 52, 49], [51, 47, 53]],  # Multiple noisy incidence series
    mean_noise = 0.12
)

# Access data
poisson_run.incidence[2]  # Second noisy incidence time series
poisson_run.mean_noise    # Noise level
```

# See Also
- [`DynamicalNoiseRun`](@ref): More complex noise model including dynamical noise
- [`NoiseRun`](@ref): Sum type that can hold either noise run type
"""
Base.@kwdef struct PoissonNoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
end

"""
    NoiseRun

Sum type for noise simulation results that can represent either dynamical or Poisson noise runs.

This sum type provides a unified interface for handling different types of noise simulation
results. It can contain either a `DynamicalNoiseRun` (with both dynamical and Poisson noise)
or a `PoissonNoiseRun` (with only Poisson noise), allowing for flexible noise modeling
approaches within the same analysis framework.

# Variants
- `DynamicalNoiseRun`: Noise simulation with both dynamical and observational noise
- `PoissonNoiseRun`: Noise simulation with only observational (Poisson) noise

# Example
```julia
# Create different noise runs
dynamical_run = DynamicalNoiseRun(...)
poisson_run = PoissonNoiseRun(...)

# Both can be stored as NoiseRun
noise1 = NoiseRun(dynamical_run)
noise2 = NoiseRun(poisson_run)

# Pattern matching to handle different types
function process_noise(noise::NoiseRun)
    @match noise begin
        DynamicalNoiseRun(args...) => process_dynamical_noise(args...)
        PoissonNoiseRun(args...) => process_poisson_noise(args...)
    end
end
```

# See Also
- [`DynamicalNoiseRun`](@ref): Dynamical noise simulation results
- [`PoissonNoiseRun`](@ref): Poisson noise simulation results
- [`AbstractNoiseRun`](@ref): Abstract base type for noise runs
"""
LightSumTypes.@sumtype NoiseRun(DynamicalNoiseRun, PoissonNoiseRun) <: AbstractNoiseRun
