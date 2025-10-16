export calculate_outbreak_thresholds,
    calculate_outbreak_thresholds!

"""
    calculate_outbreak_thresholds(
        seir_results::StructVector{SEIRRun},
        outbreak_specification::OutbreakSpecification
    )

Calculate outbreak thresholds from SEIR simulation results using outbreak specifications.

This method extracts incidence data from SEIR simulation results and delegates to the
core outbreak threshold calculation function. It provides a convenient interface when
working with structured SEIR simulation outputs.

# Arguments
- `seir_results::StructVector{SEIRRun}`: Structured vector of SEIR simulation results
- `outbreak_specification::OutbreakSpecification`: Configuration specifying outbreak criteria

# Returns
- `StructVector{OutbreakThresholds}`: Outbreak threshold information for each simulation

# See Also
- [`calculate_outbreak_thresholds(incidence_vecs, outbreak_specification)`](@ref): Core implementation
- [`SEIRRun`](@ref): Individual simulation result type
- [`OutbreakSpecification`](@ref): Outbreak detection configuration
"""
function calculate_outbreak_thresholds(
        seir_results::StructVector{SEIRRun},
        outbreak_specification::OutbreakSpecification
    )
    return calculate_outbreak_thresholds(
        seir_results.incidence,
        outbreak_specification
    )
end

"""
    calculate_outbreak_thresholds(
        incidence_vecs::VSV,
        outbreak_specification::OutbreakSpecification
    ) where {VSV <: AbstractVector{<:AbstractVector}}

Calculate outbreak thresholds from incidence timeseries using specified outbreak criteria.

This function processes multiple incidence timeseries to identify outbreak periods that
meet specified criteria for threshold, duration, and size. It applies threshold detection
to each timeseries and classifies periods as outbreaks based on the outbreak specification
parameters.

The detection process involves:
1. Identifying periods where incidence exceeds the outbreak threshold
2. Calculating bounds and durations of these periods
3. Classifying periods as outbreaks based on minimum duration and size criteria
4. Returning only those periods that qualify as true outbreaks

# Arguments
- `incidence_vecs::VSV`: Vector of incidence timeseries (each element is a vector of infection counts)
- `outbreak_specification::OutbreakSpecification`: Configuration containing:
  - `outbreak_threshold`: Minimum daily incidence to consider for outbreak detection
  - `minimum_outbreak_duration`: Minimum number of consecutive days for outbreak classification
  - `minimum_outbreak_size`: Minimum total infections during outbreak period

# Returns
- `StructVector{OutbreakThresholds}`: Structured vector where each element contains outbreak
  information for one simulation, with fields:
  - `lower_bounds`: Starting indices of outbreak periods
  - `upper_bounds`: Ending indices of outbreak periods
  - `duration`: Duration of each outbreak period
  - `num_infections_during_bounds`: Total infections during each outbreak period

# Examples
```julia
# Define outbreak criteria
outbreak_spec = OutbreakSpecification(
    outbreak_threshold = 5,
    minimum_outbreak_duration = 7,
    minimum_outbreak_size = 50
)

# Calculate outbreaks from incidence data
incidence_data = [rand(1:10, 100) for _ in 1:50]  # 50 simulations
outbreaks = calculate_outbreak_thresholds(incidence_data, outbreak_spec)

# Analyze results
println("Number of simulations with outbreaks: ",
        sum(length(outbreak.lower_bounds) > 0 for outbreak in outbreaks))
```

# See Also
- [`calculate_outbreak_thresholds!`](@ref): In-place version for performance
- [`OutbreakSpecification`](@ref): Configuration type for outbreak criteria
- [`OutbreakThresholds`](@ref): Return type containing outbreak information
- [`classify_all_outbreaks`](@ref): Function that classifies threshold periods as outbreaks
"""
function calculate_outbreak_thresholds(
        incidence_vecs::VSV,
        outbreak_specification::OutbreakSpecification
    ) where {VSV <: AbstractVector{<:AbstractVector}}

    emergent_outbreak_threshold_vecs = Vector{OutbreakThresholds}(
        undef, length(incidence_vecs)
    )

    calculate_outbreak_thresholds!(
        emergent_outbreak_threshold_vecs,
        incidence_vecs,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    return StructVector(emergent_outbreak_threshold_vecs)
end

"""
    calculate_outbreak_thresholds!(
        emergent_outbreak_threshold_vecs,
        incidence_vecs,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize
    )

In-place calculation of outbreak thresholds for performance optimization.

This function performs the core outbreak detection algorithm by processing each incidence
timeseries and storing results directly in the pre-allocated output vector. It uses a
reusable worker vector to minimize memory allocations during the threshold detection process.

The algorithm for each simulation:
1. Creates a boolean mask for incidence values above threshold
2. Applies run-length encoding to identify consecutive above-threshold periods
3. Calculates bounds and durations using [`calculate_above_threshold_bounds`](@ref)
4. Classifies periods as outbreaks using [`classify_all_outbreaks`](@ref)
5. Stores results in the pre-allocated output vector

# Arguments
- `emergent_outbreak_threshold_vecs`: Pre-allocated vector to store OutbreakThresholds results
- `incidence_vecs`: Vector of incidence timeseries to process
- `outbreakthreshold`: Minimum daily incidence for outbreak consideration
- `minoutbreakdur`: Minimum consecutive days for outbreak classification
- `minoutbreaksize`: Minimum total infections for outbreak classification

# Returns
- `nothing`: Results are stored in-place in `emergent_outbreak_threshold_vecs`

# Performance Notes
- Uses a single pre-allocated boolean worker vector to minimize allocations
- Processes simulations sequentially with `@inbounds` optimization
- Designed for high-performance batch processing of many simulations

# See Also
- [`calculate_outbreak_thresholds`](@ref): High-level interface that calls this function
- [`classify_all_outbreaks`](@ref): Function that classifies threshold periods
- [`calculate_above_threshold_bounds`](@ref): Function that calculates period bounds
"""
function calculate_outbreak_thresholds!(
        emergent_outbreak_threshold_vecs,
        incidence_vecs,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize,
    )
    above_threshold_worker_vec = Vector{Bool}(undef, length(incidence_vecs[1]))
    @inbounds for (sim, incidence_vec) in pairs(incidence_vecs)
        above_threshold_worker_vec .= incidence_vec .>= outbreakthreshold

        abovethresholdrle = StatsBase.rle(above_threshold_worker_vec)

        threshold_bounds = calculate_above_threshold_bounds(abovethresholdrle)

        local outbreak_thresholds = classify_all_outbreaks(
            incidence_vec,
            threshold_bounds,
            minoutbreakdur,
            minoutbreaksize,
        )
        emergent_outbreak_threshold_vecs[sim] = outbreak_thresholds
    end
    return nothing
end

"""
    classify_all_outbreaks(
        incidence_vec,
        incidence_thresholds::Thresholds,
        minoutbreakdur,
        minoutbreaksize
    )

Classify threshold-crossing periods as outbreaks based on duration and size criteria.

This function takes periods where incidence exceeded a threshold and determines which
of these periods qualify as true outbreaks. It applies additional criteria for minimum
duration and total infection count to filter out brief spikes or small increases that
don't constitute meaningful outbreak events.

For each threshold-crossing period, the function:
1. Calculates the total number of infections during the period
2. Applies outbreak classification criteria using [`classify_outbreak`](@ref)
3. Retains only periods that meet both duration and size requirements
4. Returns outbreak-specific threshold information

# Arguments
- `incidence_vec`: Vector of daily incidence counts for a single simulation
- `incidence_thresholds::Thresholds`: Threshold periods with bounds and durations
- `minoutbreakdur`: Minimum consecutive days required for outbreak classification
- `minoutbreaksize`: Minimum total infections required for outbreak classification

# Returns
- `OutbreakThresholds`: Filtered threshold information containing only true outbreaks:
  - `lower_bounds`: Starting indices of outbreak periods
  - `upper_bounds`: Ending indices of outbreak periods
  - `duration`: Duration of each outbreak period
  - `num_infections_during_bounds`: Total infections during each outbreak period

# Examples
```julia
# Example threshold periods from incidence analysis
thresholds = Thresholds([10, 50], [15, 60], [6, 11])
incidence = rand(1:20, 100)

# Classify with outbreak criteria
outbreaks = classify_all_outbreaks(incidence, thresholds, 7, 50)

# Check which periods qualified as outbreaks
println("Outbreaks detected: ", length(outbreaks.lower_bounds))
```

# See Also
- [`classify_outbreak`](@ref): Function that applies classification criteria to individual periods
- [`Thresholds`](@ref): Input type containing threshold period information
- [`OutbreakThresholds`](@ref): Return type with outbreak-specific information
"""
function classify_all_outbreaks(
        incidence_vec,
        incidence_thresholds::Thresholds,
        minoutbreakdur,
        minoutbreaksize,
    )
    outbreak_status_vec = similar(incidence_thresholds.lower_bounds)
    ninf_during_bounds_vec = similar(incidence_thresholds.lower_bounds)

    for i in eachindex(incidence_thresholds.lower_bounds)
        @inbounds begin
            local lower_bound = incidence_thresholds.lower_bounds[i]
            local upper_bound = incidence_thresholds.upper_bounds[i]
            local duration = incidence_thresholds.duration[i]
            local ninf_during_bounds = sum(
                @view(incidence_vec[lower_bound:upper_bound])
            )
        end

        local outbreak_status = classify_outbreak(
            ninf_during_bounds,
            duration,
            minoutbreakdur,
            minoutbreaksize,
        )

        outbreak_status_vec[i] = outbreak_status
        ninf_during_bounds_vec[i] = ninf_during_bounds

    end

    outbreak_idx = findall(isequal(1), outbreak_status_vec)
    num_outbreaks = length(outbreak_idx)

    return OutbreakThresholds(
        Vector{Int64}(incidence_thresholds.lower_bounds[outbreak_idx]),
        Vector{Int64}(incidence_thresholds.upper_bounds[outbreak_idx]),
        Vector{Int64}(incidence_thresholds.duration[outbreak_idx]),
        Vector{Int64}(ninf_during_bounds_vec[outbreak_idx])
    )
end

"""
    classify_outbreak(
        ninf_during_bounds,
        duration,
        minoutbreakdur,
        minoutbreaksize
    )

Classify a single threshold period as an outbreak based on duration and size criteria.

This function applies the binary classification logic to determine whether a period
of above-threshold incidence qualifies as a true outbreak. A period is classified
as an outbreak if it meets both the minimum duration requirement (consecutive days)
and the minimum size requirement (total infections).

# Arguments
- `ninf_during_bounds`: Total number of infections during the threshold period
- `duration`: Number of consecutive days in the threshold period
- `minoutbreakdur`: Minimum consecutive days required for outbreak classification
- `minoutbreaksize`: Minimum total infections required for outbreak classification

# Returns
- `Int`: Binary classification result:
  - `1`: Period qualifies as an outbreak (meets both duration and size criteria)
  - `0`: Period does not qualify as an outbreak

# Examples
```julia
# Example: Period with 45 infections over 10 days
is_outbreak = classify_outbreak(45, 10, 7, 40)  # Returns 1 (outbreak)

# Example: Period with 30 infections over 5 days
is_outbreak = classify_outbreak(30, 5, 7, 40)   # Returns 0 (too short)

# Example: Period with 35 infections over 8 days
is_outbreak = classify_outbreak(35, 8, 7, 40)   # Returns 0 (too few infections)
```

# See Also
- [`classify_all_outbreaks`](@ref): Function that applies this to multiple periods
- [`OutbreakSpecification`](@ref): Type that defines the classification criteria
"""
function classify_outbreak(
        ninf_during_bounds,
        duration,
        minoutbreakdur,
        minoutbreaksize
    )
    if duration >= minoutbreakdur && ninf_during_bounds >= minoutbreaksize
        return 1
    end
    return 0
end
