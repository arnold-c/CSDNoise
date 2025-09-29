using StructArrays: StructVector
using FixedSizeArrays: FixedSizeVector
using StatsBase: mean
using Bumper: @no_escape, @alloc
using Random: Random

"""
    create_testing_vecs(
        seir_results::StructVector{SEIRRun},
        noise_results::NoiseRun,
        perc_tested::Float64,
        individual_test_spec::IndividualTestSpecification,
    ) -> Vector{FixedSizeVector{Int64}}

Create diagnostic testing vectors from StructVector SEIR results and NoiseRun data.

This function processes variable-length simulation data to calculate the total number
of test-positive individuals for each simulation, accounting for both true positives
from infected individuals and false positives from noise.

# Arguments
- `seir_results`: StructVector of SEIRRun containing filtered incidence data
- `noise_results`: NoiseRun with variable-length incidence vectors
- `perc_tested`: Percentage of individuals tested (0.0 to 1.0)
- `individual_test_spec`: Test specifications including sensitivity, specificity, and lag

# Returns
- `Vector{FixedSizeVector{Int64}}`: Vector of test-positive counts, one per simulation

# Details
The function:
1. Validates that SEIR and noise incidence vectors have matching lengths for each simulation
2. Calculates the number of infected and noise individuals tested
3. Applies test sensitivity to infected individuals (true positives)
4. Applies test specificity to noise individuals (false positives)
5. Accounts for test result lag
6. Returns total test-positive individuals per simulation

Uses Bumper.jl for efficient temporary allocations to minimize heap usage.

# Example
```julia
test_results = create_testing_vecs(
    seir_results,
    noise_results,
    0.1,  # 10% testing rate
    IndividualTestSpecification(0.9, 0.95, 2)  # 90% sensitive, 95% specific, 2-day lag
)
```
"""
function create_testing_vecs(
        seir_results::StructVector{SEIRRun},
        noise_results::NoiseRun,
        perc_tested::Float64,
        individual_test_spec::IndividualTestSpecification,
    )
    nsims = length(seir_results)
    @assert nsims == length(noise_results.incidence) "Number of SEIR simulations must match noise simulations"
    @assert 0.0 <= perc_tested <= 1.0 "perc_tested must be between 0.0 and 1.0"

    # Pre-allocate result vector
    test_results = Vector{FixedSizeVector{Int64}}(undef, nsims)

    # Extract test parameters
    test_lag = individual_test_spec.test_result_lag
    sensitivity = individual_test_spec.sensitivity
    specificity = individual_test_spec.specificity

    @assert test_lag >= 0 "lag must be non-negative"
    @assert 0.0 <= sensitivity <= 1.0 "sensitivity must be between 0.0 and 1.0"
    @assert 0.0 <= specificity <= 1.0 "specificity must be between 0.0 and 1.0"

    for sim in eachindex(test_results)
        seir_incidence = seir_results.incidence[sim]
        noise_incidence = noise_results.incidence[sim]
        sim_length = length(seir_incidence)

        # Validate matching lengths
        @assert length(noise_incidence) == sim_length "SEIR and noise incidence lengths must match for simulation $sim"

        # Use Bumper for temporary allocations
        @no_escape begin
            # Temporary vectors for calculations
            seir_tested = @alloc(Int64, sim_length)
            noise_tested = @alloc(Int64, sim_length)
            true_positives = @alloc(Int64, sim_length)
            false_positives = @alloc(Int64, sim_length)

            # Calculate number tested from each source
            calculate_tested_vec!(seir_tested, seir_incidence, perc_tested)
            calculate_tested_vec!(noise_tested, noise_incidence, perc_tested)

            # Calculate true positives (from infected individuals)
            calculate_positives_vec!(true_positives, seir_tested, sim_length, test_lag, sensitivity)

            # Calculate false positives (from noise individuals)
            calculate_positives_vec!(false_positives, noise_tested, sim_length, test_lag, 1.0 - specificity)

            # Create result vector with total positives
            total_positives = FixedSizeVector{Int64}(undef, sim_length)
            @inbounds for i in 1:sim_length
                total_positives[i] = true_positives[i] + false_positives[i]
            end

            test_results[sim] = total_positives
        end
    end

    return test_results
end

"""
    calculate_tested_vec!(
        outvec::AbstractVector{Int64},
        invec::FixedSizeVector{Int64},
        perc_tested::Float64
    )

Calculate the number of individuals tested from incidence data.

Applies the testing percentage to each day's incidence, rounding to nearest integer.

# Arguments
- `outvec`: Output vector to store number tested (modified in-place)
- `invec`: Input incidence vector
- `perc_tested`: Percentage of individuals tested (0.0 to 1.0)
"""
function calculate_tested_vec!(
        outvec::AbstractVector{Int64},
        invec::FixedSizeVector{Int64},
        perc_tested::Float64
    )
    @inbounds for i in eachindex(invec)
        outvec[i] = round(Int64, invec[i] * perc_tested)
    end

    return nothing
end

"""
    calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64
    )

Calculate test-positive individuals accounting for test result lag.

Applies the test performance multiplier (sensitivity or 1-specificity) to the number
tested, with results appearing after the specified lag period.

# Arguments
- `npos_vec`: Output vector for positive test results (modified in-place)
- `tested_vec`: Input vector of individuals tested
- `sim_length`: Length of simulation
- `lag`: Test result lag in days
- `tested_multiplier`: Test performance multiplier (sensitivity or 1-specificity)
"""
function calculate_positives_vec!(
        npos_vec::AbstractVector{Int64},
        tested_vec::AbstractVector{Int64},
        sim_length::Int64,
        lag::Int64,
        tested_multiplier::Float64
    )
    # Initialize to zero
    fill!(npos_vec, 0)

    @inbounds for day in eachindex(npos_vec)
        result_day = day + lag
        if result_day <= sim_length
            npos_vec[result_day] = round(Int64, tested_vec[day] * tested_multiplier)
        end
    end

    return nothing
end
