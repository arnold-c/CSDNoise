export create_test_positive_vecs

"""
    create_test_positive_vecs(
        ensemble_results::StructVector{EnsembleSEIRRun},
        noise_results::NoiseRun,
        percent_tested::Float64,
        test_specification::IndividualTestSpecification,
    ) -> EnsembleTestResultRun

Create diagnostic testing vectors for both emergent and null scenarios from ensemble SEIR results.

This function processes ensemble simulation data containing both emergent (with outbreak)
and null (without outbreak) scenarios, generating test-positive counts for each scenario
using the same noise realization. This allows for direct comparison between scenarios
under identical noise conditions.

# Arguments
- `ensemble_results`: StructVector of EnsembleSEIRRun containing paired emergent and null simulations
- `noise_results`: NoiseRun with variable-length incidence vectors shared across scenarios
- `percent_tested`: Percentage of individuals tested (0.0 to 1.0)
- `test_specification`: Test specifications including sensitivity, specificity, and lag

# Returns
- `EnsembleTestResultRun`: Container with separate test-positive vectors for emergent and null scenarios

# Details
The function:
1. Validates that ensemble and noise simulations have matching counts
2. Validates test parameters (percent_tested, sensitivity, specificity, lag)
3. Processes each simulation pair (emergent and null) with the same noise realization
4. Calculates test-positive counts accounting for true and false positives
5. Returns paired results for direct scenario comparison

# Example
```julia
test_results = create_test_positive_vecs(
    ensemble_results,
    noise_results,
    0.1,  # 10% testing rate
    IndividualTestSpecification(0.9, 0.95, 2)  # 90% sensitive, 95% specific, 2-day lag
)
```
"""
function create_test_positive_vecs(
        ensemble_results::EnsembleSEIRRun,
        noise_results::NoiseRun,
        percent_tested::Float64,
        test_specification::IndividualTestSpecification,
    )

    nsims = length(ensemble_results.emergent_seir_run)
    @assert nsims == length(noise_results.incidence) "Number of SEIR simulations must match noise simulations"
    @assert 0.0 <= percent_tested <= 1.0 "percent_tested must be between 0.0 and 1.0"

    # Pre-allocate result vectors for both scenarios
    emergent_test_positives = Vector{Vector{Int64}}(undef, nsims)
    null_test_positives = Vector{Vector{Int64}}(undef, nsims)

    # Validate test parameters
    @assert test_specification.test_result_lag >= 0 "lag must be non-negative"
    @assert 0.0 <= test_specification.sensitivity <= 1.0 "sensitivity must be between 0.0 and 1.0"
    @assert 0.0 <= test_specification.specificity <= 1.0 "specificity must be between 0.0 and 1.0"

    # Process each simulation pair with shared noise
    for sim in eachindex(emergent_test_positives)
        emergent_seir_incidence = ensemble_results.emergent_seir_run.incidence[sim]
        null_seir_incidence = ensemble_results.null_seir_run.incidence[sim]
        noise_incidence = noise_results.incidence[sim]

        # Generate test positives for emergent scenario
        _create_test_positive_vec!(
            emergent_test_positives,
            emergent_seir_incidence,
            noise_incidence,
            percent_tested,
            test_specification,
            sim
        )

        # Generate test positives for null scenario using same noise
        _create_test_positive_vec!(
            null_test_positives,
            null_seir_incidence,
            noise_incidence,
            percent_tested,
            test_specification,
            sim
        )
    end

    return EnsembleTestResultRun(
        emergent_test_positives = emergent_test_positives,
        null_test_positives = null_test_positives
    )
end


"""
    create_test_positive_vecs(
        seir_results::StructVector{SEIRRun},
        noise_results::NoiseRun,
        percent_tested::Float64,
        test_specification::IndividualTestSpecification,
    ) -> Vector{Vector{Int64}}

Create diagnostic testing vectors from StructVector SEIR results and NoiseRun data.

This function processes variable-length simulation data to calculate the total number
of test-positive individuals for each simulation, accounting for both true positives
from infected individuals and false positives from noise.

# Arguments
- `seir_results`: StructVector of SEIRRun containing filtered incidence data
- `noise_results`: NoiseRun with variable-length incidence vectors
- `percent_tested`: Percentage of individuals tested (0.0 to 1.0)
- `test_specification`: Test specifications including sensitivity, specificity, and lag

# Returns
- `Vector{Vector{Int64}}`: Vector of test-positive counts, one per simulation

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
test_results = create_test_positive_vecs(
    seir_results,
    noise_results,
    0.1,  # 10% testing rate
    IndividualTestSpecification(0.9, 0.95, 2)  # 90% sensitive, 95% specific, 2-day lag
)
```
"""
function create_test_positive_vecs(
        seir_results::StructVector{SEIRRun},
        noise_results::NoiseRun,
        percent_tested::Float64,
        test_specification::IndividualTestSpecification,
    )
    nsims = length(seir_results)
    @assert nsims == length(noise_results.incidence) "Number of SEIR simulations must match noise simulations"
    @assert 0.0 <= percent_tested <= 1.0 "percent_tested must be between 0.0 and 1.0"

    # Pre-allocate result vector
    test_results = Vector{Vector{Int64}}(undef, nsims)

    # Extract test parameters
    @assert test_specification.test_result_lag >= 0 "lag must be non-negative"
    @assert 0.0 <= test_specification.sensitivity <= 1.0 "sensitivity must be between 0.0 and 1.0"
    @assert 0.0 <= test_specification.specificity <= 1.0 "specificity must be between 0.0 and 1.0"

    for sim in eachindex(test_results)
        seir_incidence = seir_results.incidence[sim]
        noise_incidence = noise_results.incidence[sim]

        _create_test_positive_vec!(
            test_results,
            seir_incidence,
            noise_incidence,
            percent_tested,
            test_specification,
            sim
        )
    end

    return test_results
end


"""
    _create_test_positive_vec!(
        test_results::Vector{Vector{Int64}},
        seir_incidence::Vector{Int64},
        noise_incidence::Vector{Int64},
        percent_tested::Float64,
        test_specification::IndividualTestSpecification,
        sim::Int64
    )

Internal helper function to calculate test-positive counts for a single simulation.

This function performs the core calculation of diagnostic test results by combining
true positives from infected individuals and false positives from noise individuals,
accounting for test characteristics and result lag. Uses Bumper.jl for efficient
memory management with temporary allocations.

# Arguments
- `test_results`: Pre-allocated vector to store results for all simulations
- `seir_incidence`: Daily incidence counts from SEIR simulation
- `noise_incidence`: Daily noise counts for the same time period
- `percent_tested`: Fraction of individuals tested each day (0.0 to 1.0)
- `test_specification`: Test parameters (sensitivity, specificity, lag)
- `sim`: Simulation index for storing results

# Details
The function:
1. Validates that SEIR and noise incidence vectors have matching lengths
2. Allocates temporary vectors using Bumper.jl for memory efficiency
3. Calculates number of individuals tested from each source (SEIR and noise)
4. Applies test sensitivity to infected individuals (true positives)
5. Applies false positive rate (1 - specificity) to noise individuals
6. Accounts for test result lag in both calculations
7. Combines true and false positives into total test-positive counts
8. Stores results in the appropriate simulation slot

# Performance Notes
- Uses `@no_escape` and `@alloc` from Bumper.jl for stack-allocated temporaries
- Employs `@inbounds` for optimized array access in the final summation loop
- Minimizes heap allocations by reusing temporary vectors

This is an internal function and should not be called directly by users.
"""
function _create_test_positive_vec!(
        test_results::Vector{Vector{Int64}},
        seir_incidence::Vector{Int64},
        noise_incidence::Vector{Int64},
        percent_tested::Float64,
        test_specification::IndividualTestSpecification,
        sim::Int64
    )
    sim_length = length(seir_incidence)

    # Validate matching lengths
    @assert length(noise_incidence) == sim_length "SEIR and noise incidence lengths must match for simulation $sim"

    # Use Bumper for temporary allocations
    return @no_escape begin
        # Temporary vectors for calculations
        seir_tested = @alloc(Int64, sim_length)
        noise_tested = @alloc(Int64, sim_length)
        true_positives = @alloc(Int64, sim_length)
        false_positives = @alloc(Int64, sim_length)

        # Calculate number tested from each source
        calculate_tested_vec!(seir_tested, seir_incidence, percent_tested)
        calculate_tested_vec!(noise_tested, noise_incidence, percent_tested)

        # Calculate true positives (from infected individuals)
        calculate_positives_vec!(
            true_positives,
            seir_tested,
            sim_length,
            test_specification.test_result_lag,
            test_specification.sensitivity
        )

        # Calculate false positives (from noise individuals)
        calculate_positives_vec!(
            false_positives,
            noise_tested,
            sim_length,
            test_specification.test_result_lag,
            1.0 - test_specification.specificity
        )

        # Create result vector with total positives
        total_positives = Vector{Int64}(undef, sim_length)
        @inbounds for i in 1:sim_length
            total_positives[i] = true_positives[i] + false_positives[i]
        end

        test_results[sim] = total_positives
    end
end
