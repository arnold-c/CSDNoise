export IndividualTestSpecification

"""
    IndividualTestSpecification

Diagnostic test characteristics for individual-level testing in disease surveillance.

This struct defines the performance characteristics of a diagnostic test used to detect
infected individuals in the population. These parameters are used in conjunction with
testing coverage rates to simulate realistic diagnostic testing scenarios and calculate
test-positive individuals accounting for both true and false positives.

# Fields
- `sensitivity::Float64`: Test sensitivity (true positive rate), the probability that the test correctly identifies an infected individual (0.0 to 1.0)
- `specificity::Float64`: Test specificity (true negative rate), the probability that the test correctly identifies a non-infected individual (0.0 to 1.0)
- `test_result_lag::Int64`: Number of days between test administration and result availability (≥ 0)

# Constructor
    IndividualTestSpecification(; sensitivity, specificity, test_result_lag)

# Example
```julia
# Perfect test with no lag
perfect_test = IndividualTestSpecification(
    sensitivity = 1.0,
    specificity = 1.0,
    test_result_lag = 0
)

# Realistic PCR test with 2-day processing time
pcr_test = IndividualTestSpecification(
    sensitivity = 0.95,
    specificity = 0.99,
    test_result_lag = 2
)
```

# See Also
- [`create_test_positive_vecs`](@ref): Uses test specifications to calculate test-positive individuals
"""
Base.@kwdef struct IndividualTestSpecification
    sensitivity::Float64
    specificity::Float64
    test_result_lag::Int64
end
