using Test
using CSDNoise
using StructArrays: StructVector
using StaticArrays: SVector

@testset "Diagnostic Testing StructVector Functions" begin

    @testset "calculate_tested_vec!" begin
        # Test basic functionality
        invec = Vector{Int64}([10, 20, 30, 40])
        outvec = Vector{Int64}(undef, 4)

        calculate_tested_vec!(outvec, invec, 0.1)
        @test outvec == [1, 2, 3, 4]  # 10% of each value, rounded

        # Test with zero percentage
        calculate_tested_vec!(outvec, invec, 0.0)
        @test outvec == [0, 0, 0, 0]

        # Test with 100% percentage
        calculate_tested_vec!(outvec, invec, 1.0)
        @test outvec == [10, 20, 30, 40]
    end

    @testset "calculate_positives_vec!" begin
        tested_vec = [10, 20, 30, 40]
        npos_vec = Vector{Int64}(undef, 4)

        # Test with no lag
        calculate_positives_vec!(npos_vec, tested_vec, 4, 0, 0.9)
        @test npos_vec == [9, 18, 27, 36]  # 90% of each value

        # Test with 1-day lag
        fill!(npos_vec, 0)
        calculate_positives_vec!(npos_vec, tested_vec, 4, 1, 0.9)
        @test npos_vec == [0, 9, 18, 27]  # Results shifted by 1 day

        # Test with 2-day lag (last result falls off)
        fill!(npos_vec, 0)
        calculate_positives_vec!(npos_vec, tested_vec, 4, 2, 0.9)
        @test npos_vec == [0, 0, 9, 18]  # Results shifted by 2 days

        # Test with lag longer than simulation
        fill!(npos_vec, 0)
        calculate_positives_vec!(npos_vec, tested_vec, 4, 5, 0.9)
        @test npos_vec == [0, 0, 0, 0]  # All results fall off

        # Test with zero multiplier
        fill!(npos_vec, 0)
        calculate_positives_vec!(npos_vec, tested_vec, 4, 0, 0.0)
        @test npos_vec == [0, 0, 0, 0]
    end

    @testset "create_testing_vecs" begin
        # Create test data
        nsims = 3
        sim_lengths = [5, 4, 6]  # Variable lengths

        # Create SEIR results
        seir_incidence = [
            Vector{Int64}([10, 15, 20, 25, 30]),  # sim 1: length 5
            Vector{Int64}([5, 10, 15, 20]),       # sim 2: length 4
            Vector{Int64}([8, 12, 16, 20, 24, 28]), # sim 3: length 6
        ]

        seir_states = [
            Vector{SVector{5, Int64}}([SVector(1000, 10, 5, 985, 1000) for _ in 1:5]),
            Vector{SVector{5, Int64}}([SVector(1000, 8, 4, 988, 1000) for _ in 1:4]),
            Vector{SVector{5, Int64}}([SVector(1000, 12, 6, 982, 1000) for _ in 1:6]),
        ]

        seir_Reff = [
            Vector{Float64}([1.2, 1.1, 1.0, 0.9, 0.8]),
            Vector{Float64}([1.3, 1.2, 1.1, 1.0]),
            Vector{Float64}([1.1, 1.0, 0.9, 0.8, 0.7, 0.6]),
        ]

        seir_results = StructVector{SEIRRun}(
            states = seir_states,
            incidence = seir_incidence,
            Reff = seir_Reff
        )

        # Create noise results
        noise_incidence = [
            Vector{Int64}([2, 3, 4, 5, 6]),      # sim 1: length 5
            Vector{Int64}([1, 2, 3, 4]),         # sim 2: length 4
            Vector{Int64}([3, 4, 5, 6, 7, 8]),    # sim 3: length 6
        ]

        noise_results = NoiseRun(
            noise_incidence,
            5.0,  # mean_noise
            3.0,  # mean_poisson_noise
            2.0   # mean_dynamic_noise
        )

        # Test specifications
        test_spec = IndividualTestSpecification(0.9, 0.95, 1)  # 90% sensitive, 95% specific, 1-day lag
        perc_tested = 0.1  # 10% testing rate

        # Run function
        results = create_testing_vecs(seir_results, noise_results, perc_tested, test_spec)

        # Validate results structure
        @test length(results) == nsims
        @test length(results[1]) == 5  # sim 1 length
        @test length(results[2]) == 4  # sim 2 length
        @test length(results[3]) == 6  # sim 3 length

        # Validate specific calculations for sim 1
        # Day 1: SEIR tested = 1, noise tested = 0, but results appear on day 2 due to lag
        # Day 2: True positives from day 1 = round(1 * 0.9) = 1, False positives from day 1 = round(0 * 0.05) = 0
        # Plus new tests: SEIR tested = 2, noise tested = 0
        @test results[1][1] == 0  # No results yet due to lag
        @test results[1][2] == 1  # True positives from day 1: round(1 * 0.9) + round(0 * 0.05) = 1 + 0 = 1

        # Test error conditions
        # Mismatched simulation counts
        short_noise = NoiseRun([Vector{Int64}([1, 2])], 1.0, 1.0, 0.0)
        @test_throws AssertionError create_testing_vecs(seir_results, short_noise, perc_tested, test_spec)

        # Invalid percentage
        @test_throws AssertionError create_testing_vecs(seir_results, noise_results, 1.5, test_spec)
        @test_throws AssertionError create_testing_vecs(seir_results, noise_results, -0.1, test_spec)

        # Invalid test specifications
        bad_test_spec = IndividualTestSpecification(1.5, 0.95, 1)  # Invalid sensitivity
        @test_throws AssertionError create_testing_vecs(seir_results, noise_results, perc_tested, bad_test_spec)

        bad_test_spec2 = IndividualTestSpecification(0.9, -0.1, 1)  # Invalid specificity
        @test_throws AssertionError create_testing_vecs(seir_results, noise_results, perc_tested, bad_test_spec2)

        # Mismatched incidence lengths within simulation
        bad_noise_incidence = [
            Vector{Int64}([2, 3, 4]),  # Wrong length for sim 1
            Vector{Int64}([1, 2, 3, 4]),
            Vector{Int64}([3, 4, 5, 6, 7, 8]),
        ]
        bad_noise_results = NoiseRun(bad_noise_incidence, 5.0, 3.0, 2.0)
        @test_throws AssertionError create_testing_vecs(seir_results, bad_noise_results, perc_tested, test_spec)
    end

    @testset "Edge Cases" begin
        # Test with single time point
        seir_incidence = [Vector{Int64}([10])]
        seir_states = [Vector{SVector{5, Int64}}([SVector(1000, 10, 5, 985, 1000)])]
        seir_Reff = [Vector{Float64}([1.2])]

        seir_results = StructVector{SEIRRun}(
            states = seir_states,
            incidence = seir_incidence,
            Reff = seir_Reff
        )

        noise_incidence = [Vector{Int64}([2])]
        noise_results = NoiseRun(noise_incidence, 2.0, 2.0, 0.0)

        test_spec = IndividualTestSpecification(1.0, 1.0, 0)  # Perfect test, no lag
        results = create_testing_vecs(seir_results, noise_results, 0.1, test_spec)

        @test length(results) == 1
        @test length(results[1]) == 1
        @test results[1][1] == 1  # round(10 * 0.1 * 1.0) + round(2 * 0.1 * 0.0) = 1 + 0 = 1

        # Test with zero incidence
        zero_seir_incidence = [Vector{Int64}([0, 0, 0])]
        zero_seir_states = [Vector{SVector{5, Int64}}([SVector(1000, 0, 0, 1000, 1000) for _ in 1:3])]
        zero_seir_Reff = [Vector{Float64}([0.0, 0.0, 0.0])]

        zero_seir_results = StructVector{SEIRRun}(
            states = zero_seir_states,
            incidence = zero_seir_incidence,
            Reff = zero_seir_Reff
        )

        zero_noise_incidence = [Vector{Int64}([0, 0, 0])]
        zero_noise_results = NoiseRun(zero_noise_incidence, 0.0, 0.0, 0.0)

        zero_results = create_testing_vecs(zero_seir_results, zero_noise_results, 0.1, test_spec)
        @test all(x -> all(y -> y == 0, x), zero_results)  # All zeros
    end
end
