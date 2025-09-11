#!/usr/bin/env julia
"""
Test type stability using JET to ensure no runtime dispatches
"""

using Test
using JET
using DrWatson
@quickactivate("CSDNoise")
using CSDNoise

@testset "Type Stability Tests" begin
    @testset "generate_ensemble_data type stability" begin
        # Create small test data to minimize test time
        ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(5)

        # Test that generate_ensemble_data has no runtime dispatches
        @test_opt target_modules = (CSDNoise,) generate_ensemble_data(
            ensemble_spec, null_spec, outbreak_spec
        )
    end

    @testset "calculate_beta_amp type stability" begin
        # Test individual function
        @test_opt target_modules = (CSDNoise,) calculate_beta_amp(
            1.0, 0.1, 1.0; seasonality = cos
        )
    end

    @testset "DynamicsParameters constructor type stability" begin
        # Create test specification
        spec = DynamicsParameterSpecification(
            1.0,      # beta_mean
            0.1,      # beta_force
            cos,      # seasonality
            1 / 4.0,    # sigma
            1 / 5.0,    # gamma
            1 / (70 * 365), # mu
            20,       # annual_births_per_k
            0.0,      # epsilon
            2.0,      # R_0
            0.8,      # min_burnin_vaccination_coverage
            0.9,      # max_burnin_vaccination_coverage
            0.7,      # min_vaccination_coverage
            0.8       # max_vaccination_coverage
        )

        @test_opt target_modules = (CSDNoise,) DynamicsParameters(spec; seed = 42)
    end
end
