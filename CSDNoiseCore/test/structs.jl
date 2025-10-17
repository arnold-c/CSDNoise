using Test
using CSDNoise

@testset "DynamicalNoiseSpecification" begin
    @testset "Constructor with all parameters" begin
        spec = DynamicalNoiseSpecification(
            R0 = 5.0,
            latent_period = 7,
            duration_infection = 14,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0, 1.0],
            susceptible_bounds = [0.01, 0.99],
            max_vaccination_range = 0.2
        )

        @test spec.R0 == 5.0
        @test spec.latent_period == 7
        @test spec.duration_infection == 14
        @test spec.correlation == "in-phase"
        @test spec.poisson_component == 1.0
        @test spec.vaccination_bounds == [0.0, 1.0]
        @test spec.susceptible_bounds == [0.01, 0.99]
        @test spec.max_vaccination_range == 0.2
    end

    @testset "Constructor with default parameters" begin
        spec = DynamicalNoiseSpecification(
            R0 = 12.0,
            latent_period = 8,
            duration_infection = 7,
            correlation = "out-of-phase",
            poisson_component = 0.5
        )

        @test spec.R0 == 12.0
        @test spec.latent_period == 8
        @test spec.duration_infection == 7
        @test spec.correlation == "out-of-phase"
        @test spec.poisson_component == 0.5
        @test spec.vaccination_bounds == [0.0, 1.0]
        @test spec.susceptible_bounds == [0.01, 0.99]
        @test spec.max_vaccination_range == 0.2
    end

    @testset "Constructor validates bounds" begin
        @test_throws AssertionError DynamicalNoiseSpecification(
            R0 = 5.0,
            latent_period = 7,
            duration_infection = 14,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [0.0]  # Wrong length
        )

        @test_throws AssertionError DynamicalNoiseSpecification(
            R0 = 5.0,
            latent_period = 7,
            duration_infection = 14,
            correlation = "in-phase",
            poisson_component = 1.0,
            vaccination_bounds = [1.0, 0.0]  # Wrong order
        )

        @test_throws AssertionError DynamicalNoiseSpecification(
            R0 = 5.0,
            latent_period = 7,
            duration_infection = 14,
            correlation = "in-phase",
            poisson_component = 1.0,
            susceptible_bounds = [0.99, 0.01]  # Wrong order
        )
    end
end

@testset "NoiseVaccinationOptimizationParameters" begin
    @testset "Constructor with all parameters" begin
        params = NoiseVaccinationOptimizationParameters(
            n_sobol_points = 50,
            local_algorithm = nothing,
            maxeval = 500,
            xtol_rel = 1.0e-4,
            xtol_abs = 1.0e-4,
            ftol_rel = 1.0e-5,
            verbose = true
        )

        @test params.n_sobol_points == 50
        @test isnothing(params.local_algorithm)
        @test params.maxeval == 500
        @test params.xtol_rel == 1.0e-4
        @test params.xtol_abs == 1.0e-4
        @test params.ftol_rel == 1.0e-5
        @test params.verbose == true
    end

    @testset "Constructor with default parameters" begin
        params = NoiseVaccinationOptimizationParameters()

        @test params.n_sobol_points == 100
        @test isnothing(params.local_algorithm)
        @test params.maxeval == 1000
        @test params.xtol_rel == 1.0e-3
        @test params.xtol_abs == 1.0e-3
        @test params.ftol_rel == 1.0e-4
        @test params.verbose == false
    end

    @testset "Constructor with partial parameters" begin
        params = NoiseVaccinationOptimizationParameters(
            n_sobol_points = 200,
            verbose = true
        )

        @test params.n_sobol_points == 200
        @test params.verbose == true
        @test params.maxeval == 1000  # Default value
    end
end
