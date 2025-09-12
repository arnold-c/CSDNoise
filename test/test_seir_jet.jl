using Test
using JET
using DrWatson
@quickactivate("CSDNoise")
using CSDNoise
using StaticArrays

@testset "SEIR Model Type Stability Tests" begin
    # Create minimal test data
    ensemble_spec, _, _ = create_ensemble_specs(2)
    state_params = ensemble_spec.state_parameters
    dynamics_spec = ensemble_spec.dynamics_parameter_specification
    time_params = ensemble_spec.time_parameters

    # Create dynamics parameters
    @test_opt target_modules = (CSDNoise,) DynamicsParameters(dynamics_spec; seed = 42)
    dynamics_params = DynamicsParameters(dynamics_spec; seed = 42)

    # Test data setup
    initial_states = SVector(state_params.init_states)

    @testset "seir_mod main function type stability" begin
        @test_opt target_modules = (CSDNoise,) seir_mod(
            initial_states,
            dynamics_params,
            time_params;
            seed = 42
        )
    end

    @testset "seir_mod! in-place function type stability" begin
        # Pre-allocate arrays as done in seir_mod
        state_vec = Vector{typeof(initial_states)}(undef, time_params.tlength)
        beta_vec = Vector{Float64}(undef, time_params.tlength)
        inc_vec = Vector{SVector{1, Int64}}(undef, time_params.tlength)

        @test_opt target_modules = (CSDNoise,) seir_mod!(
            state_vec,
            inc_vec,
            beta_vec,
            initial_states,
            dynamics_params,
            time_params;
            seed = 42
        )
    end

    @testset "seir_mod_loop! inner loop type stability" begin
        # Extract parameters for testing the inner loop
        beta_t = 1.5
        mu = dynamics_params.mu
        epsilon = dynamics_params.epsilon
        sigma = dynamics_params.sigma
        gamma = dynamics_params.gamma
        R_0 = dynamics_params.R_0
        vaccination_coverage = dynamics_params.vaccination_coverage
        timestep = time_params.tstep

        @test_opt target_modules = (CSDNoise,) seir_mod_loop!(
            initial_states,
            beta_t,
            mu,
            epsilon,
            sigma,
            gamma,
            R_0,
            vaccination_coverage,
            timestep
        )
    end

    @testset "calculate_beta_amp type stability" begin
        beta_mean = dynamics_params.beta_mean
        beta_force = dynamics_params.beta_force
        t = 100.0

        @test_opt target_modules = (CSDNoise,) calculate_beta_amp(
            beta_mean, beta_force, t; seasonality = cos
        )
    end

    @testset "Array conversion utilities type stability" begin
        # Test convert_svec_to_matrix
        test_svec = [SVector(1, 2, 3, 4, 5), SVector(2, 3, 4, 5, 6)]

        @test_opt target_modules = (CSDNoise,) convert_svec_to_matrix(test_svec)

        # Test convert_svec_to_matrix! in-place version
        test_matrix = Matrix{Int64}(undef, length(test_svec), length(test_svec[1]))
        @test_opt target_modules = (CSDNoise,) convert_svec_to_matrix!(test_matrix, test_svec)

        # Test convert_svec_to_array for 3D arrays
        test_svec_3d = [SVector(1, 2, 3, 4, 5) SVector(2, 3, 4, 5, 6); SVector(3, 4, 5, 6, 7) SVector(4, 5, 6, 7, 8)]
        @test_opt target_modules = (CSDNoise,) convert_svec_to_array(test_svec_3d)
    end

    @testset "Vector allocation type stability" begin
        # Test the vector allocations used in seir_mod
        tlength = time_params.tlength

        # Test state vector allocation
        @test_opt Vector{typeof(initial_states)}(undef, tlength)

        # Test beta vector allocation
        @test_opt Vector{Float64}(undef, tlength)

        # Test incidence vector allocation
        @test_opt Vector{SVector{1, Int64}}(undef, tlength)
    end

    @testset "Parameter extraction type stability" begin
        # Test that parameter extraction from structs is type stable
        @test_opt dynamics_params.mu
        @test_opt dynamics_params.epsilon
        @test_opt dynamics_params.sigma
        @test_opt dynamics_params.gamma
        @test_opt dynamics_params.R_0
        @test_opt dynamics_params.beta_mean
        @test_opt dynamics_params.beta_force
        @test_opt dynamics_params.seasonality
        @test_opt dynamics_params.burnin_vaccination_coverage
        @test_opt dynamics_params.vaccination_coverage

        @test_opt time_params.tstep
        @test_opt time_params.trange
        @test_opt time_params.burnin
        @test_opt time_params.tlength
    end

    @testset "SVector operations type stability" begin
        # Test SVector construction and indexing
        test_state = SVector(100, 10, 5, 85, 200)

        @test_opt test_state[1]
        @test_opt test_state[2]
        @test_opt test_state[3]
        @test_opt test_state[4]
        @test_opt test_state[5]

        # Test SVector arithmetic
        delta = SVector(1, -1, 0, 1, 1)
        @test_opt test_state + delta

        # Test SVector construction from components
        @test_opt SVector(test_state[1] + 1, test_state[2] - 1, test_state[3], test_state[4] + 1, test_state[5] + 1)

        # Test single-element SVector for incidence
        @test_opt SVector(5)
        @test_opt SVector(0)
    end

    @testset "Random number generation type stability" begin
        # Test the random number generation calls used in seir_mod_loop!
        using Distributions: Poisson, Binomial

        # Test Poisson random numbers
        @test_opt rand(Poisson(10.0))

        # Test Binomial random numbers
        @test_opt rand(Binomial(100, 0.1))
    end

end
