using Test
using JET
using DrWatson
@quickactivate("CSDNoise")
using CSDNoise
using StaticArrays

@testset "Type Stability Tests" begin
    # Single setup for all tests - use consistent small size for performance
    ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(3)
    state_params = ensemble_spec.state_parameters
    dynamics_spec = ensemble_spec.dynamics_parameter_specification
    time_params = ensemble_spec.time_parameters

    @testset "SEIR Model Functions" begin
        # Create test data once for this section
        dynamics_params = DynamicsParameters(dynamics_spec; seed = 42)
        initial_states = SVector(state_params.init_states)

        @testset "Main SEIR functions" begin
            @test_opt target_modules = (CSDNoise,) seir_mod(
                initial_states,
                dynamics_params,
                time_params;
                seed = 42
            )

            # Test in-place version
            state_vec = Vector{typeof(initial_states)}(undef, time_params.tlength)
            beta_vec = Vector{Float64}(undef, time_params.tlength)
            inc_vec = Vector{SVector{1, Int64}}(undef, time_params.tlength)

            @test_opt target_modules = (CSDNoise,) seir_mod!(
                state_vec,
                inc_vec,
                beta_vec,
                initial_states,
                dynamics_params,
                time_params,
                42
            )
        end

        @testset "SEIR helper functions" begin
            # Test inner loop function
            @test_opt target_modules = (CSDNoise,) seir_mod_loop!(
                initial_states,
                1.5,  # beta_t
                dynamics_params.mu,
                dynamics_params.epsilon,
                dynamics_params.sigma,
                dynamics_params.gamma,
                dynamics_params.R_0,
                dynamics_params.vaccination_coverage,
                time_params.tstep
            )

            # Test beta calculation
            @test_opt target_modules = (CSDNoise,) calculate_beta_amp(
                dynamics_params.beta_mean,
                dynamics_params.beta_force,
                100.0;
                seasonality = cos
            )
        end

        @testset "Constructor type stability" begin
            @test_opt target_modules = (CSDNoise,) DynamicsParameters(dynamics_spec; seed = 42)
        end
    end

    @testset "Ensemble Generation Functions" begin
        @testset "High-level ensemble functions" begin
            @test_opt target_modules = (CSDNoise,) generate_ensemble_data(
                ensemble_spec, null_spec, outbreak_spec
            )

            @test_opt target_modules = (CSDNoise,) generate_single_ensemble(
                ensemble_spec; seed = 42
            )
        end

        @testset "Ensemble components" begin
            # Test seir_mod! with views as used in ensemble generation
            dynamics_params = DynamicsParameters(dynamics_spec; seed = 42)
            tlength = time_params.tlength

            # Create test arrays as in generate_single_ensemble
            ensemble_seir_vecs = Array{SVector{5, Int64}, 2}(undef, tlength, 2)
            ensemble_inc_vecs = Array{SVector{1, Int64}, 2}(undef, tlength, 2)
            ensemble_beta_arr = zeros(Float64, tlength)

            @test_opt target_modules = (CSDNoise,) seir_mod!(
                @view(ensemble_seir_vecs[:, 1]),
                ensemble_inc_vecs[:, 1],
                ensemble_beta_arr,
                SVector(state_params.init_states),
                dynamics_params,
                time_params,
                42
            )
        end
    end

    @testset "Utility Functions" begin
        @testset "Array conversion functions" begin
            # Test convert_svec_to_array for 3D arrays (main use case)
            tlength = time_params.tlength
            nsims = 2
            test_svec_array = Array{typeof(state_params.init_states), 2}(undef, tlength, nsims)

            # Fill with dummy data
            for i in 1:tlength, j in 1:nsims
                test_svec_array[i, j] = state_params.init_states
            end

            @test_opt target_modules = (CSDNoise,) convert_svec_to_array(test_svec_array)

            # Test convert_svec_to_matrix for 2D arrays
            test_svec_vec = [SVector(1, 2, 3, 4, 5), SVector(2, 3, 4, 5, 6)]
            @test_opt target_modules = (CSDNoise,) convert_svec_to_matrix(test_svec_vec)

            # Test in-place version
            test_matrix = Matrix{Int64}(undef, length(test_svec_vec), length(test_svec_vec[1]))
            @test_opt target_modules = (CSDNoise,) convert_svec_to_matrix!(test_matrix, test_svec_vec)
        end

        @testset "R_eff calculation functions" begin
            dynamics_params = DynamicsParameters(dynamics_spec; seed = 42)
            tlength = time_params.tlength
            nsims = 2

            # Test calculateReffective_t! function
            ensemble_beta_arr = zeros(Float64, tlength)
            ensemble_Reff_arr = zeros(Float64, tlength, nsims)
            ensemble_seir_arr = zeros(Int64, tlength, 5, nsims)

            @test_opt target_modules = (CSDNoise,) calculateReffective_t!(
                @view(ensemble_Reff_arr[:, 1]),
                ensemble_beta_arr,
                dynamics_params,
                1,
                @view(ensemble_seir_arr[:, :, 1])
            )

            # Test Reff threshold detection
            test_Reff_vec = rand(Float64, tlength)
            @test_opt target_modules = (CSDNoise,) Reff_ge_than_one(test_Reff_vec)
        end
    end
end
