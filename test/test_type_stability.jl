using Test
using JET
using DrWatson
@quickactivate("CSDNoise")
using CSDNoise
using StaticArrays

@testset "Type Stability Tests" begin
    @testset "generate_ensemble_data type stability" begin
        # Create small test data to minimize test time
        ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(5)

        # Test that generate_ensemble_data has no runtime dispatches
        @test_opt target_modules = (CSDNoise,) generate_ensemble_data(
            ensemble_spec, null_spec, outbreak_spec
        )
    end

    @testset "generate_single_ensemble type stability" begin
        # Create small test data to minimize test time
        ensemble_spec, _, _ = create_ensemble_specs(3)  # Even smaller for focused testing

        # Test that generate_single_ensemble has no runtime dispatches
        @test_opt target_modules = (CSDNoise,) generate_single_ensemble(
            ensemble_spec; seed = 42
        )
    end

    @testset "generate_single_ensemble component functions type stability" begin
        # Create test data for component testing
        ensemble_spec, _, _ = create_ensemble_specs(2)

        # Test array allocations used in generate_single_ensemble
        @testset "Array allocations" begin
            @test_opt Vector{DynamicsParameters}(undef, ensemble_spec.nsims)
        end

        @testset "Core simulation loop components" begin
            # Test the main simulation components
            state_parameters = ensemble_spec.state_parameters
            dynamics_spec = ensemble_spec.dynamics_parameter_specification
            time_parameters = ensemble_spec.time_parameters

            # Test DynamicsParameters creation with different seeds
            @test_opt target_modules = (CSDNoise,) DynamicsParameters(
                dynamics_spec; seed = 42
            )

            # Test SVector conversion from init_states
            @test_opt SVector(state_parameters.init_states)

            # Test seir_mod! with views (as used in the ensemble loop)
            dynamics_params = DynamicsParameters(dynamics_spec; seed = 42)
            tlength = time_parameters.tlength

            # Create test arrays as in generate_single_ensemble
            init_state_type = typeof(state_parameters.init_states)
            ensemble_seir_vecs = Array{SVector{5, Int64}, 2}(undef, tlength, 2)
            ensemble_inc_vecs = Array{SVector{1, Int64}, 2}(undef, tlength, 2)
            ensemble_beta_arr = zeros(Float64, tlength)

            @test_opt target_modules = (CSDNoise,) seir_mod!(
                @view(ensemble_seir_vecs[:, 1]),
                ensemble_inc_vecs[:, 1],
                ensemble_beta_arr,
                SVector(state_parameters.init_states),
                dynamics_params,
                time_parameters,
                42
            )

        end

        @testset "Post-processing functions type stability" begin
            # Test convert_svec_to_array (used after simulation)
            tlength = ensemble_spec.time_parameters.tlength
            nsims = 2
            init_state_type = typeof(ensemble_spec.state_parameters.init_states)

            # Create test data similar to what generate_single_ensemble produces
            test_svec_array = Array{init_state_type, 2}(undef, tlength, nsims)
            # Fill with dummy data
            for i in 1:tlength, j in 1:nsims
                test_svec_array[i, j] = ensemble_spec.state_parameters.init_states
            end

            @test_opt target_modules = (CSDNoise,) convert_svec_to_array(test_svec_array)

            # Test calculateReffective_t! function
            dynamics_params = DynamicsParameters(
                ensemble_spec.dynamics_parameter_specification; seed = 42
            )
            ensemble_beta_arr = zeros(Float64, tlength)
            ensemble_Reff_arr = zeros(Float64, tlength, nsims)

            # Create dummy SEIR array for testing
            ensemble_seir_arr = zeros(Int64, tlength, 5, nsims)

            @test_opt target_modules = (CSDNoise,) calculateReffective_t!(
                @view(ensemble_Reff_arr[:, 1]),
                ensemble_beta_arr,
                dynamics_params,
                1,
                @view(ensemble_seir_arr[:, :, 1])
            )

            # Test Reff_ge_than_one function
            test_Reff_vec = rand(Float64, tlength)
            @test_opt target_modules = (CSDNoise,) Reff_ge_than_one(test_Reff_vec)
        end

        @testset "Return value construction type stability" begin
            # Test the named tuple return construction
            ensemble_spec, _, _ = create_ensemble_specs(2)
            tlength = ensemble_spec.time_parameters.tlength
            nsims = ensemble_spec.nsims

            # Create dummy data matching the return structure
            ensemble_seir_arr = zeros(Int64, tlength, 5, nsims)
            dynamics_parameters = [
                DynamicsParameters(ensemble_spec.dynamics_parameter_specification; seed = 42 + i)
                    for i in 1:nsims
            ]
            ensemble_Reff_arr = zeros(Float64, tlength, nsims)
            ensemble_Reff_thresholds_vec = [zeros(Int64, 10, 3) for _ in 1:nsims]
            ensemble_inc_vecs = Array{SVector{1, Int64}, 2}(undef, tlength, nsims)

            # Fill inc_vecs with dummy data
            for i in 1:tlength, j in 1:nsims
                ensemble_inc_vecs[i, j] = SVector(1)
            end

            # Test named tuple construction (as done in return statement)
            @test_opt (
                ensemble_seir_arr = ensemble_seir_arr,
                ensemble_spec = ensemble_spec,
                dynamics_parameters = dynamics_parameters,
                ensemble_Reff_arr = ensemble_Reff_arr,
                ensemble_Reff_thresholds_vec = ensemble_Reff_thresholds_vec,
                ensemble_inc_vecs = ensemble_inc_vecs,
            )
        end
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
