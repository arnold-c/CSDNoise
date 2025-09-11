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
    dynamics_params = DynamicsParameters(dynamics_spec; seed = 42)

    @test_opt target_modules = (CSDNoise,) seir_mod(
        SVector(state_params.init_states),
        dynamics_params,
        time_params;
        seed = 42
    )

end
