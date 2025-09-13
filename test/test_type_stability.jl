using Test
using JET
using DrWatson
@quickactivate("CSDNoise")
using CSDNoise
using StaticArrays
using Dates
using NLopt
using StructArrays

# Helper function to create test specification vectors
function create_test_specification_vectors()
    noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))]
    test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)]
    percent_tested_vec = [1.0]

    ews_method_vec = [EWSMethod(Backward())]
    ews_aggregation_vec = [Day(28)]
    ews_bandwidth_vec = [Week(52)]
    ews_lag_days_vec = [1]

    ews_metric_specification_vec = create_combinations_vec(
        EWSMetricSpecification,
        (ews_method_vec, ews_aggregation_vec, ews_bandwidth_vec, ews_lag_days_vec)
    )

    ews_metric_vec = ["autocovariance"]
    ews_enddate_type_vec = [EWSEndDateType(Reff_start())]
    ews_threshold_window_vec = [EWSThresholdWindowType(ExpandingThresholdWindow())]
    ews_threshold_percentile_vec = collect(0.5:0.02:0.99)
    ews_consecutive_thresholds_vec = collect(2:2:30)
    ews_threshold_burnin_vec = [Year(5)]

    return (;
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec,
    )
end

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
                seasonality = SeasonalityFunction(CosineSeasonality())
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

            # Test create_inc_infec_arr which is used in generate_ensemble_data
            ensemble_data = generate_single_ensemble(ensemble_spec; seed = 42)
            @test_opt target_modules = (CSDNoise,) create_inc_infec_arr(
                ensemble_data[:ensemble_inc_vecs], outbreak_spec
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

    @testset "Multistart Optimization Functions" begin
        # Create test data for optimization functions
        specification_vecs = create_test_specification_vectors()
        data_arrs = generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)

        @testset "High Priority Optimization Functions" begin
            # Create optimization scenario
            scenarios = create_optimization_scenarios(specification_vecs)
            scenario = scenarios[1]

            @test_opt target_modules = (CSDNoise,) create_cached_simulation_data(
                scenario, data_arrs
            )

            cached_data = create_cached_simulation_data(scenario, data_arrs)
            tracker = OptimizationTracker()
            test_params = [0.9, 5.0]

            @test_opt target_modules = (CSDNoise,) map_continuous_to_ews_parameters(test_params)

            @test_opt target_modules = (CSDNoise,) ews_objective_function_with_tracking(
                test_params, scenario, cached_data, tracker
            )

            @test_opt target_modules = (CSDNoise,) CSDNoise.calculate_ews_metrics_for_simulation(
                scenario.ews_metric_specification,
                cached_data.testarr[:, :, 1],
                cached_data.null_testarr[:, :, 1],
                cached_data.thresholds[1],
                scenario.ews_enddate_type
            )

        end

        @testset "Medium Priority Scenario Management Functions" begin

            @test_opt target_modules = (CSDNoise,) create_optimization_scenarios(
                specification_vecs
            )

            # Test optimize_single_scenario with minimal configuration
            scenarios = create_optimization_scenarios(specification_vecs)
            scenario = scenarios[1]
            cached_data = create_cached_simulation_data(scenario, data_arrs)

            bounds = (;
                lowers = [0.5, 2.0],
                uppers = [0.99, 30.0],
            )

            config = (;
                local_algorithm = NLopt.LN_BOBYQA,
                xtol_rel = 1.0e-3,
                xtol_abs = 1.0e-3,
                maxeval = 100,
                n_sobol_points = 10,
            )


            @test_opt target_modules = (CSDNoise,) optimize_single_scenario(
                scenario,
                data_arrs,
                bounds,
                config
            )

            # Test dataframe_row_to_scenario
            scenarios_df = CSDNoise.create_scenarios_dataframe(specification_vecs)
            test_row = scenarios_df[1, :]
            @test_opt target_modules = (CSDNoise,) CSDNoise.dataframe_row_to_scenario(test_row)
        end

        @testset "Lower Priority Data Management Functions" begin
            @test_opt target_modules = (CSDNoise,) CSDNoise.create_scenarios_dataframe(
                specification_vecs
            )

            # Test create_results_dataframe
            scenarios = create_optimization_scenarios(specification_vecs)
            dummy_results = StructVector(
                [
                    OptimizedValues(
                        0.9,
                        5,
                        0.85,
                        0.8,
                        0.9,
                    ),
                ]
            )

            @test_opt target_modules = (CSDNoise,) CSDNoise.create_results_dataframe(
                scenarios,
                dummy_results,
            )

            # Test scenarios_equal
            scenario1 = scenarios[1]
            scenario2 = scenarios[1]
            @test_opt target_modules = (CSDNoise,) CSDNoise.scenarios_equal(scenario1, scenario2)
        end
    end

    @testset "Helper Functions for create_ensemble_specs" begin
        @testset "Parameter calculation functions" begin
            # Test calculate_mu
            @test_opt target_modules = (CSDNoise,) calculate_mu(27)

            # Test calculate_beta
            @test_opt target_modules = (CSDNoise,) calculate_beta(
                R0, GAMMA, calculate_mu(27), 1, 500_000
            )

            # Test calculate_import_rate
            @test_opt target_modules = (CSDNoise,) calculate_import_rate(
                calculate_mu(27), R0, 500_000
            )

            # Test calculate_vaccination_rate_to_achieve_Reff
            @test_opt target_modules = (CSDNoise,) calculate_vaccination_rate_to_achieve_Reff(
                0.9, 10, 25_000, 500_000, R0, calculate_mu(27)
            )
        end

        @testset "Specification constructors" begin
            # Test SimTimeParameters
            @test_opt target_modules = (CSDNoise,) SimTimeParameters(;
                burnin = 365.0 * 5, tmin = 0.0, tmax = 365.0 * 20, tstep = 1.0
            )

            # Test StateParameters
            @test_opt target_modules = (CSDNoise,) StateParameters(
                500_000,
                Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
            )

            # Test DynamicsParameterSpecification
            mu = calculate_mu(27)
            beta_mean = calculate_beta(R0, GAMMA, mu, 1, 500_000)
            epsilon = calculate_import_rate(mu, R0, 500_000)

            @test_opt target_modules = (CSDNoise,) DynamicsParameterSpecification(
                beta_mean, 0.0, SeasonalityFunction(CosineSeasonality()), SIGMA, GAMMA, mu, 27.0, epsilon, R0,
                0.6, 1.0, 0.6, 0.8
            )

            # Test EnsembleSpecification
            time_spec = SimTimeParameters(; burnin = 365.0 * 5, tmin = 0.0, tmax = 365.0 * 20, tstep = 1.0)
            state_spec = StateParameters(500_000, Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95))
            dynamics_spec_test = DynamicsParameterSpecification(
                beta_mean, 0.0, SeasonalityFunction(CosineSeasonality()), SIGMA, GAMMA, mu, 27.0, epsilon, R0,
                0.6, 1.0, 0.6, 0.8
            )

            @test_opt target_modules = (CSDNoise,) EnsembleSpecification(
                state_spec, dynamics_spec_test, time_spec, 3
            )

            # Test OutbreakSpecification
            @test_opt target_modules = (CSDNoise,) OutbreakSpecification(5, 30, 500)
        end
    end

    @testset "Benchmark Functions" begin
        @testset "Ensemble and scenario calculation functions" begin
            # Test create_ensemble_specs
            @test_opt target_modules = (CSDNoise,) create_ensemble_specs(100)

            # Test calculate_scenarios
            specification_vecs = create_test_specification_vectors()
            @test_opt target_modules = (CSDNoise,) calculate_scenarios(specification_vecs)

            # Test calculate_grid_points
            @test_opt target_modules = (CSDNoise,) calculate_grid_points(specification_vecs)
        end
    end

    @testset "EWS Metrics Functions" begin
        # Test EWSMetrics constructor with various specifications
        example_ews_spec = EWSMetricSpecification(EWSMethod(Backward()), Dates.Day(1), Dates.Day(30), 1)
        example_timeseries = collect(1.0:35.0)

        @test_opt target_modules = (CSDNoise,) EWSMetrics(example_ews_spec, example_timeseries)

        # Test expanding_ews_thresholds function
        example_ews_metrics = EWSMetrics(example_ews_spec, example_timeseries)

        @test_opt target_modules = (CSDNoise,) expanding_ews_thresholds(
            example_ews_metrics,
            :mean,
            EWSThresholdWindowType(ExpandingThresholdWindow()),
            0.95,
            Dates.Day(5)
        )

    end
end
