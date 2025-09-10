using CSDNoise
using Test
using Dates

@testset "EWS Multistart Optimization" begin
    @testset "create_optimization_scenarios" begin
        # Create minimal test specification vectors
        test_specification_vecs = (
            noise_specification_vec = [
                NoiseSpecification(PoissonNoise(1.0)),
                NoiseSpecification(PoissonNoise(2.0)),
            ],
            test_specification_vec = [
                IndividualTestSpecification(0.9, 0.9, 0),
                IndividualTestSpecification(0.95, 0.95, 0),
            ],
            percent_tested_vec = [1.0],
            ews_metric_specification_vec = [
                EWSMetricSpecification(Backward, Day(28), Week(52), 1),
            ],
            ews_enddate_type_vec = [Reff_start],
            ews_threshold_window_vec = [ExpandingThresholdWindow],
            ews_threshold_burnin_vec = [Year(5)],
            ews_metric_vec = ["autocovariance", "variance"],
        )

        # Test the function
        scenarios = create_optimization_scenarios(test_specification_vecs)

        # Calculate expected number of scenarios
        expected_count = 2 * 2 * 1 * 1 * 1 * 1 * 1 * 2  # 8 scenarios

        @test length(scenarios) == expected_count
        @test scenarios isa Vector

        # Test that all scenarios are NamedTuples with correct fields
        for scenario in scenarios
            @test scenario isa NamedTuple
            @test haskey(scenario, :noise_specification)
            @test haskey(scenario, :test_specification)
            @test haskey(scenario, :percent_tested)
            @test haskey(scenario, :ews_metric_specification)
            @test haskey(scenario, :ews_enddate_type)
            @test haskey(scenario, :ews_threshold_window)
            @test haskey(scenario, :ews_threshold_burnin)
            @test haskey(scenario, :ews_metric)
        end

        # Test that we get all combinations
        noise_specs = [s.noise_specification for s in scenarios]
        test_specs = [s.test_specification for s in scenarios]
        metrics = [s.ews_metric for s in scenarios]

        # Should have both noise specifications
        @test NoiseSpecification(PoissonNoise(1.0)) in noise_specs
        @test NoiseSpecification(PoissonNoise(2.0)) in noise_specs

        # Should have both test specifications
        @test IndividualTestSpecification(0.9, 0.9, 0) in test_specs
        @test IndividualTestSpecification(0.95, 0.95, 0) in test_specs

        # Should have both metrics
        @test "autocovariance" in metrics
        @test "variance" in metrics

        # Test with single values to ensure it still works
        single_spec_vecs = (
            noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))],
            test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
            percent_tested_vec = [1.0],
            ews_metric_specification_vec = [EWSMetricSpecification(Backward, Day(28), Week(52), 1)],
            ews_enddate_type_vec = [Reff_start],
            ews_threshold_window_vec = [ExpandingThresholdWindow],
            ews_threshold_burnin_vec = [Year(5)],
            ews_metric_vec = ["autocovariance"],
        )

        single_scenarios = create_optimization_scenarios(single_spec_vecs)
        @test length(single_scenarios) == 1
        @test single_scenarios[1].noise_specification == NoiseSpecification(PoissonNoise(1.0))
        @test single_scenarios[1].ews_metric == "autocovariance"

        # Test with empty vectors (edge case)
        empty_spec_vecs = (
            noise_specification_vec = PoissonNoiseSpecification[],
            test_specification_vec = IndividualTestSpecification[],
            percent_tested_vec = Float64[],
            ews_metric_specification_vec = EWSMetricSpecification[],
            ews_enddate_type_vec = EWSEndDateType[],
            ews_threshold_window_vec = Union{Type{ExpandingThresholdWindow}, Type{RollingThresholdWindow}}[],
            ews_threshold_burnin_vec = Union{Dates.Day, Dates.Year}[],
            ews_metric_vec = String[],
        )

        empty_scenarios = create_optimization_scenarios(empty_spec_vecs)
        @test length(empty_scenarios) == 0
        @test empty_scenarios isa Vector

        # Test validation - missing required field
        incomplete_spec_vecs = (
            noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))],
            test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
            # Missing other required fields
        )
        @test_throws ArgumentError create_optimization_scenarios(incomplete_spec_vecs)

        # Test validation - extra field
        extra_field_spec_vecs = (
            noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))],
            test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
            percent_tested_vec = [1.0],
            ews_metric_specification_vec = [EWSMetricSpecification(Backward, Day(28), Week(52), 1)],
            ews_enddate_type_vec = [Reff_start],
            ews_threshold_window_vec = [ExpandingThresholdWindow],
            ews_threshold_burnin_vec = [Year(5)],
            ews_metric_vec = ["autocovariance"],
            extra_field = ["should_not_be_here"],  # This should cause an error
        )
        @test_throws ArgumentError create_optimization_scenarios(extra_field_spec_vecs)

        # Test validation - invalid metric name
        invalid_metric_spec_vecs = (
            noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))],
            test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
            percent_tested_vec = [1.0],
            ews_metric_specification_vec = [EWSMetricSpecification(Backward, Day(28), Week(52), 1)],
            ews_enddate_type_vec = [Reff_start],
            ews_threshold_window_vec = [ExpandingThresholdWindow],
            ews_threshold_burnin_vec = [Year(5)],
            ews_metric_vec = ["invalid_metric"],  # This should cause an error
        )
        @test_throws ArgumentError create_optimization_scenarios(invalid_metric_spec_vecs)

        # Test validation - invalid percent tested
        invalid_percent_spec_vecs = (
            noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))],
            test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
            percent_tested_vec = [1.5],  # Invalid: > 1.0
            ews_metric_specification_vec = [EWSMetricSpecification(Backward, Day(28), Week(52), 1)],
            ews_enddate_type_vec = [Reff_start],
            ews_threshold_window_vec = [ExpandingThresholdWindow],
            ews_threshold_burnin_vec = [Year(5)],
            ews_metric_vec = ["autocovariance"],
        )
        @test_throws ArgumentError create_optimization_scenarios(invalid_percent_spec_vecs)
    end

    @testset "_validate_specification_vectors" begin
        # Valid specification vectors for baseline
        valid_spec_vecs = (
            noise_specification_vec = [NoiseSpecification(PoissonNoise(1.0))],
            test_specification_vec = [IndividualTestSpecification(0.9, 0.9, 0)],
            percent_tested_vec = [1.0],
            ews_metric_specification_vec = [EWSMetricSpecification(Backward, Day(28), Week(52), 1)],
            ews_enddate_type_vec = [Reff_start],
            ews_threshold_window_vec = [ExpandingThresholdWindow],
            ews_threshold_burnin_vec = [Year(5)],
            ews_metric_vec = ["autocovariance"],
        )

        # Test that valid specs pass validation
        @test CSDNoise._validate_specification_vectors(valid_spec_vecs) === nothing

        # Test missing required fields
        @testset "Missing required fields" begin
            for field_to_remove in keys(valid_spec_vecs)
                incomplete_specs = NamedTuple(
                    k => v for (k, v) in pairs(valid_spec_vecs) if k != field_to_remove
                )
                @test_throws ArgumentError CSDNoise._validate_specification_vectors(incomplete_specs)
            end
        end

        # Test extra fields
        @testset "Extra fields" begin
            extra_field_specs = merge(valid_spec_vecs, (extra_field = ["not_allowed"],))
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(extra_field_specs)

            multiple_extra_specs = merge(
                valid_spec_vecs, (
                    extra_field1 = ["not_allowed1"],
                    extra_field2 = ["not_allowed2"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(multiple_extra_specs)
        end

        # Test type validation
        @testset "Type validation" begin
            # Invalid noise specification type
            invalid_noise_specs = merge(
                valid_spec_vecs, (
                    noise_specification_vec = ["not_a_noise_spec"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(invalid_noise_specs)

            # Invalid test specification type
            invalid_test_specs = merge(
                valid_spec_vecs, (
                    test_specification_vec = ["not_a_test_spec"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(invalid_test_specs)

            # Invalid EWS metric specification type
            invalid_ews_metric_specs = merge(
                valid_spec_vecs, (
                    ews_metric_specification_vec = ["not_an_ews_metric_spec"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(invalid_ews_metric_specs)

            # Invalid EWS end date type
            invalid_enddate_specs = merge(
                valid_spec_vecs, (
                    ews_enddate_type_vec = ["not_an_enddate_type"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(invalid_enddate_specs)

            # Invalid threshold window type
            invalid_window_specs = merge(
                valid_spec_vecs, (
                    ews_threshold_window_vec = ["not_a_window_type"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(invalid_window_specs)

            # Invalid burnin type
            invalid_burnin_specs = merge(
                valid_spec_vecs, (
                    ews_threshold_burnin_vec = ["not_a_date_period"],
                )
            )
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(invalid_burnin_specs)
        end

        # Test value validation
        @testset "Value validation" begin
            # Invalid percent tested values
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        percent_tested_vec = [-0.1],  # Below 0
                    )
                )
            )

            @test_throws ArgumentError CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        percent_tested_vec = [1.1],   # Above 1
                    )
                )
            )

            @test_throws ArgumentError CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        percent_tested_vec = [0.5, 1.5],  # Mixed valid/invalid
                    )
                )
            )

            # Valid boundary values should pass
            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        percent_tested_vec = [0.0],   # Boundary: 0
                    )
                )
            ) === nothing

            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        percent_tested_vec = [1.0],   # Boundary: 1
                    )
                )
            ) === nothing

            # Invalid EWS metric names
            @test_throws ArgumentError CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_metric_vec = ["invalid_metric"],
                    )
                )
            )

            @test_throws ArgumentError CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_metric_vec = ["autocovariance", "invalid_metric"],  # Mixed valid/invalid
                    )
                )
            )

            # Valid EWS metrics should pass
            valid_metrics = [
                "autocorrelation", "autocovariance", "coefficient_of_variation",
                "index_of_dispersion", "kurtosis", "mean", "skewness", "variance",
            ]
            for metric in valid_metrics
                @test CSDNoise._validate_specification_vectors(
                    merge(
                        valid_spec_vecs, (
                            ews_metric_vec = [metric],
                        )
                    )
                ) === nothing
            end
        end

        # Test valid threshold window types
        @testset "Threshold window types" begin
            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_threshold_window_vec = [ExpandingThresholdWindow],
                    )
                )
            ) === nothing

            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_threshold_window_vec = [RollingThresholdWindow],
                    )
                )
            ) === nothing

            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_threshold_window_vec = [ExpandingThresholdWindow, RollingThresholdWindow],
                    )
                )
            ) === nothing
        end

        # Test valid burnin period types
        @testset "Burnin period types" begin
            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_threshold_burnin_vec = [Day(30)],
                    )
                )
            ) === nothing

            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_threshold_burnin_vec = [Year(1)],
                    )
                )
            ) === nothing

            @test CSDNoise._validate_specification_vectors(
                merge(
                    valid_spec_vecs, (
                        ews_threshold_burnin_vec = [Day(30), Year(1)],
                    )
                )
            ) === nothing
        end

        # Test empty vectors (should warn but not error)
        @testset "Empty vectors" begin
            empty_specs = (
                noise_specification_vec = NoiseSpecification[],
                test_specification_vec = IndividualTestSpecification[],
                percent_tested_vec = Float64[],
                ews_metric_specification_vec = EWSMetricSpecification[],
                ews_enddate_type_vec = EWSEndDateType[],
                ews_threshold_window_vec = Union{Type{ExpandingThresholdWindow}, Type{RollingThresholdWindow}}[],
                ews_threshold_burnin_vec = Union{Dates.Day, Dates.Year}[],
                ews_metric_vec = String[],
            )

            # Should not throw an error, but should warn
            @test CSDNoise._validate_specification_vectors(empty_specs) === nothing
        end

        # Test mixed valid types
        @testset "Mixed valid types" begin
            mixed_specs = (
                noise_specification_vec = [
                    NoiseSpecification(PoissonNoise(1.0)),
                    NoiseSpecification(PoissonNoise(2.0)),
                    NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.8734)),
                ],
                test_specification_vec = [
                    IndividualTestSpecification(0.9, 0.9, 0),
                    IndividualTestSpecification(0.95, 0.95, 0),
                ],
                percent_tested_vec = [0.5, 1.0],
                ews_metric_specification_vec = [
                    EWSMetricSpecification(Backward, Day(28), Week(52), 1),
                    EWSMetricSpecification(Centered, Day(14), Week(26), 2),
                ],
                ews_enddate_type_vec = [Reff_start, Reff_end],
                ews_threshold_window_vec = [ExpandingThresholdWindow, RollingThresholdWindow],
                ews_threshold_burnin_vec = [Day(30), Year(1)],
                ews_metric_vec = ["autocovariance", "variance", "mean"],
            )

            @test CSDNoise._validate_specification_vectors(mixed_specs) === nothing
        end
    end
end
