using Test
using CSDNoise
using StatsBase

@testset "diag-testing-functions.jl" begin

    @testset "Moving average" begin

        daily_testpositives = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        movingavg_testpositives = calculate_movingavg(
            daily_testpositives,
            5,
        )

        @test isequal(
            movingavg_testpositives,
            [
                mean([1]),
                mean([1, 2]),
                mean([1, 2, 3]),
                mean([1, 2, 3, 4]),
                mean([1, 2, 3, 4, 5]),
                mean([2, 3, 4, 5, 6]),
                mean([3, 4, 5, 6, 7]),
                mean([4, 5, 6, 7, 8]),
                mean([5, 6, 7, 8, 9]),
                mean([6, 7, 8, 9, 10]),
            ],
        )

        @test length(movingavg_testpositives) == length(daily_testpositives)

        daily_testpositives = zeros(Int64, 10, 2)
        daily_testpositives[:, 1] .= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        calculate_movingavg!(
            @view(daily_testpositives[:, 2]),
            @view(daily_testpositives[:, 1]),
            5,
        )

        @test isequal(
            daily_testpositives[:, 2],
            Int64.(
                round.(
                    [
                        mean([1]),
                        mean([1, 2]),
                        mean([1, 2, 3]),
                        mean([1, 2, 3, 4]),
                        mean([1, 2, 3, 4, 5]),
                        mean([2, 3, 4, 5, 6]),
                        mean([3, 4, 5, 6, 7]),
                        mean([4, 5, 6, 7, 8]),
                        mean([5, 6, 7, 8, 9]),
                        mean([6, 7, 8, 9, 10]),
                    ]
                )
            ),
        )
    end

    @testset "Inferred cases" begin
        using CSDNoise, StatsBase

        daily_testpositives = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        daily_tested = repeat([10], length(daily_testpositives))

        movingavg_testpositives = calculate_test_positivity_rate(
            daily_testpositives,
            daily_tested,
            5,
        )

        @test isequal(
            movingavg_testpositives,
            [
                mean([1]) / 10,
                mean([1, 2]) / 10,
                mean([1, 2, 3]) / 10,
                mean([1, 2, 3, 4]) / 10,
                mean([1, 2, 3, 4, 5]) / 10,
                mean([2, 3, 4, 5, 6]) / 10,
                mean([3, 4, 5, 6, 7]) / 10,
                mean([4, 5, 6, 7, 8]) / 10,
                mean([5, 6, 7, 8, 9]) / 10,
                mean([6, 7, 8, 9, 10]) / 10,
            ],
        )

        @test length(movingavg_testpositives) == length(daily_testpositives)

        daily_clinic_visits = repeat([100], length(daily_testpositives))

        inferred_positives_no_lag = infer_true_positives(
            movingavg_testpositives,
            daily_testpositives,
            daily_tested,
            daily_clinic_visits,
            0,
        )

        @test isequal(
            inferred_positives_no_lag,
            daily_testpositives +
                movingavg_testpositives .* (daily_clinic_visits - daily_tested),
        )

        inferred_positives_3d_lag = infer_true_positives(
            movingavg_testpositives,
            daily_testpositives,
            daily_tested,
            daily_clinic_visits,
            3,
        )

        @test isequal(
            inferred_positives_3d_lag,
            movingavg_testpositives .* daily_clinic_visits,
        )
    end

end
