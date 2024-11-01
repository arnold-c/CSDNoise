@testset "ews-functions.jl" begin
    @testset "EWS Thresholds" begin
        using CSDNoise, StatsBase, Dates

        @testset "No NaNs" begin
            testvec = Float64.(collect(1:20))
            percentiles = [0.5, 0.9]
            daily_burn_in_int = 5
            daily_burn_in = Day(daily_burn_in_int)

            test_ewsmetrics = EWSMetrics(
                EWSMetricSpecification(
                    Backward,
                    Day(1),
                    Day(1),
                    1,
                ),
                repeat([testvec], 8)...,
                repeat([1.0], 8)...,
            )

            outer_res_percentiles, outer_exceeds_thresholds = expanding_ews_thresholds(
                test_ewsmetrics,
                :mean,
                ExpandingThresholdWindow;
                percentiles = percentiles,
                burn_in = daily_burn_in,
            )[1:2]

            res_percentiles, exceeds_thresholds = CSDNoise._expanding_ews_thresholds(
                testvec,
                ExpandingThresholdWindow,
                percentiles,
                daily_burn_in_int,
            )[1:2]

            expected_percentiles = map(
                i -> StatsBase.quantile(testvec[1:i], percentiles[1]),
                eachindex(testvec),
            )
            expected_percentiles = hcat(
                expected_percentiles,
                map(
                    i -> StatsBase.quantile(testvec[1:i], percentiles[2]),
                    eachindex(testvec),
                ),
            )
            expected_percentiles[1:daily_burn_in_int, :] .= NaN

            @test isequal(sum(isnan.(res_percentiles)), 10)
            @test isequal(sum(exceeds_thresholds), 28)
            @test isequal(res_percentiles, expected_percentiles)

            @test isequal(outer_res_percentiles, res_percentiles)
            @test isequal(outer_exceeds_thresholds, exceeds_thresholds)

            weekly_burn_in_int = 1
            weekly_burn_in = Week(weekly_burn_in_int)

            weekly_outer_res_percentiles, weekly_outer_exceeds_thresholds = expanding_ews_thresholds(
                test_ewsmetrics,
                :mean,
                ExpandingThresholdWindow;
                percentiles = percentiles,
                burn_in = weekly_burn_in,
            )[1:2]

            weekly_expected_percentiles = map(
                i -> StatsBase.quantile(testvec[1:i], percentiles[1]),
                eachindex(testvec),
            )
            weekly_expected_percentiles = hcat(
                weekly_expected_percentiles,
                map(
                    i -> StatsBase.quantile(testvec[1:i], percentiles[2]),
                    eachindex(testvec),
                ),
            )
            weekly_expected_percentiles[
                1:(Dates.days(weekly_burn_in)), :,
            ] .= NaN

            @test isequal(
                weekly_outer_res_percentiles, weekly_expected_percentiles
            )
        end

        @testset "NaNs" begin
            testvec = [NaN, NaN, collect(1:10)..., NaN, collect(11:20)...]
            percentiles = 0.5
            burn_in_int = 5
            burn_in = Day(burn_in_int)

            test_ewsmetrics = EWSMetrics(
                EWSMetricSpecification(
                    Backward,
                    Day(1),
                    Day(1),
                    1,
                ),
                repeat([testvec], 8)...,
                repeat([1.0], 8)...,
            )

            outer_res_percentiles, outer_exceeds_thresholds = expanding_ews_thresholds(
                test_ewsmetrics,
                :mean,
                ExpandingThresholdWindow;
                percentiles = percentiles,
                burn_in = burn_in,
            )[1:2]

            res_percentiles, exceeds_thresholds, worker_vec = CSDNoise._expanding_ews_thresholds(
                testvec,
                ExpandingThresholdWindow,
                percentiles,
                burn_in_int,
            )

            filtered_testvec = filter(!isnan, testvec)
            expected_percentiles = map(
                i -> StatsBase.quantile(filtered_testvec[1:i], percentiles[1]),
                eachindex(filtered_testvec),
            )
            expected_percentiles[1:(burn_in_int - 2)] .= NaN
            expected_percentiles = vcat(
                NaN,
                NaN,
                expected_percentiles[1:10],
                expected_percentiles[10],
                expected_percentiles[11:20],
            )

            expected_exceeds_thresholds = vcat(
                repeat([false], burn_in_int + 1),
                repeat([true], 6),
                false,
                repeat([true], 10),
            )

            @test isequal(sum(isnan.(res_percentiles)), 5)
            @test isequal(vec(exceeds_thresholds), expected_exceeds_thresholds)
            @test isequal(vec(res_percentiles), expected_percentiles)

            @test isequal(outer_res_percentiles, res_percentiles)
            @test isequal(outer_exceeds_thresholds, exceeds_thresholds)
        end
    end
end
