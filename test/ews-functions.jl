@testset "ews-functions.jl" begin
    @testset "EWS Thresholds" begin
        using CSDNoise, StatsBase

        @testset "No NaNs" begin
            testvec = collect(1:20)
            percentiles = [0.5, 0.9]
            burn_in = 5

            res_percentiles, exceeds_thresholds = expanding_ews_thresholds(
                testvec,
                Expanding;
                percentiles = percentiles,
                burn_in = burn_in,
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
            expected_percentiles[1:burn_in, :] .= NaN

            @test isequal(sum(isnan.(res_percentiles)), 10)
            @test isequal(sum(exceeds_thresholds), 28)
            @test isequal(res_percentiles, expected_percentiles)
        end

        @testset "NaNs" begin
            testvec = [NaN, NaN, collect(1:10)..., NaN, collect(11:20)...]
            percentiles = 0.5
            burn_in = 5

            res_percentiles, exceeds_thresholds, worker_vec = expanding_ews_thresholds(
                testvec,
                Expanding;
                percentiles = percentiles,
                burn_in = burn_in,
            )

            filtered_testvec = filter(!isnan, testvec)
            expected_percentiles = map(
                i -> StatsBase.quantile(filtered_testvec[1:i], percentiles[1]),
                eachindex(filtered_testvec),
            )
            expected_percentiles[1:(burn_in - 2)] .= NaN
            expected_percentiles = vcat(
                NaN,
                NaN,
                expected_percentiles[1:10],
                expected_percentiles[10],
                expected_percentiles[11:20],
            )

            expected_exceeds_thresholds = vcat(
                repeat([false], burn_in + 1),
                repeat([true], 6),
                false,
                repeat([true], 10),
            )

            @test isequal(sum(isnan.(res_percentiles)), 5)
            @test isequal(vec(exceeds_thresholds), expected_exceeds_thresholds)
            @test isequal(vec(res_percentiles), expected_percentiles)
        end
    end
end
