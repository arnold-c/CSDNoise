@testset "ews-metrics.jl" begin
    using CSDNoise, StatsBase

    @testset "Timeseries aggregation" begin
        testvec = collect(1:10)
        @test isequal(aggregate_timeseries(testvec, Day(1)), testvec)
        @test isequal(
            aggregate_timeseries(testvec, Day(2)),
            [sum(1:2), sum(3:4), sum(5:6), sum(7:8), sum(9:10)],
        )
        @test isequal(
            aggregate_timeseries(testvec, Day(3)),
            [sum(1:3), sum(4:6), sum(7:9)],
        )
        @test isequal(
            aggregate_timeseries(testvec, Day(3), StatsBase.mean),
            [StatsBase.mean(1:3), StatsBase.mean(4:6), StatsBase.mean(7:9)],
        )

        @test isequal(
            aggregate_timeseries(testvec, Week(1)), [sum(1:7)]
        )
        @test isequal(
            aggregate_timeseries(testvec, Week(1), StatsBase.mean),
            [StatsBase.mean(1:7)],
        )
    end

    @testset "Mean functions" begin
        testvec = collect(1:10)
        mean_vec = spaero_mean(
            Backward, testvec, Day(2)
        )
    end
end
