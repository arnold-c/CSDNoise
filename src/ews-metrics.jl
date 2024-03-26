using Statistics
using TestItems

function calculate_centered_mean(timeseries, time_step, bandwidth)
    mean_vec = zeros(Float64, length(timeseries))

    calculate_centered_mean!(mean_vec, timeseries, time_step, bandwidth)

    return mean_vec
end

function calculate_centered_mean!(mean_vec, timeseries, time_step, bandwidth)
    # TODO: implement
    delta = (bandwidth - 1) * time_step
    tlength = length(timeseries)
    for i in eachindex(timeseries)
        if i <= delta
            mean_vec[i] = mean(@view(timeseries[begin:(2i - 1)]))
            continue
        end
        if i > tlength - delta
            mean_vec[i] = mean(@view(timeseries[(2i - tlength):end]))
            continue
        end
        mean_vec[i] = mean(
            @view(timeseries[(i - delta):(i + delta)])
        )
    end
    return mean_vec
end

@testitem "Centered mean" begin
    using CSDNoise
    using Statistics

    daily_testpositives = collect(1:10)

    centered_mean_testpositives = calculate_centered_mean(
        daily_testpositives, 1, 3
    )

    @test isequal(
        centered_mean_testpositives,
        [
            mean([1]),
            mean([1, 2, 3]),
            mean([1, 2, 3, 4, 5]),
            mean([2, 3, 4, 5, 6]),
            mean([3, 4, 5, 6, 7]),
            mean([4, 5, 6, 7, 8]),
            mean([5, 6, 7, 8, 9]),
            mean([6, 7, 8, 9, 10]),
            mean([8, 9, 10]),
            mean([10]),
        ],
    )
end

function calculate_backward_mean(timeseries, time_step, bandwidth)
    mean_vec = zeros(Float64, length(timeseries))

    calculate_backward_mean!(mean_vec, timeseries, time_step, bandwidth)

    return mean_vec
end

function calculate_backward_mean!(mean_vec, timeseries, time_step, bandwidth)
    calculate_movingavg!(mean_vec, timeseries, bandwidth * time_step)
    return mean_vec
end

@testitem "Backward mean" begin
    using CSDNoise
    using Statistics

    daily_testpositives = collect(1:10)

    backward_mean_testpositives = calculate_backward_mean(
        daily_testpositives, 1, 3
    )

    @test isequal(
        backward_mean_testpositives,
        [
            mean([1]),
            mean([1, 2]),
            mean([1, 2, 3]),
            mean([2, 3, 4]),
            mean([3, 4, 5]),
            mean([4, 5, 6]),
            mean([5, 6, 7]),
            mean([6, 7, 8]),
            mean([7, 8, 9]),
            mean([8, 9, 10]),
        ],
    )
end

function calculate_autocorrelation()
end

function calculate_coefficient_of_variation()
end

function calculate_index_of_dispersion()
end

function calculate_kurtosis()
end

function calculate_skewness()
end

function calculate_variance()
end
