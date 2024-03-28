using StatsBase: StatsBase
using TestItems

function calculate_ews_metric(metric_function, timeseries, time_step, bandwidth)
    metric_vec = zeros(Float64, length(timeseries))

    metric_function(metric_vec, timeseries, time_step, bandwidth)

    return metric_vec
end

function calculate_centered_mean!(mean_vec, timeseries, time_step, bandwidth)
    calculate_centered_metric!(
        mean_vec, StatsBase.mean, timeseries, time_step, bandwidth
    )
    return mean_vec
end

function calculate_centered_metric!(
    metric_vec, metric_function, timeseries, time_step, bandwidth
)
    delta = (bandwidth - 1) * time_step
    tlength = length(timeseries)
    for i in eachindex(timeseries)
        if i <= delta
            metric_vec[i] = metric_function(@view(timeseries[begin:(2i - 1)]))
            continue
        end
        if i > tlength - delta
            metric_vec[i] = metric_function(
                @view(timeseries[(2i - tlength):end])
            )
            continue
        end
        metric_vec[i] = metric_function(
            @view(timeseries[(i - delta):(i + delta)])
        )
    end
    return metric_vec
end

function calculate_backward_mean!(mean_vec, timeseries, time_step, bandwidth)
    calculate_backward_metric!(
        mean_vec, StatsBase.mean, timeseries, time_step, bandwidth
    )
    return mean_vec
end

function calculate_backward_metric!(
    metric_vec, metric_function, timeseries, time_step, bandwidth
)
    avglag = bandwidth * time_step
    @inbounds for i in eachindex(timeseries)
        @inline moveavg_daystart = calculate_daily_movingavg_startday(i, avglag)
        metric_vec[i] = metric_function(@view(timeseries[moveavg_daystart:i]))
    end
    return nothing
end

@testitem "Mean" begin
    using CSDNoise
    using Statistics

    daily_testpositives = collect(1:10)

    centered_mean_testpositives = calculate_ews_metric(
        calculate_centered_mean!, daily_testpositives, 1, 3
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

    backward_mean_testpositives = calculate_ews_metric(
        calculate_backward_mean!, daily_testpositives, 1, 3
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

function calculate_centered_variance!(
    variance_vec, timeseries, time_step, bandwidth
)
    calculate_centered_metric!(
        variance_vec, StatsBase.var, timeseries, time_step, bandwidth
    )

    return variance_vec
end

function calculate_backward_variance!(
    variance_vec, timeseries, time_step, bandwidth
)
    return calculate_backward_metric!(
        variance_vec, StatsBase.var, timeseries, time_step, bandwidth)
end

@testitem "Variance" begin
    using CSDNoise
    using Statistics

    daily_testpositives = collect(1:10)

    centered_variance_testpositives = calculate_ews_metric(
        calculate_centered_variance!, daily_testpositives, 1, 3
    )

    @test isequal(
        centered_variance_testpositives,
        [
            var([1]),
            var([1, 2, 3]),
            var([1, 2, 3, 4, 5]),
            var([2, 3, 4, 5, 6]),
            var([3, 4, 5, 6, 7]),
            var([4, 5, 6, 7, 8]),
            var([5, 6, 7, 8, 9]),
            var([6, 7, 8, 9, 10]),
            var([8, 9, 10]),
            var([10]),
        ],
    )

    backward_variance_testpositives = calculate_ews_metric(
        calculate_backward_variance!, daily_testpositives, 1, 3
    )

    @test isequal(
        backward_variance_testpositives,
        [
            var([1]),
            var([1, 2]),
            var([1, 2, 3]),
            var([2, 3, 4]),
            var([3, 4, 5]),
            var([4, 5, 6]),
            var([5, 6, 7]),
            var([6, 7, 8]),
            var([7, 8, 9]),
            var([8, 9, 10]),
        ],
    )
end

function calculate_centered_coefficient_of_variation!(
    cov_vec, timeseries, time_step, bandwidth
)
    calculate_centered_metric!(
        cov_vec, StatsBase.variation, timeseries, time_step, bandwidth
    )
    return cov_vec
end

function calculate_backward_coefficient_of_variation!(
    cov_vec, timeseries, time_step, bandwidth
)
    return calculate_backward_metric!(
        cov_vec, StatsBase.variation, timeseries, time_step, bandwidth
    )
end

@testitem "Coefficient of Variation" begin
    using CSDNoise
    using StatsBase

    daily_testpositives = collect(1:10)

    centered_cov_testpostives = calculate_ews_metric(
        calculate_centered_coefficient_of_variation!, daily_testpositives, 1, 3
    )

    @test isequal(
        centered_cov_testpostives,
        [
            variation([1]),
            variation([1, 2, 3]),
            variation([1, 2, 3, 4, 5]),
            variation([2, 3, 4, 5, 6]),
            variation([3, 4, 5, 6, 7]),
            variation([4, 5, 6, 7, 8]),
            variation([5, 6, 7, 8, 9]),
            variation([6, 7, 8, 9, 10]),
            variation([8, 9, 10]),
            variation([10]),
        ],
    )

    backward_cov_testpostives = calculate_ews_metric(
        calculate_backward_coefficient_of_variation!, daily_testpositives, 1, 3
    )

    @test isequal(
        backward_cov_testpostives,
        [
            variation([1]),
            variation([1, 2]),
            variation([1, 2, 3]),
            variation([2, 3, 4]),
            variation([3, 4, 5]),
            variation([4, 5, 6]),
            variation([5, 6, 7]),
            variation([6, 7, 8]),
            variation([7, 8, 9]),
            variation([8, 9, 10]),
        ],
    )
end

function calculate_centered_index_of_dispersion!(
    index_of_dispersion_vec, timeseries, time_step, bandwidth
)
    calculate_centered_metric!(
        index_of_dispersion_vec,
        iod,
        timeseries,
        time_step,
        bandwidth
    )
    return nothing
end

iod(vec) = StatsBase.var(vec) / StatsBase.mean(vec)

function calculate_backward_index_of_dispersion!(
    index_of_dispersion_vec, timeseries, time_step, bandwidth
)
    return calculate_backward_metric!(
        index_of_dispersion_vec, iod, timeseries, time_step, bandwidth
    )
end

@testitem "Index of Dispersion" begin
    using CSDNoise
    using StatsBase

    daily_testpositives = collect(1:10)

    centered_iod_testpostives = calculate_ews_metric(
        calculate_centered_index_of_dispersion!, daily_testpositives, 1, 3
    )

    @test isequal(
        centered_iod_testpostives,
        [
            iod([1]),
            iod([1, 2, 3]),
            iod([1, 2, 3, 4, 5]),
            iod([2, 3, 4, 5, 6]),
            iod([3, 4, 5, 6, 7]),
            iod([4, 5, 6, 7, 8]),
            iod([5, 6, 7, 8, 9]),
            iod([6, 7, 8, 9, 10]),
            iod([8, 9, 10]),
            iod([10]),
        ],
    )

    backward_iod_testpostives = calculate_ews_metric(
        calculate_backward_index_of_dispersion!, daily_testpositives, 1, 3
    )

    @test isequal(
        backward_iod_testpostives,
        [
            iod([1]),
            iod([1, 2]),
            iod([1, 2, 3]),
            iod([2, 3, 4]),
            iod([3, 4, 5]),
            iod([4, 5, 6]),
            iod([5, 6, 7]),
            iod([6, 7, 8]),
            iod([7, 8, 9]),
            iod([8, 9, 10]),
        ],
    )
end

function calculate_centered_skewness!(
    skew_vec, timeseries, time_step, bandwidth
)
    calculate_centered_metric!(
        skew_vec, StatsBase.skewness, timeseries, time_step, bandwidth
    )
    return skew_vec
end

function calculate_backward_skewness!(
    skew_vec, timeseries, time_step, bandwidth
)
    return calculate_backward_metric!(
        skew_vec, StatsBase.skewness, timeseries, time_step, bandwidth
    )
end

@testitem "Skewness" begin
    using CSDNoise
    using StatsBase

    daily_testpositives = collect(1:10)

    centered_skewness_testpostives = calculate_ews_metric(
        calculate_centered_skewness!, daily_testpositives, 1, 3
    )

    @test isequal(
        centered_skewness_testpostives,
        [
            StatsBase.skewness([1]),
            StatsBase.skewness([1, 2, 3]),
            StatsBase.skewness([1, 2, 3, 4, 5]),
            StatsBase.skewness([2, 3, 4, 5, 6]),
            StatsBase.skewness([3, 4, 5, 6, 7]),
            StatsBase.skewness([4, 5, 6, 7, 8]),
            StatsBase.skewness([5, 6, 7, 8, 9]),
            StatsBase.skewness([6, 7, 8, 9, 10]),
            StatsBase.skewness([8, 9, 10]),
            StatsBase.skewness([10]),
        ],
    )

    backward_skewness_testpostives = calculate_ews_metric(
        calculate_backward_skewness!, daily_testpositives, 1, 3
    )

    @test isequal(
        backward_skewness_testpostives,
        [
            StatsBase.skewness([1]),
            StatsBase.skewness([1, 2]),
            StatsBase.skewness([1, 2, 3]),
            StatsBase.skewness([2, 3, 4]),
            StatsBase.skewness([3, 4, 5]),
            StatsBase.skewness([4, 5, 6]),
            StatsBase.skewness([5, 6, 7]),
            StatsBase.skewness([6, 7, 8]),
            StatsBase.skewness([7, 8, 9]),
            StatsBase.skewness([8, 9, 10]),
        ],
    )
end

function calculate_centered_kurtosis!(
    kurtosis_vec, timeseries, time_step, bandwidth
)
    calculate_centered_metric!(
        kurtosis_vec, kurtosis, timeseries, time_step, bandwidth
    )
    return kurtosis_vec
end

function calculate_backward_kurtosis!(
    kurtosis_vec, timeseries, time_step, bandwidth
)
    return calculate_backward_metric!(
        kurtosis_vec, kurtosis, timeseries, time_step, bandwidth
    )
end

#TODO: Confirm equation
function kurtosis(v)
    m = StatsBase.mean(v)
    n = length(v)
    cm2 = 0.0  # empirical 2nd centered moment (variance)
    cm4 = 0.0  # empirical 4th centered moment
    for i in 1:n
        @inbounds z = v[i] - m
        z2 = z * z
        cm2 += z2
        cm4 += z2 * z2
    end
    cm4 /= n
    cm2 /= n
    return (cm4 / (cm2 * cm2))
end

@testitem "Kurtosis" begin
    using CSDNoise

    daily_testpositives = collect(1:10)

    centered_kurtosis_testpostives = calculate_ews_metric(
        calculate_centered_kurtosis!, daily_testpositives, 1, 3
    )

    @test isequal(
        centered_kurtosis_testpostives,
        [
            kurtosis([1]),
            kurtosis([1, 2, 3]),
            kurtosis([1, 2, 3, 4, 5]),
            kurtosis([2, 3, 4, 5, 6]),
            kurtosis([3, 4, 5, 6, 7]),
            kurtosis([4, 5, 6, 7, 8]),
            kurtosis([5, 6, 7, 8, 9]),
            kurtosis([6, 7, 8, 9, 10]),
            kurtosis([8, 9, 10]),
            kurtosis([10]),
        ],
    )

    backward_kurtosis_testpostives = calculate_ews_metric(
        calculate_backward_kurtosis!, daily_testpositives, 1, 3
    )

    @test isequal(
        backward_kurtosis_testpostives,
        [
            kurtosis([1]),
            kurtosis([1, 2]),
            kurtosis([1, 2, 3]),
            kurtosis([2, 3, 4]),
            kurtosis([3, 4, 5]),
            kurtosis([4, 5, 6]),
            kurtosis([5, 6, 7]),
            kurtosis([6, 7, 8]),
            kurtosis([7, 8, 9]),
            kurtosis([8, 9, 10]),
        ],
    )
end

#TODO: Implement
function calculate_centered_autocorrelation!(
    autocorrelation_vec, timeseries, time_step, bandwidth
)
    return nothing
end

#TODO: Implement
function calculate_backward_autocorrelation!(
        autocorrelation_vec, timeseries, time_step, bandwidth
)
    return nothing
end

#TODO: Implement
function autocorrelation(v)
    return nothing
end

@testitem "Autocorrelation" begin
    using CSDNoise
    using StatsBase

    daily_testpositives = collect(1:10)

    centered_autocorrelation_testpostives = calculate_ews_metric(
        calculate_centered_autocorrelation!, daily_testpositives, 1, 3
    )

    @test isequal(
        centered_autocorrelation_testpostives,
        [
            autocorrelation([1]),
            autocorrelation([1, 2, 3]),
            autocorrelation([1, 2, 3, 4, 5]),
            autocorrelation([2, 3, 4, 5, 6]),
            autocorrelation([3, 4, 5, 6, 7]),
            autocorrelation([4, 5, 6, 7, 8]),
            autocorrelation([5, 6, 7, 8, 9]),
            autocorrelation([6, 7, 8, 9, 10]),
            autocorrelation([8, 9, 10]),
            autocorrelation([10]),
        ],
    )

    backward_autocorrelation_testpostives = calculate_ews_metric(
        calculate_backward_autocorrelation!, daily_testpositives, 1, 3
    )

    @test isequal(
        backward_autocorrelation_testpostives,
        [
            autocorrelation([1]),
            autocorrelation([1, 2]),
            autocorrelation([1, 2, 3]),
            autocorrelation([2, 3, 4]),
            autocorrelation([3, 4, 5]),
            autocorrelation([4, 5, 6]),
            autocorrelation([5, 6, 7]),
            autocorrelation([6, 7, 8]),
            autocorrelation([7, 8, 9]),
            autocorrelation([8, 9, 10]),
        ],
    )
end
