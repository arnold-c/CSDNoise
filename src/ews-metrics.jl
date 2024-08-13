using SumTypes
using DataFrames: DataFrames
using Bumper
using StrideArrays
using AllocCheck
using Printf: @sprintf
using TestItems

function EWSMetrics(
    ews_spec::EWSMetricSpecification, timeseries
)
    aggregated_timeseries = aggregate_timeseries(
        timeseries, ews_spec.aggregation
    )

    mean_vec = spaero_mean(
        ews_spec.method, aggregated_timeseries, ews_spec.bandwidth
    )
    var_vec = spaero_var(
        ews_spec.method, mean_vec, aggregated_timeseries, ews_spec.bandwidth
    )
    var2_vec = var_vec .^ 2
    sd_vec = sqrt.(var_vec)
    sd3_vec = sd_vec .^ 3
    m3_vec = _spaero_moment(
        ews_spec.method, mean_vec, aggregated_timeseries, 3, ews_spec.bandwidth
    )
    m4_vec = _spaero_moment(
        ews_spec.method, mean_vec, aggregated_timeseries, 4, ews_spec.bandwidth
    )
    autocov_vec = spaero_autocov(
        ews_spec.method,
        mean_vec,
        aggregated_timeseries,
        ews_spec.bandwidth;
        lag = ews_spec.lag,
    )
    return EWSMetrics(
        ews_spec,
        mean_vec,
        var_vec,
        spaero_cov(sd_vec, mean_vec),
        spaero_iod(var_vec, mean_vec),
        spaero_skew(m3_vec, sd3_vec),
        spaero_kurtosis(m4_vec, var2_vec),
        autocov_vec,
        spaero_autocor(autocov_vec, sd_vec),
    )
end

function aggregate_timeseries(timeseries, aggregation)
    if aggregation == 1
        return timeseries
    end
    return _aggregate_timeseries(timeseries, aggregation)
end

function _aggregate_timeseries(timeseries, aggregation)
    aggregate_timeseries = zeros(
        eltype(timeseries), length(timeseries) ÷ aggregation
    )
    for i in eachindex(aggregate_timeseries)
        aggregate_timeseries[i] = sum(
            @view(timeseries[((i - 1) * aggregation + 1):(i * aggregation)])
        )
    end
    return aggregate_timeseries
end

@testitem "Timeseries aggregation" begin
    using CSDNoise

    testvec = collect(1:10)
    @test isequal(aggregate_timeseries(testvec, 1), testvec)
    @test isequal(
        aggregate_timeseries(testvec, 2),
        [sum(1:2), sum(3:4), sum(5:6), sum(7:8), sum(9:10)],
    )
    @test isequal(
        aggregate_timeseries(testvec, 3),
        [sum(1:3), sum(4:6), sum(7:9)],
    )
end

function spaero_mean(method::EWSMethod, timeseries, bandwidth)
    mean_vec = zeros(Float64, length(timeseries))
    spaero_mean!(mean_vec, method, timeseries, bandwidth)
    return mean_vec
end

function spaero_mean!(
    mean_vec, method::EWSMethod, timeseries, bandwidth
)
    mean_func! = _get_mean_func(method)
    mean_func!(mean_vec, timeseries, bandwidth)
    return nothing
end

function _get_mean_func(method::EWSMethod)
    return @cases method begin
        Backward => _spaero_backward_mean!
        Centered => _spaero_centered_mean!
    end
end

function _spaero_centered_mean!(mean_vec, timeseries, bw)
    tlength = length(timeseries)
    @inbounds for i in eachindex(timeseries)
        if i < bw
            mean_vec[i] = mean(@view(timeseries[begin:(i + bw - 1)]))
        elseif i + bw > tlength
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):end]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):(i + bw - 1)]))
        end
    end
    return nothing
end

function _spaero_backward_mean!(
    mean_vec, timeseries, bw
)
    @inbounds for i in eachindex(timeseries)
        if i < bw
            mean_vec[i] = mean(@view(timeseries[begin:i]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):i]))
        end
    end
    return nothing
end

function _spaero_moment(
    method::EWSMethod, timeseries, moment, bandwidth
)
    return _spaero_moment(
        method,
        spaero_mean(method, timeseries, bandwidth),
        timeseries,
        moment,
        bandwidth,
    )
end

function _spaero_moment(
    method::EWSMethod, mean_timeseries, timeseries, moment, bandwidth
)
    @no_escape begin
        diff = @alloc(Float64, length(timeseries))
        diff .= (timeseries .- mean_timeseries) .^ moment
        spaero_mean(method, diff, bandwidth)
    end
end

function spaero_var(method::EWSMethod, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    return _spaero_moment(method, mean_vec, timeseries, 2, bandwidth)
end

function spaero_var(method::EWSMethod, mean_vec, timeseries, bandwidth)
    return _spaero_moment(method, mean_vec, timeseries, 2, bandwidth)
end

function spaero_cov(method::EWSMethod, timeseries, bandwidth)
    @no_escape begin
        mean_vec = @alloc(Float64, length(timeseries))
        var_vec = @alloc(Float64, length(timeseries))
        mean_vec = spaero_mean(method, timeseries, bandwidth)
        var_vec = spaero_var(method, mean_vec, timeseries, bandwidth)
        sd_vec = @alloc(Float64, length(var_vec))
        sd_vec .= sqrt.(var_vec)
        spaero_cov(sd_vec, mean_vec)
    end
end

function spaero_cov(sd_vec, mean_vec)
    return sd_vec ./ mean_vec
end

function spaero_iod(method::EWSMethod, timeseries, bandwidth)
    var_vec = spaero_var(method, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    return spaero_iod(var_vec, mean_vec)
end

function spaero_iod(var_vec, mean_vec)
    return var_vec ./ mean_vec
end

function spaero_skew(method::EWSMethod, timeseries, bandwidth)
    @no_escape begin
        # use for both mean and sd vec as only one needed at a time
        worker_vec = @alloc(Float64, length(timeseries))
        m3_vec = @alloc(Float64, length(timeseries))
        sd3_vec = @alloc(Float64, length(timeseries))
        worker_vec .= spaero_mean(method, timeseries, bandwidth)
        m3_vec .= _spaero_moment(method, worker_vec, timeseries, 3, bandwidth)
        sd3_vec .=
            sqrt.(spaero_var(method, worker_vec, timeseries, bandwidth)) .^ 3
        spaero_skew(m3_vec, sd3_vec)
    end
end

function spaero_skew(m3_vec, sd3_vec)
    return m3_vec ./ sd3_vec
end

function spaero_kurtosis(method::EWSMethod, timeseries, bandwidth)
    @no_escape begin
        mean_vec = @alloc(Float64, length(timeseries))
        m4_vec = @alloc(Float64, length(timeseries))
        var2_vec = @alloc(Float64, length(timeseries))
        mean_vec .= spaero_mean(method, timeseries, bandwidth)
        m4_vec .= _spaero_moment(method, mean_vec, timeseries, 4, bandwidth)
        var2_vec .= spaero_var(method, mean_vec, timeseries, bandwidth) .^ 2
        spaero_kurtosis(m4_vec, var2_vec)
    end
end

function spaero_kurtosis(m4_vec, var2_vec)
    return m4_vec ./ var2_vec
end

function spaero_autocov(method::EWSMethod, timeseries, bandwidth; lag = 1)
    @no_escape begin
        mean_vec = @alloc(Float64, length(timeseries))
        mean_vec .= spaero_mean(method, timeseries, bandwidth)
        spaero_autocov(
            method, mean_vec, timeseries, bandwidth; lag = lag
        )
    end
end

function spaero_autocov(
    method::EWSMethod,
    mean_vec,
    timeseries,
    bandwidth;
    lag = 1,
)
    meandiff = timeseries .- mean_vec
    autocov_vec = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i <= lag
            continue
        end
        autocov_vec[i] = meandiff[i] * meandiff[i - lag]
    end
    autocov_vec = spaero_mean(method, autocov_vec, bandwidth)
    autocov_vec[begin:lag] .= NaN
    return autocov_vec
end

function spaero_autocor(method::EWSMethod, timeseries, bandwidth; lag = 1)
    sd2 = zeros(Float64, length(timeseries))
    @no_escape begin
        var_vec = @alloc(Float64, length(timeseries))
        sd_vec = @alloc(Float64, length(timeseries))
        autocov_vec = @alloc(Float64, length(timeseries))
        var_vec .= spaero_var(method, timeseries, bandwidth)
        sd_vec .= sqrt.(var_vec)
        autocov_vec .= spaero_autocov(method, timeseries, bandwidth; lag = lag)
        lagged_sd = @alloc(Float64, length(sd_vec))
        _lagged_vector(lagged_sd, sd_vec, lag)
        sd2 .= sd_vec .* lagged_sd
    end
    return autocov_vec ./ sd2
end

function spaero_autocor(autocov_vec, sd_vec; lag = 1)
    sd2 = zeros(Float64, length(sd_vec))
    @no_escape begin
        lagged_sd = @alloc(Float64, length(sd_vec))
        _lagged_vector(lagged_sd, sd_vec, lag)
        sd2 .= sd_vec .* lagged_sd
    end
    return autocov_vec ./ sd2
end

function _lagged_vector(lagged_vec, vec, lag)
    @inbounds for i in eachindex(lagged_vec)
        if i <= lag
            lagged_vec[i] = NaN
            continue
        end
        lagged_vec[i] = vec[i - lag]
    end
    return nothing
end

function compare_against_spaero(
    spaero_ews::T1, my_ews::T2;
    ews = [
        :autocorrelation,
        :autocovariance,
        :coefficient_of_variation,
        :index_of_dispersion,
        :kurtosis,
        :mean,
        :skewness,
        :variance,
    ],
    tolerance = 1e-13,
    showwarnings = true,
    showdiffs = false,
) where {T1<:DataFrames.DataFrame,T2<:EWSMetrics}
    warnings = 0
    maxabsdiff = 0.0
    warning_metric = :none
    for metric in ews
        spaero = getproperty(spaero_ews, metric)
        my = getproperty(my_ews, metric)

        diff = abs.(spaero .- my)

        filtered_diff = filter(x -> x > tolerance, skipmissing(diff))

        if length(filtered_diff) > 0
            warnings += 1
            if maximum(filtered_diff) > maxabsdiff
                maxabsdiff = maximum(filtered_diff)
                warning_metric = metric
            end
            if showdiffs
                println()
                @warn "There are differences in the spaero and my implementation of the $metric EWS."
                println(
                    filter_spaero_comparison(
                        DataFrames.DataFrame(
                            [spaero, my, diff], [:spaero, :mine, :absdiff]
                        );
                        tolerance = tolerance,
                        warn = false,
                    ),
                )
                showwarnings = false
            end
            if showwarnings
                println()
                @warn "There are differences in the spaero and my implementation of the $metric EWS."
            end
        end
    end
    println()
    if warnings > 0
        maxabsdiff = @sprintf("%.4E", maxabsdiff)
        @warn "🟡 There were warnings in $warnings metrics for the $(my_ews.ews_specification.method) 🟡\n The max absolute difference was $maxabsdiff and occured in $warning_metric"
    else
        @info "✅ There were no warnings - all metrics are within tolerance ✅"
    end
end

function compare_against_spaero(spaero_ews, my_ews)
    df = DataFrames.DataFrame([spaero_ews my_ews], [:spaero, :mine])
    df.absdiff = abs.(df.spaero .- df.mine)
    return df
end

function filter_spaero_comparison(spaero_ews, my_ews; tolerance = 1e-13)
    df = compare_against_spaero(spaero_ews, my_ews)
    return filter_spaero_comparison(df; tolerance = tolerance)
end

function filter_spaero_comparison(df; tolerance = 1e-13, warn = true)
    subsetted = DataFrames.subset(
        df, :absdiff => x -> x .> tolerance; skipmissing = true
    )
    if DataFrames.nrow(subsetted) > 0 && warn
        println()
        @warn "There are differences in the metrics between spaero and my implementation of EWS."
    end
    return subsetted
end
