using SumTypes
using DataFrames: DataFrames

@sum_type EWSMethod begin
    Backward
    Centered
end

struct SpaeroEWSMetricSpecification
    method::EWSMethod
    bandwidth::Int
    lag::Int
end

struct SpaeroEWSMetrics{
    T1<:AbstractFloat,T2<:SpaeroEWSMetricSpecification,
    T3<:AbstractArray{<:AbstractFloat},
}
    timestep::T1
    ews_specification::T2
    mean::T3
    variance::T3
    coefficient_of_variation::T3
    index_of_dispersion::T3
    skewness::T3
    kurtosis::T3
    autocovariance::T3
    autocorrelation::T3
end

function SpaeroEWSMetrics(
    ews_spec::SpaeroEWSMetricSpecification, timeseries, time_step
)
    mean_vec = spaero_mean(ews_spec.method, timeseries, ews_spec.bandwidth)
    var_vec = spaero_var(
        ews_spec.method, mean_vec, timeseries, ews_spec.bandwidth
    )
    sd_vec = sqrt.(var_vec)
    m3_vec = _spaero_moment(
        ews_spec.method, mean_vec, timeseries, 3, ews_spec.bandwidth
    )
    m4_vec = _spaero_moment(
        ews_spec.method, mean_vec, timeseries, 4, ews_spec.bandwidth
    )
    autocov_vec = spaero_autocov(
        ews_spec.method,
        mean_vec,
        timeseries,
        ews_spec.bandwidth;
        lag = ews_spec.lag,
    )
    return SpaeroEWSMetrics(
        time_step,
        ews_spec,
        mean_vec,
        var_vec,
        spaero_cov(var_vec, mean_vec),
        spaero_iod(var_vec, mean_vec),
        spaero_skew(m3_vec, sd_vec),
        spaero_kurtosis(m4_vec, var_vec),
        autocov_vec,
        spaero_autocor(autocov_vec, sd_vec),
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

function _spaero_moment(method::EWSMethod, timeseries, moment, bandwidth)
    return _spaero_moment(
        spaero_mean(method, timeseries, bandwidth), timeseries, moment,
        bandwidth,
    )
end

function _spaero_moment(
    method::EWSMethod, mean_timeseries, timeseries, moment, bandwidth
)
    diff = (timeseries .- mean_timeseries) .^ moment
    return spaero_mean(method, diff, bandwidth)
end

function spaero_var(method::EWSMethod, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    return _spaero_moment(method, mean_vec, timeseries, 2, bandwidth)
end

function spaero_var(method::EWSMethod, mean_vec, timeseries, bandwidth)
    return _spaero_moment(method, mean_vec, timeseries, 2, bandwidth)
end

function spaero_cov(method::EWSMethod, timeseries, bandwidth)
    var_vec = spaero_var(method, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    return spaero_cov(var_vec, mean_vec)
end

function spaero_cov(var_vec, mean_vec)
    return sqrt.(var_vec) ./ mean_vec
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
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    m3_vec = _spaero_moment(method, mean_vec, timeseries, 3, bandwidth)
    var_vec = spaero_var(method, timeseries, bandwidth)
    return spaero_skew(m3_vec, sqrt.(var_vec))
end

function spaero_skew(m3_vec, sd_vec)
    return m3_vec ./ sd_vec .^ 3
end

function spaero_kurtosis(method::EWSMethod, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    m4_vec = _spaero_moment(method, mean_vec, timeseries, 4, bandwidth)
    var_vec = spaero_var(method, timeseries, bandwidth)
    return spaero_kurtosis(m4_vec, var_vec)
end

function spaero_kurtosis(m4_vec, var_vec)
    return m4_vec ./ var_vec .^ 2
end

function spaero_autocov(method::EWSMethod, timeseries, bandwidth; lag = 1)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    autocov_vec = spaero_autocov(
        method, mean_vec, timeseries, bandwidth; lag = lag
    )
    return autocov_vec
end

function spaero_autocov(
    method::EWSMethod,
    mean_vec,
    timeseries,
    bandwidth;
    lag = 1
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
    var_vec = spaero_var(method, timeseries, bandwidth)
    sd_vec = sqrt.(var_vec)
    autocov_vec = spaero_autocov(method, timeseries, bandwidth; lag = lag)
    return spaero_autocor(autocov_vec, sd_vec; lag = lag)
end

function spaero_autocor(autocov_vec, sd_vec; lag = 1)
    lagged_sd = zeros(Float64, length(sd_vec))
    @inbounds for i in eachindex(sd_vec)
        if i <= lag
            lagged_sd[i] = NaN
            continue
        end
        lagged_sd[i] = sd_vec[i - lag]
    end
    return autocov_vec ./ (sd_vec .* lagged_sd)
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
    showdiffs = true,
) where {T1<:DataFrames.DataFrame,T2<:SpaeroEWSMetrics}
    warnings = 0
    for metric in ews
        spaero = getproperty(spaero_ews, metric)
        my = getproperty(my_ews, metric)

        diff = abs.(spaero .- my)

        filtered_diff = filter(x -> x > tolerance, skipmissing(diff))

        if length(filtered_diff) > 0
            warnings += 1
            if showwarnings
                println()
                @warn "There are differences in the spaero and my implementation of the $metric EWS."
                if showdiffs
                    println(
                        filter_spaero_comparison(
                            DataFrames.DataFrame(
                                [spaero, my, diff], [:spaero, :mine, :absdiff]
                            );
                            tolerance = tolerance,
                            warn = false,
                        ),
                    )
                end
            end
        end
    end
    println()
    if warnings > 0
        @warn "🟡 There were warnings in $warnings metrics 🟡"
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
