using SumTypes

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
    return SpaeroEWSMetrics(
        time_step,
        ews_spec,
        spaero_mean(ews_spec.method, timeseries, ews_spec.bandwidth),
        spaero_var(ews_spec.method, timeseries, ews_spec.bandwidth),
        spaero_cov(ews_spec.method, timeseries, ews_spec.bandwidth),
        spaero_iod(ews_spec.method, timeseries, ews_spec.bandwidth),
        spaero_skew(ews_spec.method, timeseries, ews_spec.bandwidth),
        spaero_kurtosis(ews_spec.method, timeseries, ews_spec.bandwidth),
        spaero_autocov(
            ews_spec.method,
            timeseries,
            ews_spec.bandwidth;
            lag = ews_spec.lag
        ),
        spaero_autocor(
            ews_spec.method,
            timeseries,
            ews_spec.bandwidth;
            lag = ews_spec.lag
        ),
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

function spaero_cov(method::EWSMethod, timeseries, bandwidth)
    var = spaero_var(method, timeseries, bandwidth)
    mean = spaero_mean(method, timeseries, bandwidth)
    return sqrt.(var) ./ mean
end

# function spaero_cov(method::EWSMethod, timeseries, bandwidth)
#     mean_func, moment_func = window_functions(method)
#     return sqrt.(moment_func(timeseries, 2, bandwidth)) ./
#            mean_func(timeseries, bandwidth)
# end

function spaero_iod(method::EWSMethod, timeseries, bandwidth)
    var = spaero_var(method, timeseries, bandwidth)
    mean = spaero_mean(method, timeseries, bandwidth)
    return var ./ mean
end

function spaero_skew(method::EWSMethod, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    m3 = _spaero_moment(method, mean_vec, timeseries, 3, bandwidth)
    var = spaero_var(method, timeseries, bandwidth)
    return m3 ./ var .^ 1.5
end

function spaero_kurtosis(method::EWSMethod, timeseries, bandwidth)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    m4 = _spaero_moment(method, mean_vec, timeseries, 4, bandwidth)
    var = spaero_var(method, timeseries, bandwidth)

    return m4 ./ var .^ 2
end

function spaero_autocov(method::EWSMethod, timeseries, bandwidth; lag = 1)
    mean_vec = spaero_mean(method, timeseries, bandwidth)
    meandiff = timeseries .- mean_vec
    worker = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i <= lag
            continue
        end
        worker[i] = meandiff[i] * meandiff[i - lag]
    end
    worker = spaero_mean(method, worker, bandwidth)
    worker[begin:lag] .= NaN
    return worker
end

function spaero_autocor(method::EWSMethod, timeseries, bandwidth; lag = 1)
    var = spaero_var(method, timeseries, bandwidth)
    sd = sqrt.(var)
    lagged_sd = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i <= 1
            lagged_sd[i] = NaN
            continue
        end
        lagged_sd[i] = sd[i - lag]
    end
    return spaero_autocov(method, timeseries, bandwidth) ./ (sd .* lagged_sd)
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
    show = true,
) where {T1<:DataFrame,T2<:SpaeroEWSMetrics}
    for metric in ews
        spaero = getproperty(spaero_ews, metric)
        my = getproperty(my_ews, metric)

        diff = abs.(spaero .- my)

        filtered_diff = filter(x -> x > tolerance, skipmissing(diff))

        if length(filtered_diff) > 0
            println()
            @warn "There are differences in the spaero and my implementation of the $metric EWS."
            if show
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
