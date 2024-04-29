using SumTypes

function spaero_centered_mean(timeseries, bw)
    tlength = length(timeseries)
    mean_vec = zeros(Float64, length(timeseries))
    spaero_centered_mean!(mean_vec, tlength, timeseries, bw)
    return mean_vec
end

function spaero_centered_mean!(mean_vec, tlength, timeseries, bw)
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

function spaero_centered_moment(timeseries, moment, bw)
    return spaero_centered_moment(
        spaero_centered_mean(timeseries, bw), timeseries, moment, bw
    )
end

function spaero_centered_moment(mean_timeseries, timeseries, moment, bw)
    diff = (timeseries .- mean_timeseries) .^ moment
    return spaero_centered_mean(diff, bw)
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

function filter_spaero_comparison(df; tolerance = 1e-13)
    return DataFrames.subset(
        df, :absdiff => x -> x .> tolerance; skipmissing = true
    )
end

function spaero_backward_mean(timeseries, bw)
    tlength = length(timeseries)
    mean_vec = zeros(Float64, tlength)
    @inbounds for i in eachindex(timeseries)
        if i < bw
            mean_vec[i] = mean(@view(timeseries[begin:i]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):i]))
        end
    end
    return mean_vec
end

function spaero_backward_moment(timeseries, moment, bw)
    return spaero_backward_moment(
        timeseries, spaero_backward_mean(timeseries, bw), moment, bw
    )
end

function spaero_backward_moment(timeseries, mean_timeseries, moment, bw)
    diff = (timeseries .- mean_timeseries) .^ moment
    return spaero_backward_mean(diff, bw)
end

function spaero_cov(mean_func, moment_func, timeseries, bandwidth)
    return sqrt.(moment_func(timeseries, 2, bandwidth)) ./
           mean_func(timeseries, bandwidth)
end

@sum_type EWSMethod begin
    Backward
    Centered
end

function window_functions(
    method::EWSMethod
)::Tuple{
    Union{typeof(spaero_centered_mean),typeof(spaero_backward_mean)},
    Union{typeof(spaero_centered_moment),typeof(spaero_backward_moment)},
}
    mean_func = @cases method begin
        Backward => spaero_backward_mean
        Centered => spaero_centered_mean
    end
    moment_func = @cases method begin
        Backward => spaero_backward_moment
        Centered => spaero_centered_moment
    end
    return mean_func, moment_func
end

function spaero_cov(method::EWSMethod, timeseries, bandwidth)
    mean_func, moment_func = window_functions(method)
    return sqrt.(moment_func(timeseries, 2, bandwidth)) ./
           mean_func(timeseries, bandwidth)
end

function spaero_iod(method::EWSMethod, timeseries, bandwidth)
    mean_func, moment_func = window_functions(method)
    return (moment_func(timeseries, 2, bandwidth)) ./
           mean_func(timeseries, bandwidth)
end

function spaero_skew(method::EWSMethod, timeseries, bandwidth)
    moment_func = window_functions(method)[2]
    return moment_func(timeseries, 3, bandwidth) ./
           moment_func(timeseries, 2, bandwidth) .^ 1.5
end

function spaero_kurtosis(method::EWSMethod, timeseries, bandwidth)
    moment_func = window_functions(method)[2]
    return moment_func(timeseries, 4, bandwidth) ./
           moment_func(timeseries, 2, bandwidth) .^ 2
end

function spaero_autocov(method::EWSMethod, timeseries, bandwidth)
    mean_func = window_functions(method)[1]
    meandiff = timeseries .- mean_func(timeseries, bandwidth)
    worker = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i == 1
            continue
        end
        worker[i] = meandiff[i] * meandiff[i - 1]
    end
    worker = mean_func(worker, bandwidth)
    worker[1] = NaN
    return worker
end

function spaero_autocor(method::EWSMethod, timeseries, bandwidth)
    mean_func, moment_func = window_functions(method)
    sd = sqrt.(moment_func(timeseries, 2, bandwidth))
    lagged_sd = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i == 1
            lagged_sd[i] = NaN
            continue
        end
        lagged_sd[i] = sd[i - 1]
    end
    return spaero_autocov(method, timeseries, bandwidth) ./ (sd .* lagged_sd)
end
