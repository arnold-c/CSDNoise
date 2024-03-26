using Statistics

function calculate_centered_mean(timeseries, time_step, bandwidth)
    mean_vec = zeros(Float64, length(timeseries))

    calculate_centered_mean!(mean_vec, time_step, bandwidth, timeseries)

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

function calculate_backward_mean(timeseries, time_step, bandwidth)
    mean_vec = zeros(Float64, length(timeseries))

    calculate_backward_mean!(mean_vec, time_step, bandwidth, timeseries)

    return mean_vec
end

function calculate_backward_mean!(mean_vec, timeseries, time_step, bandwidth)
    calculate_movingavg!(mean_vec, timeseries, bandwidth * time_step)
    return mean_vec
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
