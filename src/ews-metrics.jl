using Statistics

function calculate_centered_mean()
end

function calculate_backward_mean(time_step, bandwidth, timeseries)
    mean_vec = zeros(Float64, length(timeseries))

    calculate_backward_mean!(mean_vec, time_step, bandwidth, timeseries)

    return mean_vec
end

function calculate_backward_mean!(mean_vec, time_step, bandwidth, timeseries)
    # TODO: implement
    delta = (bandwidth - 1) * time_step
    for i in eachindex(timeseries)
        mean_vec[i] = calculate_mean(
            @view(timeseries[(i - delta):(i + delta)])
        )
    end
    return mean_vec
end

function calculate_mean(timeseries)
    # TODO: implement
    return nothing
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
