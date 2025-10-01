export calculate_above_threshold_bounds

function calculate_above_threshold_bounds(timeseries_rle)
    # Calculate upper and lower indices of consecutive days of infection
    timeseries_rle_accum = accumulate(+, timeseries_rle[2])
    upperbound_indices = findall(isequal(1), timeseries_rle[1])

    upper_bounds = Vector{Int64}(undef, length(upperbound_indices))
    lower_bounds = similar(upper_bounds)
    duration = similar(upper_bounds)

    @inbounds upper_bounds .= @view(
        timeseries_rle_accum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : timeseries_rle_accum[x - 1] + 1,
        lower_bounds,
        upperbound_indices,
    )
    duration .= upper_bounds .- lower_bounds .+ 1

    return Thresholds(lower_bounds, upper_bounds, duration)
end
