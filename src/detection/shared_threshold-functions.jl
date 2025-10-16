export calculate_above_threshold_bounds

"""
    calculate_above_threshold_bounds(timeseries_rle)

Calculate the bounds and durations of consecutive periods where a condition is true.

This function takes a run-length encoded (RLE) representation of a boolean timeseries
and identifies all consecutive periods where the condition is true (value = 1). For each
such period, it calculates the lower bound (start index), upper bound (end index), and
duration of the period.

The function processes the RLE data by:
1. Computing cumulative positions using the run lengths
2. Finding all positions where the condition is true (value = 1)
3. Calculating upper bounds from cumulative positions
4. Computing lower bounds as the start of each consecutive period
5. Determining durations as the difference between bounds

# Arguments
- `timeseries_rle`: A tuple `(values, lengths)` from `StatsBase.rle()` where:
  - `values`: Vector of boolean values (0 or 1) indicating condition status
  - `lengths`: Vector of integers indicating run lengths for each value

# Returns
- `Thresholds`: A struct containing:
  - `lower_bounds`: Vector of starting indices for each above-threshold period
  - `upper_bounds`: Vector of ending indices for each above-threshold period
  - `duration`: Vector of durations (in time steps) for each above-threshold period

# Examples
```julia
# Example with Reff >= 1 detection
Reff_vec = [0.8, 0.9, 1.1, 1.2, 1.0, 0.9, 1.3, 1.1]
Reff_rle = StatsBase.rle(Reff_vec .>= 1)
thresholds = calculate_above_threshold_bounds(Reff_rle)
# Returns periods where Reff >= 1
```

# See Also
- [`Thresholds`](@ref): Return type containing threshold bounds and durations
- [`Reff_ge_than_one`](@ref): Function that uses this for Reff >= 1 detection
- [`StatsBase.rle`](@ref): Function for creating run-length encoded input
"""
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
