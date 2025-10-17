export calculate_ews_trigger_index

"""
    calculate_ews_trigger_index(ews_thresholds::AbstractVector{Bool}, consecutive_thresholds=2)

Core implementation for calculating EWS trigger index from vector of threshold exceedances.

This method implements the main algorithm for detecting consecutive threshold
exceedances using an efficient cumulative sum approach. It processes a boolean
vector representing threshold exceedances and identifies the first occurrence
of the required number of consecutive exceedances.

# Arguments
- `ews_thresholds::AbstractVector{Bool}`: Vector of boolean threshold exceedance indicators
- `consecutive_thresholds::Int=2`: Number of consecutive exceedances required

# Returns
- `Try.Ok{Int}`: Index where consecutive exceedances first occur
- `Try.Err{Nothing}`: Error if insufficient consecutive exceedances found

# Algorithm Details
1. Compute cumulative sum of boolean exceedances
2. For each index i > consecutive_thresholds:
   - Check if cumulative_sum[i] - cumulative_sum[i - consecutive_thresholds] >= consecutive_thresholds
   - This efficiently detects consecutive_thresholds consecutive `true` values
3. Return first index meeting the condition

# Performance
- Time complexity: O(n) where n is the length of ews_thresholds
- Space complexity: O(n) for the cumulative sum vector

# Example
```julia
# Detect 3 consecutive exceedances
exceedances = [false, true, false, true, true, true, false]
result = calculate_ews_trigger_index(exceedances, 3)
# Returns Try.Ok(6) - first index after 3 consecutive trues at positions 4,5,6
```

# See Also
- [`calculate_ews_trigger_index`](@ref): Matrix method that delegates to this implementation
"""
function calculate_ews_trigger_index(
        ews_thresholds::T1,
        consecutive_thresholds = 2,
    ) where {T1 <: AbstractVector{<:Bool}}
    cumulative_thresholds = cumsum(ews_thresholds)
    for (i, v) in pairs(cumulative_thresholds)
        if i > consecutive_thresholds &&
                v >=
                cumulative_thresholds[i - consecutive_thresholds] +
                consecutive_thresholds
            return Try.Ok(i)
        end
    end
    return Try.Err(nothing)
end
