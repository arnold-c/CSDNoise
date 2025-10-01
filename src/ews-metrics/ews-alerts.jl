export calculate_ews_trigger_index

function calculate_ews_trigger_index(
        ews_thresholds::T1,
        consecutive_thresholds = 2,
    ) where {T1 <: AbstractMatrix{<:Bool}}
    reshaped_ews_thresholds = reshape(ews_thresholds, :)

    @assert length(reshaped_ews_thresholds) == length(ews_thresholds[:, 1])

    return calculate_ews_trigger_index(
        reshaped_ews_thresholds,
        consecutive_thresholds,
    )
end

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
