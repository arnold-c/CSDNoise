export calculate_ews_lead_time

function calculate_ews_lead_time(
        ews_thresholds;
        week_aggregation = 1,
        consecutive_thresholds = 2,
        output_type = :days,
    )
    threshold_index = Try.@? calculate_ews_trigger_index(
        ews_thresholds, consecutive_thresholds
    )

    return calculate_ews_lead_time(
        ews_thresholds,
        threshold_index;
        week_aggregation = week_aggregation,
        output_type = output_type,
    )
end

function calculate_ews_lead_time(
        ews_thresholds, threshold_index;
        week_aggregation = 1,
        output_type = :days,
    )
    output_multiplier = if output_type == :days
        7
    elseif output_type == :weeks
        1
    elseif output_type == :months
        7 / 30.5
    elseif output_type == :years
        7 / 365
    else
        error(
            "Unknown output type: $(output_type).\nChoose between :days, :weeks, :months, or :years."
        )
    end

    if isnothing(threshold_index)
        @warn "No ews trigger index found. Returning nothing."
        return nothing
    end

    return (length(ews_thresholds) - threshold_index) * week_aggregation *
        output_multiplier
end
