export calculate_bandwidth_and_return_ews_metric_spec

function calculate_bandwidth_and_return_ews_metric_spec(
        ews_metric_spec_components...
    )
    @assert length(ews_metric_spec_components) == 4

    bandwidth = calculate_bandwidth(
        ews_metric_spec_components[3],
        ews_metric_spec_components[2],
    )

    return EWSMetricSpecification(
        ews_metric_spec_components[1],
        ews_metric_spec_components[2],
        bandwidth,
        ews_metric_spec_components[4],
    )
end

function calculate_bandwidth(bandwidth_days, aggregation_days)
    return bandwidth_days ÷ aggregation_days
end
