export get_ews_metric_specification_description,
    method_string

function get_ews_metric_specification_description(ews_metric_specification)
    return "Method: $(method_string(ews_metric_specification.method)), Aggregation: $(ews_metric_specification.aggregation), Bandwidth: $(ews_metric_specification.bandwidth), Lag: $(ews_metric_specification.lag)"
end

method_string(method::EWSMethod) = lowercase(split(string(method), "::")[1])
