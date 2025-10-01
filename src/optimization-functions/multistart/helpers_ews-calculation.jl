export calculate_ews_metrics_for_simulation

# TODO: Remove as not needed with new structvector approach with filtering
"""
    calculate_ews_metrics_for_simulation(ews_metric_specification, testarr_view, null_testarr_view, threshold, ews_enddate_type)

Calculate EWS metrics for a single simulation given the test array views and threshold.
Returns a tuple of (ews_metrics, null_ews_metrics) wrapped in Try types.
"""
function calculate_ews_metrics_for_simulation(
        ews_metric_specification,
        testarr_view,
        threshold,
        ews_enddate_val
    )
    # Calculate EWS metrics
    ews_vals = EWSMetrics(
        ews_metric_specification,
        @view(testarr_view[1:ews_enddate_val, 5])
    )

    return ews_vals
end
