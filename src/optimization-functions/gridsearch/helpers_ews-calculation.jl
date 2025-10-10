export generate_ensemble_ews_metrics

function generate_ensemble_ews_metrics(
        ews_metric_specification,
        test_positive_results::EnsembleTestResultRun,
    )
    emergent_test_positives = test_positive_results.emergent_test_positives
    null_test_positives = test_positive_results.null_test_positives
    nsims = length(emergent_test_positives)

    emergent_ews_metrics = Vector{EWSMetrics}(undef, nsims)
    null_ews_metrics = similar(emergent_ews_metrics)

    for i in eachindex(emergent_test_positives)
        emergent_ews_metrics[i] = EWSMetrics(
            ews_metric_specification,
            emergent_test_positives[i],
        )
        null_ews_metrics[i] = EWSMetrics(
            ews_metric_specification,
            null_test_positives[i],
        )
    end

    return EnsembleEWSMetrics(
        emergent_ews_metrics = StructVector(emergent_ews_metrics),
        null_ews_metrics = StructVector(null_ews_metrics)
    )
end
