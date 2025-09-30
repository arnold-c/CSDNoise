export generate_ensemble_ews_metrics

function generate_ensemble_ews_metrics(
        data_arrs::NamedTuple{T},
        testarr,
        null_testarr,
        ews_metric_specification,
        ews_enddate_type,
    ) where {T}
    thresholds = get_enddate_thresholds(data_arrs, ews_enddate_type)

    ensemble_nsims = size(testarr, 3)
    ews_metrics = Vector{EWSMetrics}()
    null_ews_metrics = Vector{EWSMetrics}()

    for sim in 1:ensemble_nsims
        ews_enddate = calculate_ews_enddate(thresholds[sim], ews_enddate_type)
        if Try.isok(ews_enddate)
            ews_enddate_val = Try.unwrap(ews_enddate)

            ews_metric = calculate_ews_metrics_for_simulation(
                ews_metric_specification,
                @view(testarr[:, :, sim]),
                thresholds[sim],
                ews_enddate_val
            )
            push!(ews_metrics, ews_metric)

            null_ews_metric = calculate_ews_metrics_for_simulation(
                ews_metric_specification,
                @view(null_testarr[:, :, sim]),
                thresholds[sim],
                ews_enddate_val
            )
            push!(null_ews_metrics, null_ews_metric)
        end
    end

    return ews_metrics, null_ews_metrics
end
