export create_cached_simulation_data

"""
    create_cached_simulation_data(scenario, data_arrs)

Pre-compute expensive simulation data once per scenario to avoid repeated computation
during parameter optimization. This includes noise arrays, test arrays, and EWS metrics.
"""
function create_cached_simulation_data(
        scenario::OptimizationScenario,
        data_arrs::NamedTuple
    )
    @unpack noise_specification, test_specification, percent_tested,
        ews_enddate_type, ews_metric_specification = scenario

    @unpack emergent_incidence_arr, null_incidence_arr,
        ensemble_specification, ensemble_single_Reff_thresholds_vec,
        emergent_outbreak_threshold_vecs = data_arrs

    # Create noise array once per scenario (expensive operation)
    noisearr = create_noise_arr(
        noise_specification,
        emergent_incidence_arr,
        ensemble_specification;
        seed = 1234,
    )[1]

    # Create test arrays once per scenario (expensive operation)
    testarr = create_testing_arrs(
        emergent_incidence_arr,
        noisearr,
        percent_tested,
        test_specification,
    )

    null_testarr = create_testing_arrs(
        null_incidence_arr,
        noisearr,
        percent_tested,
        test_specification,
    )

    # Must be the same size as need to pair each null series to a
    # emergent series when calculating the end date for the EWS
    # metric calculation
    @assert size(testarr, 3) == size(null_testarr, 3)

    # Select appropriate thresholds
    thresholds = get_enddate_thresholds(data_arrs, ews_enddate_type)

    # Pre-compute EWS metrics for all simulations
    ensemble_nsims = size(emergent_incidence_arr, 3)
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

    # Ensure we have at least some valid simulations
    if isempty(ews_metrics)
        error("No valid EWS metrics could be computed for any simulation")
    end

    return CachedSimulationData(
        testarr,
        null_testarr,
        thresholds,
        ews_metrics,
        null_ews_metrics
    )
end
