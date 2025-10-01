export create_testing_arrs,
    create_testing_arrs!

function create_testing_arrs(
        incarr,
        noisearr,
        perc_tested,
        individual_test_spec::IndividualTestSpecification,
    )
    testarr = zeros(Int64, size(incarr, 1), 5, size(incarr, 3))

    create_testing_arrs!(
        testarr,
        incarr,
        noisearr,
        perc_tested,
        individual_test_spec.test_result_lag,
        individual_test_spec.sensitivity,
        individual_test_spec.specificity,
    )

    return testarr
end

function create_testing_arrs!(
        testarr,
        incarr,
        noisearr,
        perc_tested,
        testlag,
        testsens,
        testspec,
    )
    tlength = size(testarr, 1)

    for sim in axes(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(
            @view(testarr[:, 1, sim]), @view(incarr[:, 1, sim]), perc_tested
        )

        # Number of noise individuals tested
        calculate_tested!(
            @view(testarr[:, 2, sim]), @view(noisearr[:, sim]), perc_tested
        )

        # Number of test positive INFECTED individuals
        calculate_true_positives!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            tlength,
            testlag,
            testsens,
        )

        # Number of test positive NOISE individuals
        calculate_noise_positives!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            tlength,
            testlag,
            testspec,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])
    end

    return nothing
end

function create_testing_arrs(
        incarr,
        noisearr,
        outbreak_detect_spec::OutbreakDetectionSpecification,
        individual_test_spec::IndividualTestSpecification,
        time_specification::SimTimeParameters,
        ews_specification::EWSMetricSpecification,
    )
    testarr = zeros(Int64, size(incarr, 1), 5, size(incarr, 3))
    ewsvec = Vector{EWSMetrics}(undef, size(incarr, 3))
    test_movingavg_arr = zeros(Int64, size(incarr, 1), size(incarr, 3))
    inferred_positives_arr = zeros(Float64, size(incarr, 1), size(incarr, 3))
    ntested_worker_vec = Vector{Int64}(undef, size(incarr, 1))
    test_positivity_rate_worker_vec = Vector{Float64}(undef, size(incarr, 1))

    create_testing_arrs!(
        testarr,
        ewsvec,
        test_movingavg_arr,
        inferred_positives_arr,
        ntested_worker_vec,
        test_positivity_rate_worker_vec,
        incarr,
        noisearr,
        outbreak_detect_spec.alert_method.method_name,
        outbreak_detect_spec.alert_threshold,
        outbreak_detect_spec.moving_average_lag,
        outbreak_detect_spec.percent_tested,
        individual_test_spec.test_result_lag,
        individual_test_spec.sensitivity,
        individual_test_spec.specificity,
        time_specification.tstep,
        ews_specification,
    )

    return testarr, ewsvec, test_movingavg_arr, inferred_positives_arr
end

function create_testing_arrs!(
        testarr,
        ewsvec,
        test_movingavg_arr,
        inferred_positives_arr,
        ntested_worker_vec,
        test_positivity_rate_worker_vec,
        incarr,
        noisearr,
        alert_method,
        alertthreshold,
        moveavglag,
        perc_tested,
        testlag,
        testsens,
        testspec,
        tstep,
        ews_specification,
    )
    tlength = size(testarr, 1)

    for sim in axes(incarr, 3)
        # Number of infectious individuals tested
        calculate_tested!(
            @view(testarr[:, 1, sim]), @view(incarr[:, 1, sim]), perc_tested
        )

        # Number of noise individuals tested
        calculate_tested!(
            @view(testarr[:, 2, sim]), @view(noisearr[:, sim]), perc_tested
        )

        # Number of TOTAL individuals tested
        @. @views ntested_worker_vec .=
            testarr[:, 1, sim] + testarr[:, 2, sim]

        # Number of test positive INFECTED individuals
        calculate_true_positives!(
            @view(testarr[:, 3, sim]),
            @view(testarr[:, 1, sim]),
            tlength,
            testlag,
            testsens,
        )

        # Number of test positive NOISE individuals
        calculate_noise_positives!(
            @view(testarr[:, 4, sim]),
            @view(testarr[:, 2, sim]),
            tlength,
            testlag,
            testspec,
        )

        # Number of test positive TOTAL individuals
        @. testarr[:, 5, sim] =
            @view(testarr[:, 3, sim]) + @view(testarr[:, 4, sim])

        # Calculate moving average of TOTAL test positives
        calculate_movingavg!(
            @view(test_movingavg_arr[:, sim]),
            @view(testarr[:, 5, sim]),
            moveavglag,
        )

        calculate_test_positivity_rate!(
            test_positivity_rate_worker_vec,
            @view(testarr[:, 5, sim]),
            ntested_worker_vec,
            moveavglag,
        )

        infer_true_positives!(
            @view(inferred_positives_arr[:, sim]),
            test_positivity_rate_worker_vec,
            @view(testarr[:, 5, sim]),
            ntested_worker_vec,
            @view(incarr[:, 1, sim]) + @view(noisearr[:, sim]),
            testlag,
        )

        ewsvec[sim] = if alert_method == "movingavg"
            EWSMetrics(
                ews_specification,
                @view(test_movingavg_arr[:, sim]),
            )
        elseif alert_method == "dailythreshold_movingavg"
            EWSMetrics(
                ews_specification,
                @view(test_movingavg_arr[:, sim]),
            )
        elseif alert_method == "inferred_movingavg"
            EWSMetrics(
                ews_specification,
                @view(inferred_positives_arr[:, sim]),
            )
        else
            error("alert_method must be one of \"movingavg\", \"dailythreshold_movingavg\", \"inferred_movingavg\"")
        end
    end

    return nothing
end
