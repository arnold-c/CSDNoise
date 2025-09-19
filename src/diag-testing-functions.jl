# module DiagTestingFunctions
#
# export create_testing_arr, create_testing_arr!, calculate_tested!,
#     calculate_pos, calculate_pos!, calculate_movingavg, calculate_movingavg!,
#     detectoutbreak, detectoutbreak!, calculate_ot_characterstics,
#     calculate_noutbreaks, calculate_OutbreakThresholdChars,
#     run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation

using DrWatson
using StatsBase
using FLoops
using NaNMath: NaNMath

# include("detection-thresholds.jl")
# # using .DetectionThresholds
#
# include("structs.jl")
# using .ODStructs

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

function calculate_tested!(outvec, invec, perc_tested)
    return @. outvec = round(invec * perc_tested)
end

function calculate_noise_positives!(outvec, tested_vec, tlength, lag, spec)
    tested_multiplier = 1.0 - spec
    calculate_positives!(outvec, tested_vec, tlength, lag, tested_multiplier)
    return nothing
end

function calculate_true_positives!(outvec, tested_vec, tlength, lag, sens)
    calculate_positives!(outvec, tested_vec, tlength, lag, sens)
    return nothing
end

function calculate_positives!(
        npos_vec, tested_vec, tlength, lag, tested_multiplier
    )
    @inbounds for day in eachindex(tested_vec)
        if day + lag <= tlength
            npos_vec[day + lag] = Int64(
                round(tested_vec[day] * tested_multiplier)
            )
        end
    end
    return nothing
end

function calculate_movingavg(invec, avglag)
    outvec = zeros(Float64, size(invec, 1))

    calculate_movingavg!(outvec, invec, avglag)

    return outvec
end

function calculate_movingavg(
        invec::T1, avglag
    ) where {T1 <: AbstractArray{Integer}}
    outvec = zeros(eltype(invec), size(invec, 1))

    calculate_movingavg!(outvec, invec, avglag)

    return outvec
end

function calculate_movingavg!(outvec, invec, avglag)
    if avglag == 0
        outvec .= invec
        return nothing
    end

    @inbounds for day in eachindex(invec)
        outvec[day] = calculate_float_daily_movingavg(invec, day, avglag)
    end
    return nothing
end

function calculate_movingavg!(
        outvec::T1, invec, avglag
    ) where {T1 <: AbstractArray{<:Integer}}
    if avglag == 0
        outvec .= invec
        return nothing
    end

    @inbounds for day in eachindex(invec)
        outvec[day] = calculate_int_daily_movingavg(invec, day, avglag)
    end
    return nothing
end

function calculate_float_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = calculate_daily_movingavg_startday(day, avglag)
    return mean(@view(invec[moveavg_daystart:day]))
end

function calculate_int_daily_movingavg(invec, day, avglag)
    @inline moveavg_daystart = calculate_daily_movingavg_startday(day, avglag)
    return Int64(round(mean(@view(invec[moveavg_daystart:day]))))
end

function calculate_daily_movingavg_startday(day, avglag)
    if day < avglag
        moveavg_daystart = 1
    else
        moveavg_daystart = day - avglag + 1
    end
    return moveavg_daystart
end

function infer_true_positives(
        test_positivity_rate_vec,
        test_positive_vec,
        total_test_vec,
        total_clinic_visit_vec,
        test_result_lag,
    )
    inferred_positives = zeros(Float64, length(total_clinic_visit_vec))

    if test_result_lag == 0
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            test_positive_vec,
            total_test_vec,
            total_clinic_visit_vec,
        )
    else
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            total_clinic_visit_vec,
        )
    end

    return inferred_positives
end

function infer_true_positives!(
        inferred_positives,
        test_positivity_rate_vec,
        test_positive_vec,
        total_test_vec,
        total_clinic_visit_vec,
        test_result_lag,
    )
    if test_result_lag == 0
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            test_positive_vec,
            total_test_vec,
            total_clinic_visit_vec,
        )
    else
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            total_clinic_visit_vec,
        )
    end

    return nothing
end
function infer_true_positives!(
        inferred_positives,
        test_positivity_rate_vec,
        test_positive_vec,
        total_test_vec,
        total_clinic_visit_vec,
    )
    @. inferred_positives =
        test_positive_vec +
        (total_clinic_visit_vec - total_test_vec) * test_positivity_rate_vec

    return nothing
end

function infer_true_positives!(
        inferred_positives,
        test_positivity_rate_vec,
        total_clinic_visit_vec,
    )
    return @. inferred_positives = total_clinic_visit_vec * test_positivity_rate_vec
end

function calculate_test_positivity_rate(
        test_positive_vec, total_test_vec, agg_days
    )
    test_positivity_rate_vec = zeros(Float64, length(total_test_vec))

    calculate_test_positivity_rate!(
        test_positivity_rate_vec, test_positive_vec, total_test_vec, agg_days
    )

    return test_positivity_rate_vec
end

function calculate_test_positivity_rate!(
        test_positivity_rate_vec, test_positive_vec, total_test_vec, agg_days
    )
    for i in axes(test_positivity_rate_vec, 1)
        if i < agg_days
            test_positivity_rate_vec[i] =
                sum(test_positive_vec[begin:i]) / sum(total_test_vec[begin:i])
            continue
        end

        test_positivity_rate_vec[i] =
            sum(test_positive_vec[(i - agg_days + 1):i]) /
            sum(total_test_vec[(i - agg_days + 1):i])
    end
    return
end

function calculate_test_positivity(
        true_positive_vec, total_test_vec, alert_vec, agg_days
    )
    @views outvec = zeros(Float64, length(true_positive_vec) ÷ agg_days, 2)
    @inbounds for i in axes(outvec, 1)
        start_ind = 1 + (i - 1) * agg_days
        end_ind = start_ind + (agg_days - 1)

        @views total_test_sum = sum(total_test_vec[start_ind:end_ind])
        @views true_positive_sum = sum(true_positive_vec[start_ind:end_ind])
        @views num_outbreak_days = sum(alert_vec[start_ind:end_ind])
        agg_outbreak_status = num_outbreak_days >= agg_days / 2 ? 1 : 0

        outvec[i, 1] = true_positive_sum / total_test_sum
        outvec[i, 2] = agg_outbreak_status
    end
    return outvec
end

# end
