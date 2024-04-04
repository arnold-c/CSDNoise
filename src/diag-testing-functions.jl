# module DiagTestingFunctions
#
# export create_testing_arr, create_testing_arr!, calculate_tested!,
#     calculate_pos, calculate_pos!, calculate_movingavg, calculate_movingavg!,
#     detectoutbreak, detectoutbreak!, calculate_ot_characterstics,
#     calculate_noutbreaks, calculate_OutbreakThresholdChars,
#     run_OutbreakThresholdChars_creation, OutbreakThresholdChars_creation

using DrWatson
using Statistics
using FreqTables
using ThreadsX
using FLoops
using StructArrays
using NaNMath: NaNMath
using TestItems

# include("detection-thresholds.jl")
# # using .DetectionThresholds
#
# include("structs.jl")
# using .ODStructs

function create_testing_arrs(
    incarr,
    noisearr,
    outbreak_detect_spec::OutbreakDetectionSpecification,
    individual_test_spec::IndividualTestSpecification,
    time_specification::SimTimeParameters,
    ews_specification::EWSMetricSpecification,
)
    testarr = zeros(Int64, size(incarr, 1), 7, size(incarr, 3))
    ewsarr = Array{EWSMetrics}(undef, size(incarr, 1), size(incarr, 3))
    test_movingavg_arr = zeros(Int64, size(incarr, 1), size(incarr, 3))
    inferred_positives_arr = zeros(Float64, size(incarr, 1), size(incarr, 3))
    ntested_worker_vec = Vector{Int64}(undef, size(incarr, 1))
    test_positivity_rate_worker_vec = Vector{Float64}(undef, size(incarr, 1))

    create_testing_arrs!(
        testarr,
        ewsarr,
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

    return testarr, ewsarr, test_movingavg_arr, inferred_positives_arr
end

function create_testing_arrs!(
    testarr,
    ewsarr,
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
        @turbo @. @views ntested_worker_vec .=
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

        detectoutbreak_args = @match alert_method begin
            "movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(test_movingavg_arr[:, sim]),
                alertthreshold,
            )
            "dailythreshold_movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(testarr[:, 5, sim]),
                @view(test_movingavg_arr[:, sim]),
                alertthreshold,
            )
            "inferred_movingavg" => (
                @view(testarr[:, 6, sim]),
                @view(inferred_positives_arr[:, sim]),
                alertthreshold,
            )
        end

        # TOTAL Test positive individuals trigger outbreak response
        detectoutbreak!(detectoutbreak_args...)

        ewsarr[:, sim] .= @match alert_method begin
            "movingavg" => EWSMetrics(
                ews_specification,
                @view(test_movingavg_arr[:, sim]),
                tstep
            )

            "dailythreshold_movingavg" => EWSMetrics(
                ews_specification,
                @view(test_movingavg_arr[:, sim]),
                tstep
            )

            "inferred_movingavg" => EWSMetrics(
                ews_specification,
                @view(inferred_positives_arr[:, sim]),
                tstep,
            )
        end

        # Triggered outbreak equal to actual outbreak status
        @. testarr[:, 7, sim] =
            @view(testarr[:, 6, sim]) == @view(incarr[:, 3, sim])
    end

    return nothing
end

function calculate_tested!(outvec, invec, perc_tested)
    @. outvec = round(invec * perc_tested)
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

@testitem "Moving average" begin
    using CSDNoise, Statistics

    daily_testpositives = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    movingavg_testpositives = calculate_movingavg(
        daily_testpositives,
        5,
    )

    @test isequal(
        movingavg_testpositives,
        [
            mean([1]),
            mean([1, 2]),
            mean([1, 2, 3]),
            mean([1, 2, 3, 4]),
            mean([1, 2, 3, 4, 5]),
            mean([2, 3, 4, 5, 6]),
            mean([3, 4, 5, 6, 7]),
            mean([4, 5, 6, 7, 8]),
            mean([5, 6, 7, 8, 9]),
            mean([6, 7, 8, 9, 10]),
        ],
    )

    @test length(movingavg_testpositives) == length(daily_testpositives)

    daily_testpositives = zeros(Int64, 10, 2)
    daily_testpositives[:, 1] .= [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    calculate_movingavg!(
        @view(daily_testpositives[:, 2]),
        @view(daily_testpositives[:, 1]),
        5
    )

    @test isequal(
        daily_testpositives[:, 2],
        Int64.(
            round.([
                mean([1]),
                mean([1, 2]),
                mean([1, 2, 3]),
                mean([1, 2, 3, 4]),
                mean([1, 2, 3, 4, 5]),
                mean([2, 3, 4, 5, 6]),
                mean([3, 4, 5, 6, 7]),
                mean([4, 5, 6, 7, 8]),
                mean([5, 6, 7, 8, 9]),
                mean([6, 7, 8, 9, 10]),
            ])
        ),
    )
end

function calculate_movingavg(invec, avglag)
    outvec = zeros(Float64, size(invec, 1))

    calculate_movingavg!(outvec, invec, avglag)

    return outvec
end

function calculate_movingavg(
    invec::T1, avglag
) where {T1<:AbstractArray{Integer}}
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
) where {T1<:AbstractArray{<:Integer}}
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

@testitem "Detect outbreak" begin
    using CSDNoise

    incvec = [1, 3, 10, 15, 20, 3, 1]
    avgvec = calculate_movingavg(incvec, 3)
    threshold = 5

    outbreakvec = detectoutbreak(incvec, threshold)
    @test isequal(outbreakvec, [0, 0, 1, 1, 1, 0, 0])

    outbreakvec = detectoutbreak(avgvec, threshold)
    @test isequal(outbreakvec, [0, 0, 0, 1, 1, 1, 1])

    outbreakvec = detectoutbreak(incvec, avgvec, threshold)
    @test isequal(outbreakvec, [0, 0, 1, 1, 1, 1, 1])
end

function detectoutbreak(incvec, avgvec, threshold)
    outbreak = zeros(Int64, length(incvec))

    detectoutbreak!(outbreak, incvec, avgvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, incvec, avgvec, threshold)

Determines whether an outbreak has been detected based on an infection daily and moving average timeseries and a threshold.
`incvec` should be the daily incidence timeseries, and `avgvec` should be the moving average of the daily incidence timeseries.
"""
function detectoutbreak!(outbreakvec, incvec, avgvec, threshold)
    @. outbreakvec = ifelse(incvec >= threshold || avgvec >= threshold, 1, 0)

    return nothing
end

function detectoutbreak(infectionsvec, threshold)
    outbreak = zeros(Int64, length(infectionsvec))

    detectoutbreak!(outbreak, infectionsvec, threshold)

    return outbreak
end

"""
    detectoutbreak!(outbreakvec, infectionsvec, threshold)

Determines whether an outbreak has been detected based on an infection timeseries and a threshold.
The `infectionsvec` can either be a vector of daily infections or a vector of the moving average of daily infections.
"""
function detectoutbreak!(outbreakvec, infectionsvec, threshold)
    @. outbreakvec = ifelse(infectionsvec >= threshold, 1, 0)

    return nothing
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
            total_clinic_visit_vec
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
            total_clinic_visit_vec
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
    total_clinic_visit_vec
)
    @. inferred_positives = total_clinic_visit_vec * test_positivity_rate_vec
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
end

@testitem "Inferred cases" begin
    using CSDNoise, Statistics

    daily_testpositives = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    daily_tested = repeat([10], length(daily_testpositives))

    movingavg_testpositives = calculate_test_positivity_rate(
        daily_testpositives,
        daily_tested,
        5
    )

    @test isequal(
        movingavg_testpositives,
        [
            mean([1]) / 10,
            mean([1, 2]) / 10,
            mean([1, 2, 3]) / 10,
            mean([1, 2, 3, 4]) / 10,
            mean([1, 2, 3, 4, 5]) / 10,
            mean([2, 3, 4, 5, 6]) / 10,
            mean([3, 4, 5, 6, 7]) / 10,
            mean([4, 5, 6, 7, 8]) / 10,
            mean([5, 6, 7, 8, 9]) / 10,
            mean([6, 7, 8, 9, 10]) / 10,
        ],
    )

    @test length(movingavg_testpositives) == length(daily_testpositives)

    daily_clinic_visits = repeat([100], length(daily_testpositives))

    inferred_positives_no_lag = infer_true_positives(
        movingavg_testpositives,
        daily_testpositives,
        daily_tested,
        daily_clinic_visits,
        0,
    )

    @test isequal(
        inferred_positives_no_lag,
        daily_testpositives +
        movingavg_testpositives .* (daily_clinic_visits - daily_tested),
    )

    inferred_positives_3d_lag = infer_true_positives(
        movingavg_testpositives,
        daily_testpositives,
        daily_tested,
        daily_clinic_visits,
        3,
    )

    @test isequal(
        inferred_positives_3d_lag,
        movingavg_testpositives .* daily_clinic_visits,
    )
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

function calculate_OutbreakThresholdChars(
    testarr, infecarr, thresholds_vec, noise_rubella_prop
)
    OT_chars = map(axes(infecarr, 3)) do sim
        dailychars = calculate_daily_detection_characteristics(
            @view(testarr[:, 6, sim]), @view(infecarr[:, 3, sim])
        )
        alertrle = rle(@view(testarr[:, 6, sim]))
        outbreakbounds = thresholds_vec[sim][
            (@view(thresholds_vec[sim][:, 4]) .== 1), 1:3
        ]
        alertbounds = calculate_outbreak_thresholds(alertrle; ncols = 2)

        detectionchars = calculate_outbreak_detection_characteristics(
            outbreakbounds, alertbounds
        )

        first_matched_bounds = filter_first_matched_bounds(
            detectionchars.matched_bounds
        )
        delay_vec = calculate_delay_vec(first_matched_bounds)

        cases_before_alert_vec, caseperc_before_alert,
        cases_after_alert_vec, caseperc_after_alert_vec = calculate_cases_before_after_alert(
            @view(infecarr[:, 1, sim]), first_matched_bounds, delay_vec
        )

        unavoidable_cases = calculate_unavoidable_cases(
            detectionchars.missed_outbreak_size, cases_before_alert_vec
        )

        avoidable_cases = calculate_avoidable_cases(cases_after_alert_vec)

        n_outbreak_cases = unavoidable_cases + avoidable_cases

        n_tests = calculate_n_tests(
            @view(testarr[:, 1, sim]), @view(testarr[:, 2, sim])
        )

        OutbreakThresholdChars(
            dailychars...,
            detectionchars...,
            delay_vec,
            cases_before_alert_vec,
            caseperc_before_alert,
            cases_after_alert_vec,
            caseperc_after_alert_vec,
            unavoidable_cases,
            avoidable_cases,
            n_outbreak_cases,
            n_tests,
            NaNMath.mean(noise_rubella_prop),
        )
    end

    return StructArray(OT_chars)
end

function calculate_n_tests(infectious_tested_vec, noise_tested_vec)
    return sum(infectious_tested_vec) + sum(noise_tested_vec)
end

function calculate_unavoidable_cases(
    missed_outbreak_vec, cases_before_alert_vec
)
    return sum(missed_outbreak_vec) + sum(cases_before_alert_vec)
end

function calculate_avoidable_cases(cases_after_alert_vec)
    return sum(cases_after_alert_vec)
end

function calculate_daily_detection_characteristics(testvec, infecvec)
    crosstab = freqtable(testvec, infecvec)

    tp = freqtable_error_default_zero(crosstab, 1, 1)
    fp = freqtable_error_default_zero(crosstab, 1, 0)
    tn = freqtable_error_default_zero(crosstab, 0, 0)
    fn = freqtable_error_default_zero(crosstab, 0, 1)

    sens = tp / (tp + fn)
    spec = tn / (tn + fp)

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    return sens, spec, ppv, npv
end

# TODO: write tests that check the function works when crosstab produces a 2x2 table,
# a 2x1, and a 1x1 table
function freqtable_error_default_zero() end

function freqtable_error_default_zero(
    freqtable, testing_class::T, actual_class::T
) where {T<:Integer}
    return try
        freqtable[Name(testing_class), Name(actual_class)]
    catch
        0
    end
end

function calculate_noutbreaks(outbreakrle)
    return length(findall(==(1), outbreakrle[1]))
end

function calculate_outbreak_detection_characteristics(
    outbreakbounds, alertbounds
)
    filtered_matched_bounds, periodssum_vec, alerts_per_outbreak_vec = match_outbreak_detection_bounds(
        outbreakbounds, alertbounds
    )

    noutbreaks = size(outbreakbounds, 1)
    nalerts = size(alertbounds, 1)

    detected_outbreak_size = periodssum_vec[(alerts_per_outbreak_vec .> 0)]
    missed_outbreak_size = periodssum_vec[(alerts_per_outbreak_vec .== 0)]

    n_true_outbreaks_detected = length(
        Set(@view(filtered_matched_bounds[:, 1]))
    )

    n_missed_outbreaks = noutbreaks - n_true_outbreaks_detected

    n_correct_alerts = size(filtered_matched_bounds, 1)

    n_false_alerts = sum(nalerts - n_correct_alerts)

    perc_true_outbreaks_detected = n_true_outbreaks_detected / noutbreaks
    perc_true_outbreaks_missed = n_missed_outbreaks / noutbreaks
    falsealert_trueoutbreak_prop = n_false_alerts / noutbreaks
    correctalert_trueoutbreak_prop = n_correct_alerts / noutbreaks # c.f. sensitivity

    trueoutbreak_alerts_prop = noutbreaks / nalerts
    outbreaksmissed_alerts_prop = n_missed_outbreaks / nalerts
    perc_alerts_false = n_false_alerts / nalerts
    perc_alerts_correct = n_correct_alerts / nalerts # c.f. PPV

    accuracy = NaNMath.mean([perc_true_outbreaks_detected, perc_alerts_correct])

    return (
        accuracy = accuracy,
        matched_bounds = filtered_matched_bounds,
        noutbreaks = noutbreaks,
        nalerts = nalerts,
        detected_outbreak_size = detected_outbreak_size,
        missed_outbreak_size = missed_outbreak_size,
        n_true_outbreaks_detected = n_true_outbreaks_detected,
        n_missed_outbreaks = n_missed_outbreaks,
        n_correct_alerts = n_correct_alerts,
        n_false_alerts = n_false_alerts,
        alertsperoutbreak = alerts_per_outbreak_vec,
        periodsumvec = periodssum_vec,
        perc_true_outbreaks_detected = perc_true_outbreaks_detected,
        perc_true_outbreaks_missed = perc_true_outbreaks_missed,
        falsealert_trueoutbreak_prop = falsealert_trueoutbreak_prop,
        correctalert_trueoutbreak_prop = correctalert_trueoutbreak_prop,
        trueoutbreak_alerts_prop = trueoutbreak_alerts_prop,
        outbreaksmissed_alerts_prop = outbreaksmissed_alerts_prop,
        perc_alerts_false = perc_alerts_false,
        perc_alerts_correct = perc_alerts_correct,
    )
end

function match_outbreak_detection_bounds(outbreakbounds, alertbounds)
    all_matched_bounds = zeros(
        Int64, size(outbreakbounds, 1) + size(alertbounds, 1), 5
    )
    alerts_per_outbreak_vec = zeros(Int64, size(outbreakbounds, 1))
    periodssum_vec = zeros(Int64, size(outbreakbounds, 1))

    outbreak_number = 1
    alert_rownumber = 1
    for (outbreak_number, (outbreaklower, outbreakupper, periodsum)) in
        pairs(eachrow(outbreakbounds))
        periodssum_vec[outbreak_number] = periodsum
        for (alertlower, alertupper) in
            eachrow(@view(alertbounds[alert_rownumber:end, :]))
            if alertlower > outbreakupper
                break
            end
            if alertlower >= outbreaklower
                all_matched_bounds[alert_rownumber, :] .= outbreaklower,
                outbreakupper, alertlower, alertupper, periodsum

                alerts_per_outbreak_vec[outbreak_number] += 1

                alert_rownumber += 1
                continue
            end
            if alertlower <= outbreaklower &&
                alertupper > outbreaklower
                all_matched_bounds[alert_rownumber, :] .= outbreaklower,
                outbreakupper, alertlower, alertupper, periodsum

                alerts_per_outbreak_vec[outbreak_number] += 1

                alert_rownumber += 1
                continue
            end
        end
        outbreak_number += 1
    end

    filtered_matched_bounds = @view(
        all_matched_bounds[(all_matched_bounds[:, 2] .> 0), :]
    )
    return filtered_matched_bounds, periodssum_vec, alerts_per_outbreak_vec
end

function calculate_delay_vec(first_matchedbounds)
    return @views first_matchedbounds[:, 3] .- first_matchedbounds[:, 1]
end

@testitem "Matched bounds" begin
    using CSDNoise

    matched_bounds = [
        10 60 5 15 600
        10 60 17 40 600
        10 60 50 80 600
        100 180 90 105 700
        100 180 110 160 700
        380 410 390 420 900
        500 540 495 550 1000
    ]
    filtered_bounds = [
        10 60 5 15 600
        100 180 90 105 700
        380 410 390 420 900
        500 540 495 550 1000
    ]

    @test isequal(
        filter_first_matched_bounds(matched_bounds), filtered_bounds
    )
end

function filter_first_matched_bounds(matchedbounds)
    indices = calculate_first_matched_bounds_index(matchedbounds)
    return matchedbounds[indices, :]
end

function calculate_first_matched_bounds_index(matchedbounds)
    return map(
        outbreaklower -> findfirst(
            isequal(outbreaklower),
            @view(matchedbounds[:, 1])
        ),
        unique(@view(matchedbounds[:, 1])),
    )
end

function calculate_cases_before_after_alert(
    incvec, first_matchedbounds, delay_vec
)
    casesaftervec = zeros(Int64, length(delay_vec))
    casesafterpercvec = zeros(Float64, length(casesaftervec))
    casesbeforevec = zeros(Int64, length(delay_vec))
    casesbeforepercvec = zeros(Float64, length(casesaftervec))
    calculate_cases_before_after_alert!(
        casesbeforevec,
        casesbeforepercvec,
        casesaftervec,
        casesafterpercvec,
        incvec,
        first_matchedbounds,
        delay_vec,
    )
    return casesbeforevec, casesbeforepercvec, casesaftervec, casesafterpercvec
end

function calculate_cases_before_after_alert!(
    casesbeforevec,
    casesbeforepercvec,
    casesaftervec,
    casesafterpercvec,
    invec,
    first_matchedbounds,
    delay_vec,
)
    for (
        alertnumber,
        (outbreaklower, outbreakupper, alertlower, alertupper, periodsum),
    ) in
        pairs(eachrow(first_matchedbounds))
        if delay_vec[alertnumber] <= 0
            casesaftervec[alertnumber, 1] = periodsum
            casesbeforevec[alertnumber, 1] = 0
        else
            @views casesaftervec[alertnumber] = sum(
                invec[alertlower:outbreakupper]
            )
            casesbeforevec[alertnumber] = periodsum - casesaftervec[alertnumber]
        end
        casesbeforepercvec[alertnumber] =
            casesbeforevec[alertnumber] / periodsum
        casesafterpercvec[alertnumber] = casesaftervec[alertnumber] / periodsum
    end
    return nothing
end

@testitem "Cases before and after alert" begin
    using CSDNoise
    incarr = [
        repeat([1], 9)...,
        repeat([12], 51)...,
        repeat([1], 39)...,
        repeat([8], 81)...,
        repeat([1], 199)...,
        repeat([30], 31)...,
        repeat([1], 89)...,
        repeat([25], 41)...,
    ]

    matched_bounds = [
        10 60 5 15 612
        100 180 90 105 648
        380 410 390 420 930
        500 540 495 550 1025
    ]

    delay_vec = [-5, -10, 10, -5]

    @test isequal(
        calculate_cases_before_after_alert(
            incarr, matched_bounds, delay_vec
        ),
        (
            [0, 0, 300, 0],
            [0.0, 0.0, 300 / 930, 0.0],
            [612, 648, 630, 1025],
            [1.0, 1.0, 630 / 930, 1.0],
        ),
    )
end

@testitem "Outbreak detection characteristics" begin
    using CSDNoise
    using Statistics
    outbreakbounds = [
        2 4 500 100
        10 60 600 10
        100 180 700 2000
        300 340 800 20
        380 410 900 400
        500 540 1000 5000
    ]
    detectionbounds = [
        5 15
        17 40
        50 80
        90 105
        110 160
        390 420
        495 550
        590 595
    ]

    @test isequal(
        calculate_outbreak_detection_characteristics(
            outbreakbounds, detectionbounds
        ),
        (
            accuracy = mean([4 / 6, 7 / 8]),
            matched_bounds = [
                10 60 5 15 600
                10 60 17 40 600
                10 60 50 80 600
                100 180 90 105 700
                100 180 110 160 700
                380 410 390 420 900
                500 540 495 550 1000
            ],
            noutbreaks = 6,
            nalerts = 8,
            detected_outbreak_size = [600, 700, 900, 1000],
            missed_outbreak_size = [500, 800],
            n_true_outbreaks_detected = 4,
            n_missed_outbreaks = 2,
            n_correct_alerts = 7,
            n_false_alerts = 1,
            alertsperoutbreak = [0, 3, 2, 0, 1, 1],
            periodsumvec = [500, 600, 700, 800, 900, 1000],
            perc_true_outbreaks_detected = 4 / 6,
            perc_true_outbreaks_missed = 2 / 6,
            falsealert_trueoutbreak_prop = 1 / 6,
            correctalert_trueoutbreak_prop = 7 / 6,
            trueoutbreak_alerts_prop = 6 / 8,
            outbreaksmissed_alerts_prop = 2 / 8,
            perc_alerts_false = 1 / 8,
            perc_alerts_correct = 7 / 8,
        ),
    )
end

# end
