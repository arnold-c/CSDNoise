function plot_all_threshold_comparisons(percent_clinic_tested, base_parameters)
    ensemble_chars_vec = collect_threshold_char_vec(
        percent_clinic_tested, base_parameters
    )

    @unpack ensemble_specification,
    noise_specification, alertmethod,
    moving_avg_detection_lag =
        base_parameters

    baseplotdirpath = joinpath(
        plotsdir("ensemble/testing-comparison"),
        "R0_$(ensemble_specification.dynamics_parameters.R_0)",
        getdirpath(noise_specification),
        alertmethod,
    )

    if alertmethod != "dailythreshold"
        baseplotdirpath = joinpath(
            baseplotdirpath,
            "moveavglag_$(moving_avg_detection_lag)"
        )
    end

    clinic_tested_dir = joinpath(
        baseplotdirpath,
        "clinic-tested_$(percent_clinic_tested)"
    )

    accuracy_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_accuracy_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :accuracy,
                color = (ACCURACY_COLOR, 0.7),
            ),
        ],
        percent_clinic_tested;
        legend = false,
        xlabel = "Accuracy",
        columnfacetchar_label = "Alert Threshold",
        bins = 0.0:0.01:1.01,
        plotname = accuracy_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Accuracy plot saved"

    sens_spec_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_sens-spec_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :daily_sensitivity,
                label = "Sensitivity",
                color = (DAILY_SENSITIVITY_COLOR, 0.7),
            ),
            (
                char = :daily_specificity,
                label = "Specificity",
                color = (DAILY_SPECIFICITY_COLOR, 0.7),
            ),
        ],
        percent_clinic_tested;
        columnfacetchar_label = "Alert Threshold",
        bins = 0.0:0.01:1.01,
        plotname = sens_spec_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Sensitivity and specificity plot saved"

    ppv_npv_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_ppv-npv_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :daily_ppv,
                label = "PPV",
                color = (DAILY_PPV_COLOR, 0.7),
            ),
            (
                char = :daily_npv,
                label = "NPV",
                color = (DAILY_NPV_COLOR, 0.7),
            ),
        ],
        percent_clinic_tested;
        columnfacetchar_label = "Alert Threshold",
        bins = 0.0:0.01:1.01,
        plotname = ppv_npv_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "PPV and NPV plot saved"

    detectiondelay_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_detectiondelay_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :detectiondelays,
            color = (DETECTION_DELAY_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "Detection Delay (days)",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 5.0,
        meanlines = true,
        legend = false,
        plotname = detectiondelay_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Detection delay plot saved"

    nalertsperoutbreak_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_nalertsperoutbreak_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :n_alerts_per_outbreak,
            color = (N_ALERTS_PER_OUTBREAK_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "Alerts per Outbreak",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
        plotname = nalertsperoutbreak_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of alerts per outbreak plot saved"

    nfalsealerts_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_nfalsealerts_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :n_false_alerts,
            color = (N_FALSE_ALERTS_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "# False Alerts",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
        plotname = nfalsealerts_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of false alerts plot saved"

    nalerts_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_nalerts_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :nalerts,
            color = (N_ALERTS_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "# Alerts",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 10.0,
        meanlines = true,
        legend = false,
        plotname = nalerts_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of alerts plot saved"

    noutbreaks_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_noutbreaks_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :noutbreaks,
            color = (N_OUTBREAKS_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "# Outbreaks",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
        plotname = noutbreaks_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of outbreaks plot saved"

    nmissedoutbreaks_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_nmissedoutbreaks_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :n_missed_outbreaks,
            color = (N_MISSED_OUTBREAKS_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "# Missed Outbreaks",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 1.0,
        meanlines = true,
        legend = false,
        plotname = nmissedoutbreaks_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of missed outbreaks plot saved"

    size_outbreaks_detectmissed_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_size-outbreaks-detected-missed_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :detected_outbreak_size,
                label = "Size of Outbreaks Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7),
            ),
            (
                char = :missed_outbreak_size,
                label = "Size of Outbreaks Missed",
                color = (PERC_OUTBREAKS_MISSED_COLOR, 0.7),
            ),
        ],
        percent_clinic_tested;
        columnfacetchar_label = "Alert Threshold",
        binwidth = 200,
        normalization = :pdf,
        plotname = size_outbreaks_detectmissed_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Size of outbreaks detected and missed plot saved"

    perc_detectmissed_outbreak_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_percent-outbreaks-detected-missed_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :perc_true_outbreaks_detected,
                label = "Percent Outbreaks Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7),
            ),
            (
                char = :perc_true_outbreaks_missed,
                label = "Percent Outbreaks Missed",
                color = (PERC_OUTBREAKS_MISSED_COLOR, 0.7)),
        ],
        percent_clinic_tested;
        columnfacetchar_label = "Alert Threshold",
        binwidth = 0.02,
        plotname = perc_detectmissed_outbreak_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Percent outbreaks detected and missed plot saved"

    for i in eachindex(ensemble_chars_vec)
        perc_alerts_sum = sum(
            ensemble_chars_vec[i].OT_chars.perc_true_outbreaks_detected .+
            ensemble_chars_vec[i].OT_chars.perc_true_outbreaks_missed .- 1.0,
        )
        if perc_alerts_sum != 0.0
            @error "Warning. Sum of perc_alerts_false + perc_alerts_correct != 1.0, for i = $i. perc_alerts_sum = $perc_alerts_sum"
        end
    end

    perc_alerts_correctfalse_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_percent-alerts-correct-false_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :perc_alerts_correct,
                label = "Percent Alerts\nThat Are Correct",
                color = (PERC_ALERTS_CORRECT_COLOR, 0.7),
            ),
            (
                char = :perc_alerts_false,
                label = "Percent Alerts\nThat Are False",
                color = (PERC_ALERTS_FALSE_COLOR, 0.7)),
        ],
        percent_clinic_tested;
        columnfacetchar_label = "Alert Threshold",
        bins = -0.01:0.02:1.01,
        plotname = perc_alerts_correctfalse_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Percent alerts correct and false plot saved"

    for i in eachindex(ensemble_chars_vec)
        perc_alerts_sum = sum(
            ensemble_chars_vec[i].OT_chars.perc_alerts_false .+
            ensemble_chars_vec[i].OT_chars.perc_alerts_correct .- 1.0,
        )
        if perc_alerts_sum != 0.0
            @error "Warning. Sum of perc_alerts_false + perc_alerts_correct != 1.0, for i = $i. perc_alerts_sum = $perc_alerts_sum,\nTimes no outbreaks detected = $(length(findall(==(0), ensemble_chars_vec[i].OT_chars.nalerts)))"
            nan_perc_alerts_sum = NaNMath.sum(
                ensemble_chars_vec[i].OT_chars.perc_alerts_false .+
                ensemble_chars_vec[i].OT_chars.perc_alerts_correct .- 1.0,
            )
            if nan_perc_alerts_sum != 0.0
                @error "Ignoring NaN values in the percentage of alerts doesn't correct the issue"
                continue
            end
            @info "Ignoring NaN values in the percentage of alerts does correct the issue"
        end
    end

    perc_alertscorrect_outbreaksdetected_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_percent-alerts-correct-outbreaks-detected_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
            (
                char = :perc_alerts_correct,
                label = "Percent Alerts\nThat Are Correct",
                color = (PERC_ALERTS_CORRECT_COLOR, 0.7),
            ),
            (
                char = :perc_true_outbreaks_detected,
                label = "Percent Outbreaks\nThat Are Detected",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 0.7)),
        ],
        percent_clinic_tested;
        columnfacetchar_label = "Alert Threshold",
        bins = -0.01:0.02:1.01,
        plotname = perc_alertscorrect_outbreaksdetected_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Percent alerts correct and outbreaks detected plot saved"

    navoidablecases_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_n-avoidable-cases_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :avoidable_cases,
            color = (PERC_ALERTS_CORRECT_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "Number of Avoidable Cases",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 50.0,
        meanlines = true,
        legend = false,
        plotname = navoidablecases_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of avoidable cases plot saved"

    perc_casesbeforealerts_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_percent-cases-before-alerts_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :cases_perc_before_alerts,
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "Percentage of Outbreak\nBefore Alerts",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 0.01,
        meanlines = true,
        legend = false,
        plotname = perc_casesbeforealerts_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Percentage of cases before alerts plot saved"

    nunavoidablecases_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_n-unavoidable-cases_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :unavoidable_cases,
            color = (PERC_OUTBREAKS_MISSED_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "Number of Unavoidable Cases",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 50.0,
        meanlines = true,
        legend = false,
        plotname = nunavoidablecases_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Number of unavoidable cases plot saved"

    perc_casesafteralerts_plotname = "compare-outbreak_clinic-tested-$(percent_clinic_tested)_percent-cases-after-alerts_plot"
    save_compare_ensemble_OTchars_plot(
        ensemble_chars_vec,
        :alert_threshold,
        [
        (
            char = :cases_perc_after_alerts,
            color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
        )
    ],
        percent_clinic_tested;
        xlabel = "Percentage of Outbreak\nAfter Alerts",
        columnfacetchar_label = "Alert Threshold",
        binwidth = 0.01,
        meanlines = true,
        legend = false,
        plotname = perc_casesafteralerts_plotname,
        plotdirpath = clinic_tested_dir,
    )
    @info "Percentage of cases after alerts plot saved"

    @info "Finished saving all plots for % clinic tested = $percent_clinic_tested"
    @info "=============================================="
    println()

    return nothing
end

function collect_threshold_char_vec(percent_clinic_tested, base_parameters)
    @unpack test_spec_vec,
    alertthreshold_vec,
    ensemble_specification,
    noise_specification,
    outbreak_specification,
    moving_avg_detection_lag,
    percent_visit_clinic,
    alertmethod = base_parameters

    non_clinical_case_test_spec_vec = filter(
        spec -> spec != CLINICAL_CASE_TEST_SPEC,
        test_spec_vec
    )

    ensemble_scenario_spec_vec = Vector{ScenarioSpecification}(
        undef,
        length(non_clinical_case_test_spec_vec) * length(alertthreshold_vec) +
        length(alertthreshold_vec),
    )
    ensemble_chars_vec = Vector(
        undef, length(ensemble_scenario_spec_vec)
    )

    outbreak_detect_spec_vec = map(
        threshold -> OutbreakDetectionSpecification(
            threshold,
            moving_avg_detection_lag,
            percent_visit_clinic,
            percent_clinic_tested,
            alertmethod,
        ),
        alertthreshold_vec,
    )
    clinical_case_outbreak_detect_spec_vec = map(
        threshold -> OutbreakDetectionSpecification(
            threshold,
            moving_avg_detection_lag,
            percent_visit_clinic,
            1.0,
            alertmethod,
        ),
        alertthreshold_vec,
    )

    for (i, (ind_test_spec, outbreak_detect_spec)) in enumerate(
        Iterators.product(
            non_clinical_case_test_spec_vec, outbreak_detect_spec_vec
        ),
    )
        ensemble_scenario_spec = ScenarioSpecification(
            ensemble_specification,
            outbreak_specification,
            noise_specification,
            outbreak_detect_spec,
            ind_test_spec,
        )

        ensemble_scenario_spec_vec[i] = ensemble_scenario_spec
    end

    ensemble_scenario_spec_vec[(end - length(clinical_case_outbreak_detect_spec_vec) + 1):end] .= create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_specification],
            [outbreak_specification],
            [noise_specification],
            clinical_case_outbreak_detect_spec_vec,
            # TODO: update this to calculate for all detection thresholds
            [CLINICAL_CASE_TEST_SPEC],
        ),
    )

    @floop for (i, ensemble_scenario_spec) in pairs(ensemble_scenario_spec_vec)
        ensemble_chars_file = get_ensemble_file(ensemble_scenario_spec)

        ensemble_chars_vec[i] = (
            OT_chars = ensemble_chars_file["OT_chars"],
            outbreak_detect_spec = ensemble_scenario_spec.outbreak_detection_specification,
            ind_test_spec = ensemble_scenario_spec.individual_test_specification,
            noise_specification = noise_specification,
        )
    end

    return sort!(
        ensemble_chars_vec;
        by = x -> (
            x.outbreak_detect_spec.alert_threshold,
            x.ind_test_spec.specificity,
            x.ind_test_spec.test_result_lag,
        ),
    )
end
