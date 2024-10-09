using StructArrays

function plot_all_single_scenarios(
    noisevec,
    noisedir,
    incvec,
    outbreak_status_vec,
    testvec,
    test_movingvg_vec,
    Reff_vec,
    Reff_thresholds,
    outbreak_thresholds,
    ewsvec,
    ewsdir,
    ews_enddate,
    test_specification,
    outbreak_detection_specification,
    time_specification;
    ews_metrics = [
        "autocorrelation",
        "autocovariance",
        "coefficient_of_variation",
        "index_of_dispersion",
        "kurtosis",
        "mean",
        "skewness",
        "variance",
    ],
    aggregation = 1,
    sim = 1,
    force = false,
)
    ensemble_noise_plotpath = joinpath(
        plotsdir(),
        "ensemble",
        "single-scenario",
        noisedir,
        "sim-$(sim)",
    )
    mkpath(ensemble_noise_plotpath)

    noise_plottitle = "Sens: $(test_specification.sensitivity), Spec: $(test_specification.specificity), Lag: $(test_specification.test_result_lag),\nThreshold: $(outbreak_detection_specification.alert_threshold), Perc Clinic Tested: $(outbreak_detection_specification.percent_clinic_tested)\nNoise: $(noisedir), Alert Method: $(outbreak_detection_specification.alert_method.method_name)"

    plotpath = joinpath(
        ensemble_noise_plotpath,
        "ensemble-sim_single-scenario_incidence-testing_sim-$(sim).png",
    )

    if !isfile(plotpath) || force
        ensemble_single_scenario_incidence_testing_plot = incidence_testing_plot(
            incvec,
            outbreak_status_vec,
            noisevec,
            testvec,
            test_movingvg_vec,
            time_specification;
            aggregation = aggregation,
            plottitle = noise_plottitle,
        )

        save(
            plotpath,
            ensemble_single_scenario_incidence_testing_plot,
        )

        Makie.empty!(ensemble_single_scenario_incidence_testing_plot)
    end

    ews_plotpath = joinpath(
        ensemble_noise_plotpath,
        ewsdir,
        "enddate-$(ews_enddate)",
    )

    mkpath(ews_plotpath)

    for ewsmetric in ews_metrics
        plotpath = joinpath(
            ews_plotpath,
            "ensemble-sim_single-scenario_$(ewsmetric)_ews_agg-$(aggregation)_sim-$(sim).png",
        )

        if !isfile(plotpath) || force
            ensemble_single_Reff_ews_plot = Reff_ews_plot(
                incvec,
                Reff_vec,
                Reff_thresholds,
                ewsvec,
                Symbol(ewsmetric),
                outbreak_thresholds,
                time_specification,
            )

            save(
                plotpath,
                ensemble_single_Reff_ews_plot,
            )
        end
    end

    return nothing
end
