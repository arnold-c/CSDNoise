function plot_all_single_scenarios(
    noisearr,
    poisson_noise_prop,
    noisedir,
    incarr,
    testarr,
    test_movingvg_arr,
    test_specification,
    outbreak_detection_specification,
    time_specification,
)
    ensemble_noise_plotpath = joinpath(
        plotsdir(),
        "ensemble",
        "single-scenario",
        noisedir,
    )
    mkpath(ensemble_noise_plotpath)

    ensemble_noise_fig = visualize_ensemble_noise(
        noisearr,
        poisson_noise_prop,
        time_specification,
        noisedir,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_noise.png",
        ),
        ensemble_noise_fig; size = (2200, 1600),
    )

    Makie.empty!(ensemble_noise_fig)

    noise_plottitle = "Sens: $(test_specification.sensitivity), Spec: $(test_specification.specificity), Lag: $(test_specification.test_result_lag),\nThreshold: $(outbreak_detection_specification.alert_threshold), Perc Clinic Tested: $(outbreak_detection_specification.percent_clinic_tested)\nNoise: $(noisedir), Alert Method: $(outbreak_detection_specification.alert_method.method_name)"

    ensemble_single_scenario_incidence_testing_plot = incidence_testing_plot(
        incarr,
        noisearr,
        testarr,
        test_movingvg_arr,
        outbreak_detection_specification,
        time_specification;
        sim = 1,
        plottitle = noise_plottitle,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_incidence-testing.png",
        ),
        ensemble_single_scenario_incidence_testing_plot,
    )

    Makie.empty!(ensemble_single_scenario_incidence_testing_plot)

    ensemble_single_scenario_testing_timeseries_plot = testing_plot(
        testarr,
        time_specification;
        plottitle = noise_plottitle,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_testing-timeseries.png",
        ),
        ensemble_single_scenario_testing_timeseries_plot,
    )

    Makie.empty!(ensemble_single_scenario_testing_timeseries_plot)

    ensemble_single_scenario_outbreak_dist_plot = ensemble_outbreak_distribution_plot(
        testarr,
        incarr;
        plottitle = noise_plottitle,
    )

    save(
        joinpath(
            ensemble_noise_plotpath,
            "ensemble-sim_single-scenario_outbreak-distribution.png",
        ),
        ensemble_single_scenario_outbreak_dist_plot,
    )

    Makie.empty!(ensemble_single_scenario_outbreak_dist_plot)

    return nothing
end
