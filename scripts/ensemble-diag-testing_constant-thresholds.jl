#%%
using DrWatson
@quickactivate "OutbreakDetection"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using CSV: CSV
using DataFrames
using DataFramesMeta
using Statistics

using OutbreakDetection

includet(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
constant_threshold_test_spec_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.7, 0.7, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    CLINICAL_TEST_SPECS...,
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 3),
    IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

constant_threshold_alertthreshold_vec = [4]

R_0_vec = collect(8.0:4.0:20.0)

ensemble_dynamics_spec_vec = create_combinations_vec(
    DynamicsParameters,
    (
        [ensemble_state_specification.init_states.N],
        [27],
        [0.2],
        [SIGMA],
        [GAMMA],
        R_0_vec,
        [0.8],
    ),
)

ensemble_spec_vec = create_combinations_vec(
    EnsembleSpecification,
    (
        [ensemble_model_type],
        [ensemble_state_specification],
        ensemble_dynamics_spec_vec,
        [ensemble_time_specification],
        [ensemble_nsims],
    ),
)

alert_method_vec = ["movingavg", "dailythreshold_movingavg"]

#%%
for (
    ensemble_noise_specification,
    ensemble_specification,
    alertmethod,
    alertthreshold,
) in
    Iterators.product(
    ensemble_noise_specification_vec, ensemble_spec_vec, alert_method_vec,
    constant_threshold_alertthreshold_vec,
)
    @info "Creating tables for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification)), $(alertmethod), constant threshold: $(alertthreshold)"
    println("==============================================")

    optimal_threshold_comparison_params = (
        alertthreshold_vec = [alertthreshold],
        ensemble_specification = ensemble_specification,
        noise_specification = ensemble_noise_specification,
        outbreak_specification = ensemble_outbreak_specification,
        moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
        percent_visit_clinic = ensemble_percent_visit_clinic,
        alertmethod = alertmethod,
    )

    constant_thresholds_vec = calculate_OptimalThresholdCharacteristics(
        ensemble_percent_clinic_tested_vec,
        constant_threshold_test_spec_vec,
        optimal_threshold_comparison_params,
    )

    noise_specification_path = getdirpath(ensemble_noise_specification)
    noisespec_alertmethod_path = joinpath(noise_specification_path, alertmethod)

    if alertmethod != "dailythreshold"
        noisespec_alertmethod_path = joinpath(
            noisespec_alertmethod_path,
            "moveavglag_$(ensemble_moving_avg_detection_lag)",
        )
    end

    noisespec_alertmethod_filename = replace(
        noisespec_alertmethod_path,
        "/" => "_"
    )

    basedirpath = joinpath(
        "R0_$(ensemble_specification.dynamics_parameters.R_0)",
        noisespec_alertmethod_path,
    )

    baseplotdirpath = joinpath(
        plotsdir("ensemble/constant-thresholds"),
        "threshold_$(alertthreshold)",
        basedirpath,
    )

    clinictested_plotdirpath = joinpath(baseplotdirpath, "clinic-tested")

    compare_optimal_thresholds_chars_plot(
        constant_thresholds_vec,
        [
            (
                char = :accuracy,
                label = "Accuracy",
                color = (ACCURACY_COLOR, 0.7),
                binwidth = 0.01,
            ),
            (
                char = :detectiondelays,
                label = "Detection Delay (Days)",
                color = (DETECTION_DELAY_COLOR, 1.0),
                binwidth = 10,
            ),
            (
                char = :unavoidable_cases,
                label = "Unavoidable Cases",
                color = (PERC_OUTBREAKS_MISSED_COLOR, 1.0),
                binwidth = 500,
            ),
            (
                char = :avoidable_cases,
                label = "Avoidable Cases",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
                binwidth = 500,
            ),
        ];
        plotdirpath = clinictested_plotdirpath,
    )

    test_plotdirpath = joinpath(baseplotdirpath, "tests")

    compare_optimal_thresholds_test_chars_plot(
        constant_thresholds_vec,
        [
            (
                char = :accuracy,
                label = "Accuracy",
                color = (ACCURACY_COLOR, 0.7),
                binwidth = 0.01,
            ),
            (
                char = :detectiondelays,
                label = "Detection Delay (Days)",
                color = (DETECTION_DELAY_COLOR, 1.0),
                binwidth = 10,
            ),
            (
                char = :unavoidable_cases,
                label = "Unavoidable Cases",
                color = (PERC_OUTBREAKS_MISSED_COLOR, 1.0),
                binwidth = 500,
            ),
            (
                char = :avoidable_cases,
                label = "Avoidable Cases",
                color = (PERC_OUTBREAKS_DETECTED_COLOR, 1.0),
                binwidth = 500,
            ),
        ];
        plotdirpath = test_plotdirpath,
    )

    cfr_df = CSV.read(
        datadir("CFR_2022.csv"),
        DataFrame; delim = ',',
        header = true,
        types = [String, Int64, Float64],
    )
    dropmissing!(cfr_df)

    gha_cfr = only(cfr_df[cfr_df.country .== "GHA", :CFR])

    population_df = CSV.read(
        datadir("input-populations.csv"),
        DataFrame; delim = ',',
        header = true
    )

    gha_2022_pop = only(
        population_df[population_df.ISO3_code .== "GHA", "2022"]
    )
    gha_2022_scale_population =
        gha_2022_pop / ensemble_state_specification.init_states.N

    countries = [
        (;
            name = "Ghana",
            code = "GHA",
            cfr = gha_cfr,
            year = "2022",
            population_size = gha_2022_pop,
            scale_population = gha_2022_scale_population,
        ),
    ]

    tabledirpath = joinpath(
        outdir("ensemble/constant-thresholds"),
        "threshold_$(alertthreshold)",
        basedirpath,
    )

    tablefilename = "constant-threshold_threshold-$(alertthreshold)_$(noisespec_alertmethod_filename)_thresholds"

    create_and_save_xlsx_optimal_threshold_summaries(
        constant_thresholds_vec;
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        constant_thresholds_vec, :detectiondelays;
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        constant_thresholds_vec,
        :unavoidable_cases;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        constant_thresholds_vec, :avoidable_cases;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        constant_thresholds_vec, :n_outbreak_cases;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        constant_thresholds_vec, :n_tests;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    @info "All plots and tables saved for $(ensemble_noise_specification.noise_type)"
    println("==============================================")
end
