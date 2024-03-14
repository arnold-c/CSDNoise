#%%
using DrWatson
@quickactivate "CSDNoise"

using ProgressMeter
using FLoops
using NaNMath: NaNMath
using CSV: CSV
using DataFrames
using DataFramesMeta
using Statistics

using CSDNoise

includet(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
optimal_threshold_test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 3),
    IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

optimal_threshold_alertthreshold_vec = collect(10:4:50)

R_0_vec = collect(16.0)

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

alert_method_vec = ["inferred_movingavg"]

#%%
for (ensemble_noise_specification, ensemble_specification, alertmethod) in
    Iterators.product(
    ensemble_noise_specification_vec, ensemble_spec_vec, alert_method_vec
)
    @info "Creating plots and tables for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification)), $(alertmethod)"
    println("==============================================")

    optimal_threshold_comparison_params = (
        alertthreshold_vec = optimal_threshold_alertthreshold_vec,
        ensemble_specification = ensemble_specification,
        noise_specification = ensemble_noise_specification,
        outbreak_specification = ensemble_outbreak_specification,
        moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
        percent_visit_clinic = ensemble_percent_visit_clinic,
        alertmethod = alertmethod,
    )

    optimal_thresholds_vec = calculate_OptimalThresholdCharacteristics(
        ensemble_percent_clinic_tested_vec,
        optimal_threshold_test_spec_vec,
        optimal_threshold_comparison_params,
    )

    noise_specification_path = getdirpath(ensemble_noise_specification)
    noisespec_alertmethod_path = joinpath(noise_specification_path, alertmethod)

    if alertmethod == "inferred_movingavg"
        noisespec_alertmethod_path = joinpath(
            noisespec_alertmethod_path,
            "inferred_movingavg_$(ensemble_moving_avg_detection_lag)",
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
        plotsdir("ensemble/optimal-thresholds"),
        basedirpath
    )

    clinictested_plotdirpath = joinpath(baseplotdirpath, "clinic-tested")

    compare_optimal_thresholds_chars_plot(
        optimal_thresholds_vec,
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
        optimal_thresholds_vec,
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
        silencewarnings = true,
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
        outdir("ensemble/optimal-threshold-results"),
        basedirpath
    )

    tablefilename = "optimal-threshold_$(noisespec_alertmethod_filename)_thresholds"

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec;
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec, :detectiondelays;
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec,
        :unavoidable_cases;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec, :avoidable_cases;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec, :n_outbreak_cases;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec, :n_tests;
        scale_annual = 1 / nyears,
        countries = countries,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec, :noutbreaks;
        scale_annual = 1 / nyears,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    create_and_save_xlsx_optimal_threshold_summaries(
        optimal_thresholds_vec, :nalerts;
        scale_annual = 1 / nyears,
        tabledirpath = tabledirpath,
        filename = tablefilename,
    )

    GC.gc(true)

    @info "All plots and tables saved for R0: $(ensemble_specification.dynamics_parameters.R_0), $(getdirpath(ensemble_noise_specification))"
    println()
    println("==============================================")
    println()
end
