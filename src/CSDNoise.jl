"""
Note that each file is a module that defines the exported functions.
They are just listed here for convenience of sourcing one file.
"""
module CSDNoise

# using Reexport
include("DrWatson-helpers.jl")
export outdir

include("transmission-functions.jl")
export calculate_beta, calculate_beta_amp,
    calculateReffective, calculateReffective_t!, calculateR0,
    calculate_import_rate, calculate_mu
# @reexport using .TransmissionFunctions

include("dynamics-constants.jl")
export POPULATION_N, LATENT_PER_DAYS, DUR_INF_DAYS, R0, SIGMA, GAMMA,
    LIFE_EXPECTANCY_YEARS, ANNUAL_BIRTHS_PER_K,
    MU, BETA_MEAN, BETA_FORCE, EPSILON

# @reexport using .EWSMetrics

include("structs.jl")
export SimTimeParameters, EnsembleSpecification,
    DynamicsParameterSpecification, DynamicsParameters,
    calculate_vaccination_rate_to_achieve_Reff,
    StateParameters, OutbreakDetectionSpecification,
    OutbreakSpecification, IndividualTestSpecification, get_test_description,
    PoissonNoiseSpecification, DynamicalNoiseSpecification, NoiseSpecification,
    calculate_min_max_vaccination_range,
    get_noise_description, get_noise_magnitude, get_noise_magnitude_description,
    getdirpath,
    EWSMethod, Backward, Centered, method_string,
    EWSMetricSpecification,
    get_ews_metric_specification_description,
    AbstractEWSThresholdWindow,
    ExpandingThresholdWindow,
    RollingThresholdWindow,
    ScenarioSpecification,
    EWSMetrics
# @reexport using .ODStructs

include("ews-metrics.jl")
export aggregate_Reff_vec, aggregate_thresholds_vec, aggregate_timeseries,
    spaero_mean, spaero_mean!,
    spaero_var,
    spaero_cov,
    spaero_iod,
    spaero_skew,
    spaero_kurtosis,
    spaero_autocov,
    spaero_autocor,
    spaero_corkendall,
    compare_against_spaero, filter_spaero_comparison,
    ews_as_df

include("ews-functions.jl")
export calculate_bandwidth, calculate_bandwidth_and_return_ews_metric_spec,
    expanding_ews_thresholds,
    EWSEndDateType, Reff_start, Reff_end, Outbreak_start, Outbreak_end,
    Outbreak_middle,
    calculate_ews_enddate,
    tycho_testing_plots,
    simulation_tau_heatmap_df!,
    tycho_tau_heatmap_df,
    calculate_ews_lead_time, calculate_ews_trigger_index,
    ews_lead_time_df!

include("ews-hyperparam-optimization.jl")
export ews_hyperparam_optimization,
    ews_hyperparam_gridsearch,
    ews_hyperparam_gridsearch!,
    check_missing_ews_hyperparameter_simulations,
    load_most_recent_hyperparam_file,
    get_most_recent_hyperparam_filepath,
    optimal_ews_heatmap_df,
    optimal_ews_heatmap_plot,
    simulate_and_plot_ews_survival,
    ews_survival_plot,
    simulate_ews_survival_data,
    create_ews_survival_data

include(
    "tycho-cleaning.jl"
)
export cdc_week_to_date,
    calculate_aggregation_cases, fill_aggregation_values,
    calculate_ews_enddate

include("test-constants.jl")
export CLINICAL_CASE_TEST_SPEC, EPI_LINKED_CASE_TEST_SPEC, CLINICAL_TEST_SPECS

include("SEIR-model.jl")
export seir_mod, seir_mod!, seir_mod_loop!,
    convert_svec_to_matrix, convert_svec_to_matrix!, convert_svec_to_array
# @reexport using .SEIRModel

include("detection-thresholds.jl")
export create_inc_infec_arr,
    create_inc_infec_arr!,
    calculate_outbreak_thresholds,
    Reff_ge_than_one
# @reexport using .DetectionThresholds

include("diag-testing-functions.jl")
export create_testing_arrs, create_testing_arrs!, calculate_tested!,
    calculate_positives!,
    calculate_true_positives!, calculate_noise_positives!,
    infer_true_positives, calculate_test_positivity_rate,
    calculate_movingavg, calculate_movingavg!,
    calculate_test_positivity
# @reexport using .DiagTestingFunctions

include(
    "ensemble-functions.jl"
)
export create_combinations_vec, create_ensemble_spec_combinations,
    run_ensemble_jump_prob, run_jump_prob,
    get_ensemble_file
# @reexport using .EnsembleFunctions

include("noise-functions.jl")
export create_noise_arr, add_poisson_noise_arr!
# @reexport using .NoiseFunctions

include("plotting-functions/helpers_plots.jl")
export BASE_COLOR, OUTBREAK_COLOR, REFF_GT_ONE_COLOR,
    line_and_hline!

include("plotting-functions/single-simulation_plots.jl")
export incidence_prevalence_plot,
    Reff_plot,
    incidence_testing_plot,
    ensemble_incarr_Reff_plot

include(
    "plotting-functions/tycho_plots.jl"
)
export tycho_epicurve, tycho_noise_components_epicurve,
    tycho_test_positive_components_epicurve,
    tycho_tau_distribution,
    tycho_tau_heatmap_plot,
    ews_lead_time_plot

include(
    "plotting-functions/simulation-ews_plots.jl"
)
export Reff_ews_plot, simulation_tau_distribution

include("ensemble-sim_single-scenario_plots.jl")
export plot_all_single_scenarios

include("plotting-functions/simulation-optimal-ews_plots.jl")
export create_optimal_ews_plots

include("plotting-functions/hyperparam-debugging_plots.jl")
export hyperparam_debugging_Reff_plot

@static if false
    include("../scripts/ensemble-sim.jl")
    include("../scripts/ensemble-sim_inferred-scenario-visualizations.jl")
    include("../scripts/tycho-visualization.jl")
    include("../scripts/ensemble-sim_ews-optimization.jl")
    include("../manuscript/scripts/optimal-thresholds.jl")
end

end
