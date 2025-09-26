#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack

using CSDNoise
using StructArrays
using SumTypes
using Try
using DataFrames
using StatsBase: StatsBase
using Random: Random
using Distributions: Distributions
using StyledStrings
using Dates
using NLopt
using MultistartOptimization
using FLoops: FLoops

using ProgressMeter

#%%
burnin_years = 5
nyears = 20
burnin_time = 365.0 * burnin_years
ensemble_time_specification = SimTimeParameters(;
    burnin = 365.0 * burnin_years, tmin = 0.0, tmax = 365.0 * nyears,
    tstep = 1.0,
)

ensemble_state_specification = StateParameters(
    500_000,
    Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
)

mu = calculate_mu(27)
beta_mean = calculate_beta(
    R0, GAMMA, mu, 1, ensemble_state_specification.init_states.N
)
epsilon = calculate_import_rate(
    mu, R0, ensemble_state_specification.init_states.N
)

min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
    0.9,
    burnin_years * 2,
    ensemble_state_specification.init_states.S,
    ensemble_state_specification.init_states.N,
    R0,
    mu,
)

max_burnin_vaccination_coverage = 1.0
min_vaccination_coverage = 0.0
max_vaccination_coverage = 0.8

ensemble_dynamics_specification = DynamicsParameterSpecification(
    beta_mean,
    0.0,
    cos,
    SIGMA,
    GAMMA,
    mu,
    27,
    epsilon,
    R0,
    min_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage,
    0.6,
    0.8,
)

null_dynamics_specification = DynamicsParameterSpecification(
    map(
        pn -> getproperty(ensemble_dynamics_specification, pn),
        filter(
            name ->
            name != :min_vaccination_coverage &&
                name != :max_vaccination_coverage,
            propertynames(ensemble_dynamics_specification),
        ),
    )...,
    nothing,
    nothing,
)

#%%
ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

null_specification = EnsembleSpecification(
    ensemble_state_specification,
    null_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

ensemble_outbreak_specification = OutbreakSpecification(
    5, 30, 500
)

#%%
ensemble_dynamics_parameters_sa = StructVector(
    get_ensemble_file(
        ensemble_specification
    )["dynamics_parameters"],
)

null_dynamics_parameters_sa = StructVector(
    get_ensemble_file(
        null_specification
    )["dynamics_parameters"],
)

@assert ensemble_dynamics_parameters_sa.burnin_vaccination_coverage ==
    null_dynamics_parameters_sa.burnin_vaccination_coverage ==
    null_dynamics_parameters_sa.vaccination_coverage

#%%
ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

null_single_scenario_inc_file = get_ensemble_file(
    null_specification, ensemble_outbreak_specification
)

emergent_incidence_arr = ensemble_single_scenario_inc_file["emergent_incidence_arr"]
emergent_outbreak_threshold_vecs = ensemble_single_scenario_inc_file["emergent_outbreak_threshold_vecs"]

null_incidence_arr = null_single_scenario_inc_file["emergent_incidence_arr"]
null_outbreak_threshold_vecs = null_single_scenario_inc_file["emergent_outbreak_threshold_vecs"]

ensemble_single_Reff_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_arr"]

ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

null_single_Reff_arr = get_ensemble_file(
    null_specification
)["ensemble_Reff_arr"]

null_single_Reff_thresholds_vec = get_ensemble_file(
    null_specification
)["ensemble_Reff_thresholds_vec"]

#%%
logfilepath = scriptsdir("ensemble-sim_ews-multistart-optimization.log.txt")

noise_specification_vec = [
    NoiseSpecification(PoissonNoise(1.0)),
    NoiseSpecification(PoissonNoise(7.0)),
    NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.8734)),
    NoiseSpecification(DynamicalNoise(5.0, 7, 14, "in-phase", 0.15, 0.102)),
]

test_specification_vec = [
    IndividualTestSpecification(0.5, 0.5, 0),
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(0.95, 0.95, 0),
    IndividualTestSpecification(0.96, 0.96, 0),
    IndividualTestSpecification(0.97, 0.97, 0),
    IndividualTestSpecification(0.98, 0.98, 0),
    IndividualTestSpecification(0.99, 0.99, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

percent_tested_vec = [1.0]

#%%
ews_method_vec = [
    # Centered,
    Backward,
]
ews_aggregation_vec = [
    # Day(7),
    # Day(14),
    Day(28),
]
ews_bandwidth_vec = [Week(52)]
ews_lag_days_vec = [1]

ews_metric_specification_vec = create_combinations_vec(
    EWSMetricSpecification,
    (
        ews_method_vec,
        ews_aggregation_vec,
        ews_bandwidth_vec,
        ews_lag_days_vec,
    ),
)

ews_metric_vec = [
    "autocorrelation",
    "autocovariance",
    "coefficient_of_variation",
    "index_of_dispersion",
    "kurtosis",
    "mean",
    "skewness",
    "variance",
]

#%%
ews_enddate_type_vec = [
    Reff_start,
    # Reff_end
]
ews_threshold_window_vec = [ExpandingThresholdWindow]
ews_threshold_burnin_vec = [
    # Day(50),
    Year(5),
]

#%%
# Multistart optimization configuration
# Define parameter bounds for continuous optimization
multistart_config = (
    # Parameter bounds (will be optimized continuously)
    quantile_bounds = (0.5, 0.99),           # Threshold quantile range
    consecutive_bounds = (1.0, 30.0),          # Consecutive thresholds range

    # Optimization settings
    n_sobol_points = 100,                      # Number of Sobol sequence starting points
    maxeval = 1000,                            # Maximum evaluations per local optimization
    local_algorithm = NLopt.LN_BOBYQA,         # Local optimization algorithm
    xtol_rel = 1.0e-3,                           # Relative tolerance for parameters
    xtol_abs = 1.0e-3,                           # Absolute tolerance for parameters
    ftol_rel = 1.0e-4,                           # Relative tolerance for function values
    # executor = FLoops.SequentialEx(),          # Executor for parallel processing
    executor = FLoops.ThreadedEx(),          # Executor for parallel processing
)

#%%
specification_vecs = (;
    noise_specification_vec,
    test_specification_vec,
    percent_tested_vec,
    ews_metric_specification_vec,
    ews_enddate_type_vec,
    ews_threshold_window_vec,
    ews_threshold_burnin_vec,
    ews_metric_vec,
)

#%%
println(styled"{green:Starting Multistart EWS Hyperparameter Optimization}")
println(styled"Configuration:")
println(styled"  Sobol points per scenario: {blue:$(multistart_config.n_sobol_points)}")
println(styled"  Max evaluations per local opt: {yellow:$(multistart_config.maxeval)}")
println(styled"  Parameter bounds:")
println(styled"    Quantile: {cyan:$(multistart_config.quantile_bounds)}")
println(styled"    Consecutive: {cyan:$(multistart_config.consecutive_bounds)}")
println(styled"  Burnin (fixed): {cyan:$(ews_threshold_burnin_vec)}")
println(styled"  Executor: {magenta:$(multistart_config.executor)}")

#%%
multistart_optimal_ews_df = ews_multistart_optimization(
    specification_vecs,
    (;
        ensemble_specification,
        emergent_incidence_arr,
        null_incidence_arr,
        ensemble_single_Reff_thresholds_vec,
        emergent_outbreak_threshold_vecs,
    );
    # File management
    filedir = outdir("ensemble", "ews-multistart-optimization"),
    optimization_filename_base = "ews-multistart-optimization.jld2",

    # Multistart configuration
    n_sobol_points = multistart_config.n_sobol_points,
    maxeval = multistart_config.maxeval,
    local_algorithm = multistart_config.local_algorithm,
    xtol_rel = multistart_config.xtol_rel,
    xtol_abs = multistart_config.xtol_abs,
    ftol_rel = multistart_config.ftol_rel,
    executor = multistart_config.executor,

    # Parameter bounds
    quantile_bounds = multistart_config.quantile_bounds,
    consecutive_bounds = multistart_config.consecutive_bounds,

    # Control options
    force = true,
    return_df = true,
    save_results = true,
    verbose = true,
)

#%%
println(styled"\n{green:Multistart Optimization Results Summary}")
println(styled"Total scenarios optimized: {yellow:$(nrow(optimal_ews_df))}")
println(styled"Best overall accuracy: {blue:$(round(maximum(optimal_ews_df.accuracy), digits=4))}")
println(styled"Mean accuracy: {cyan:$(round(mean(optimal_ews_df.accuracy), digits=4))}")

# Show parameter distribution summary
println(styled"\n{green:Optimal Parameter Distributions}")
println(styled"Quantile - Min: {cyan:$(round(minimum(optimal_ews_df.ews_threshold_quantile), digits=3))}, Max: {cyan:$(round(maximum(optimal_ews_df.ews_threshold_quantile), digits=3))}, Mean: {cyan:$(round(mean(optimal_ews_df.ews_threshold_quantile), digits=3))}")
println(styled"Consecutive - Min: {cyan:$(minimum(optimal_ews_df.ews_consecutive_thresholds))}, Max: {cyan:$(maximum(optimal_ews_df.ews_consecutive_thresholds))}, Mean: {cyan:$(round(mean(optimal_ews_df.ews_consecutive_thresholds), digits=1))}")

burnin_days = Dates.days.(optimal_ews_df.ews_threshold_burnin)
println(styled"Burnin (fixed) - Values: {cyan:$(unique(burnin_days))} days")

#%%
create_optimal_ews_plots(
    optimal_multistart_results,
    ensemble_specification,
    emergent_incidence_arr,
    null_incidence_arr,
    ensemble_single_Reff_thresholds_vec,
    ["specificity", "speed"];
    optimal_grouping_parameters = [
        :noise_specification,
        :test_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_metric,
    ],
    output_format = "png",
    force_heatmap = true,
    force_survival = false,
)

# #%%
# # Compare with grid search results if available
# println(styled"\n{green:Comparison with Grid Search (if available)}")
#
# try
#     # Try to load most recent grid search results
#     grid_search_filepath = get_most_recent_hyperparam_filepath(
#         "ews-hyperparam-optimization.jld2",
#         outdir("ensemble", "ews-hyperparam-optimization")
#     )
#
#     if Try.isok(grid_search_filepath)
#         grid_results = load(Try.unwrap(grid_search_filepath))["optimal_ews_df"]
#
#         println(styled"Grid search results found!")
#         println(styled"Grid search scenarios: {yellow:$(nrow(grid_results))}")
#         println(styled"Grid search best accuracy: {blue:$(round(maximum(grid_results.accuracy), digits=4))}")
#
#         # Compare best accuracies by noise type
#         println(styled"\n{green:Accuracy Comparison by Noise Type}")
#         for noise_spec in unique(optimal_ews_df.noise_specification)
#             ms_subset = filter(r -> r.noise_specification == noise_spec, optimal_ews_df)
#             grid_subset = filter(r -> r.noise_specification == noise_spec, grid_results)
#
#             if nrow(grid_subset) > 0
#                 ms_best = maximum(ms_subset.accuracy)
#                 grid_best = maximum(grid_subset.accuracy)
#                 improvement = ms_best - grid_best
#
#                 println(styled"$(get_noise_description(noise_spec)):")
#                 println(styled"  Multistart: {blue:$(round(ms_best, digits=4))}")
#                 println(styled"  Grid search: {cyan:$(round(grid_best, digits=4))}")
#                 println(styled"  Improvement: {$(improvement >= 0 ? "green" : "red"):$(round(improvement, digits=4))}")
#             end
#         end
#     else
#         println(styled"{yellow:No grid search results found for comparison}")
#     end
# catch e
#     println(styled"{yellow:Could not load grid search results: $e}")
# end
#
#%%
# Analysis of optimization efficiency
println(styled"\n{green:Optimization Efficiency Analysis}")

scenarios_count = nrow(optimal_ews_df)
n_sobol_per_scenario = multistart_config.n_sobol_points
total_multistart_evaluations = scenarios_count * n_sobol_per_scenario

println(styled"Scenarios optimized: {blue:$scenarios_count}")
println(styled"Sobol points per scenario: {cyan:$n_sobol_per_scenario}")
println(styled"Total multistart evaluations: {yellow:$total_multistart_evaluations}")

# Estimate equivalent grid search size
# Grid search would need: n_quantiles × n_consecutive evaluations per scenario
# With the original bounds: (0.5:0.01:0.99) × (2:1:30) = 50 × 29 = 1450 per scenario
equivalent_grid_evaluations = scenarios_count * 1450  # Conservative estimate
efficiency_gain = equivalent_grid_evaluations / total_multistart_evaluations

println(styled"\n{green:Efficiency Comparison}")
println(styled"Equivalent grid search evaluations: {red:$equivalent_grid_evaluations}")
println(styled"Multistart evaluations: {green:$total_multistart_evaluations}")
println(styled"Efficiency gain: {blue:$(round(efficiency_gain, digits=1))x fewer evaluations}")

#%%
debug_multistart_plots = true

if debug_multistart_plots
    # Select best performing configuration for detailed analysis
    best_idx = argmax(optimal_ews_df.accuracy)
    best_config = optimal_ews_df[best_idx, :]

    println(styled"\n{green:Best Configuration Analysis}")
    println(styled"Best accuracy: {blue:$(round(best_config.accuracy, digits=4))}")
    println(styled"Noise: {cyan:$(get_noise_description(best_config.noise_specification))}")
    println(styled"Test: {yellow:$(get_test_description(best_config.test_specification))}")
    println(styled"EWS metric: {magenta:$(best_config.ews_metric)}")
    println(styled"Optimal parameters:")
    println(styled"  Quantile: {cyan:$(round(best_config.ews_threshold_quantile, digits=3))}")
    println(styled"  Consecutive: {cyan:$(best_config.ews_consecutive_thresholds)}")
    println(styled"  Burnin (fixed): {cyan:$(best_config.ews_threshold_burnin)}")

    # Parameter distribution analysis
    println(styled"\n{green:Parameter Distribution Analysis}")

    # Group by noise type and show parameter preferences
    for noise_spec in unique(optimal_ews_df.noise_specification)
        subset = filter(r -> r.noise_specification == noise_spec, optimal_ews_df)

        println(styled"\n$(get_noise_description(noise_spec)):")
        println(styled"  Mean quantile: {cyan:$(round(mean(subset.ews_threshold_quantile), digits=3))}")
        println(styled"  Mean consecutive: {cyan:$(round(mean(subset.ews_consecutive_thresholds), digits=1))}")
        println(styled"  Burnin (fixed): {cyan:$(unique(subset.ews_threshold_burnin))}")
        println(styled"  Best accuracy: {blue:$(round(maximum(subset.accuracy), digits=4))}")
    end
end
#
# #%%
# println(styled"\n{green:✅ Multistart EWS Hyperparameter Optimization Complete!}")
# println(styled"Results saved to: {cyan:$(outdir("ensemble", "ews-multistart-optimization"))}")
# println(styled"Plots saved to: {cyan:$(plotsdir())}")
