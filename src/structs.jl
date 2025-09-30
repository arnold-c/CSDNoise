# module ODStructs
#
# export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
#     StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
#     OutbreakSpecification, IndividualTestSpecification, NoiseSpecification

using StaticArrays
using LabelledArrays
using StructArrays
using Distributions: Distributions
using Random: Random
using UnPack: @unpack
using Dates: Dates
using LightSumTypes: @sumtype, variant
using NLopt: NLopt

# include("transmission-functions.jl")
# using .TransmissionFunctions

struct SimTimeParameters
    burnin::Float64
    tmin::Float64
    tmax::Float64
    tstep::Float64
    trange::StepRangeLen{Float64, Float64, Float64, Int64}
    tspan::Tuple{Float64, Float64}
    tlength::Int64
end

function SimTimeParameters(;
        burnin = 0.0, tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0
    )
    @assert burnin <= tmax
    return SimTimeParameters(
        burnin, tmin, tmax, tstep, tmin:tstep:tmax, (tmin, tmax),
        length(tmin:tstep:tmax),
    )
end

# Define seasonality function variants
abstract type AbstractSeasonalityFunction end
struct CosineSeasonality end
struct SineSeasonality end
@sumtype SeasonalityFunction(CosineSeasonality, SineSeasonality) <: AbstractSeasonalityFunction


struct DynamicsParameterSpecification
    contact_matrix::Matrix{Int64}
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    min_burnin_vaccination_coverage::Float64
    max_burnin_vaccination_coverage::Float64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
end

# function DynamicsParameterSpecification(
#         beta_mean,
#         beta_force,
#         seasonality,
#         sigma,
#         gamma,
#         mu,
#         annual_births_per_k,
#         epsilon,
#         R_0,
#         burnin_vaccination_params::Tuple{
#             Nothing, Float64, Float64, Int64,
#         },
#         vaccination_params::Union{
#             Tuple{Float64, Float64}, Tuple{Nothing, Nothing},
#         },
#         initial_states,
#     )
#     return DynamicsParameterSpecification(
#         beta_mean,
#         beta_force,
#         seasonality,
#         sigma,
#         gamma,
#         mu,
#         annual_births_per_k,
#         epsilon,
#         R_0,
#         calculate_vaccination_rate_to_achieve_Reff(
#             burnin_vaccination_params[3],
#             burnin_vaccination_params[4],
#             initial_states.S,
#             initial_states.N,
#             R_0,
#             mu,
#         ),
#         burnin_vaccination_params[2],
#         vaccination_params[1],
#         vaccination_params[2],
#     )
# end

function DynamicsParameterSpecification(
        contact_matrix,
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_vaccination_coverage::Nothing,
        max_vaccination_coverage::Nothing,
    )
    return DynamicsParameterSpecification(
        contact_matrix,
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
    )
end

"""
    calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, initial_states, dynamics_specification
    )

Assuming that in the burnin period the rate of infections is negligible, the time to reach a target Reff can be calculated by the difference in the rates in and out of the Susceptible group.

μN(1-ρ) -> S -> μS
dS/dt = μ(N - Nρ - S)
"""
function calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, initial_states, R_0, mu
    )
    @unpack S, N = initial_states

    return calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, S, N, R_0, mu
    )
end

function calculate_vaccination_rate_to_achieve_Reff(
        target_Reff, target_years, S, N, R_0, mu
    )
    @assert target_Reff > 0
    @assert target_Reff < 1.2

    vaccination_coverage =
        1 -
        (((target_Reff * N) / R_0 - S) / (365 * target_years * mu) + S) / N

    if vaccination_coverage > 1 || vaccination_coverage < 0
        return error(
            "Target Reff cannot be reached in the burn-in period. Initial Reff = $(R_0 * S / N). Try a longer burn-in period or a smaller target Reff."
        )
    end
    return round(vaccination_coverage; digits = 4)
end

struct DynamicsParameters
    contact_matrix::Matrix{Int64}
    beta_mean::Float64
    beta_force::Float64
    seasonality::SeasonalityFunction
    sigma::Float64
    gamma::Float64
    mu::Float64
    annual_births_per_k::Float64
    epsilon::Float64
    R_0::Float64
    min_burnin_vaccination_coverage::Float64
    max_burnin_vaccination_coverage::Float64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
    burnin_vaccination_coverage::Float64
    vaccination_coverage::Float64
end

function DynamicsParameters(
        dynamic_parameter_specification::DynamicsParameterSpecification; seed = 1234
    )
    Random.seed!(seed)

    burnin_vaccination_coverage =
    if dynamic_parameter_specification.min_burnin_vaccination_coverage ==
            dynamic_parameter_specification.max_burnin_vaccination_coverage
        dynamic_parameter_specification.min_burnin_vaccination_coverage
    else
        sample_vaccination_coverage(
            dynamic_parameter_specification.min_burnin_vaccination_coverage,
            dynamic_parameter_specification.max_burnin_vaccination_coverage,
        )
    end

    vaccination_coverage =
    if dynamic_parameter_specification.min_burnin_vaccination_coverage ==
            dynamic_parameter_specification.min_vaccination_coverage &&
            dynamic_parameter_specification.max_burnin_vaccination_coverage ==
            dynamic_parameter_specification.max_vaccination_coverage
        burnin_vaccination_coverage
    else
        sample_vaccination_coverage(
            dynamic_parameter_specification.min_vaccination_coverage,
            dynamic_parameter_specification.max_vaccination_coverage,
        )
    end

    dynamics_parameters = DynamicsParameters(
        dynamic_parameter_specification.contact_matrix,
        dynamic_parameter_specification.beta_mean,
        dynamic_parameter_specification.beta_force,
        dynamic_parameter_specification.seasonality,
        dynamic_parameter_specification.sigma,
        dynamic_parameter_specification.gamma,
        dynamic_parameter_specification.mu,
        dynamic_parameter_specification.annual_births_per_k,
        dynamic_parameter_specification.epsilon,
        dynamic_parameter_specification.R_0,
        dynamic_parameter_specification.min_burnin_vaccination_coverage,
        dynamic_parameter_specification.max_burnin_vaccination_coverage,
        dynamic_parameter_specification.min_vaccination_coverage,
        dynamic_parameter_specification.max_vaccination_coverage,
        burnin_vaccination_coverage,
        vaccination_coverage,
    )
    return dynamics_parameters
end

function sample_vaccination_coverage(
        min_coverage,
        max_coverage,
        digits = 4
    )
    return round(
        rand(
            Distributions.Uniform(
                min_coverage,
                max_coverage
            )
        );
        digits = digits
    )
end


struct StateParameters
    init_states::SLArray{Tuple{5}, Int64, 1, 5, (:S, :E, :I, :R, :N)}
    init_state_props::SLArray{Tuple{4}, Float64, 1, 4, (:s_prop, :e_prop, :i_prop, :r_prop)}
end

function StateParameters(N::Int64, init_state_props::Dict)
    return StateParameters(;
        N = N,
        s_prop = init_state_props[:s_prop],
        e_prop = init_state_props[:e_prop],
        i_prop = init_state_props[:i_prop],
    )
end

function StateParameters(;
        N = 500_00, s_prop = 0.1, e_prop = 0.01, i_prop = 0.01
    )
    r_prop = 1 - (s_prop + e_prop + i_prop)

    states = SLVector(;
        S = Int64(round(s_prop * N)),
        E = Int64(round(e_prop * N)),
        I = Int64(round(i_prop * N)),
        R = Int64(round(r_prop * N)),
        N = N,
    )
    state_props = SLVector(;
        s_prop = s_prop,
        e_prop = e_prop,
        i_prop = i_prop,
        r_prop = r_prop,
    )

    return StateParameters(
        states, state_props
    )
end

struct EnsembleSpecification
    state_parameters::StateParameters
    dynamics_parameter_specification::DynamicsParameterSpecification
    time_parameters::SimTimeParameters
    nsims::Int64
    dirpath::String
end

function EnsembleSpecification(
        state_parameters::StateParameters,
        dynamics_parameter_specification::DynamicsParameterSpecification,
        time_parameters::SimTimeParameters,
        nsims::Int64,
    )
    dirpath = outdir(
        "ensemble",
        "seasonal-infectivity-import",
        "tau-leaping",
        "N_$(state_parameters.init_states.N)",
        "r_$(state_parameters.init_state_props.r_prop)",
        "nsims_$(nsims)",
        "R0_$(dynamics_parameter_specification.R_0)",
        "latent_period_$(round(1 / dynamics_parameter_specification.sigma; digits = 2))",
        "infectious_period_$(round(1 / dynamics_parameter_specification.gamma; digits = 2))",
        "min_burnin_vaccination_coverage_$(dynamics_parameter_specification.min_burnin_vaccination_coverage)",
        "max_burnin_vaccination_coverage_$(dynamics_parameter_specification.max_burnin_vaccination_coverage)",
        "min_vaccination_coverage_$(dynamics_parameter_specification.min_vaccination_coverage)",
        "max_vaccination_coverage_$(dynamics_parameter_specification.max_vaccination_coverage)",
        "births_per_k_$(dynamics_parameter_specification.annual_births_per_k)",
        "beta_force_$(dynamics_parameter_specification.beta_force)",
        "burnin_$(time_parameters.burnin)",
        "tmax_$(time_parameters.tmax)",
        "tstep_$(time_parameters.tstep)",
    )

    return EnsembleSpecification(
        state_parameters,
        dynamics_parameter_specification,
        time_parameters,
        nsims,
        dirpath,
    )
end

struct OutbreakSpecification
    outbreak_threshold::Int64
    minimum_outbreak_duration::Int64
    minimum_outbreak_size::Int64
    dirpath::String
end

function OutbreakSpecification(
        outbreak_threshold, minimum_outbreak_duration, minimum_outbreak_size
    )
    dirpath = joinpath(
        "min_outbreak_dur_$(minimum_outbreak_duration)",
        "min_outbreak_size_$(minimum_outbreak_size)",
        "outbreak_threshold_$(outbreak_threshold)",
    )

    return OutbreakSpecification(
        outbreak_threshold,
        minimum_outbreak_duration,
        minimum_outbreak_size,
        dirpath,
    )
end

struct AlertMethod
    method_name::String
    function AlertMethod(method_name::String)
        available_test_methods = [
            "dailythreshold", "movingavg", "dailythreshold_movingavg",
            "inferred_movingavg",
        ]
        if !in(method_name, available_test_methods)
            error(
                "$(method_name) is not a valid test method. It must be one of $(available_test_methods)"
            )
        end
        return new(method_name)
    end
end

struct OutbreakDetectionSpecification
    alert_threshold::Int64
    moving_average_lag::Int64
    percent_visit_clinic::Float64
    percent_clinic_tested::Float64
    percent_tested::Float64
    alert_method::AlertMethod
    dirpath::String
end

function OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        alert_method,
    )
    alertdirpath = joinpath(
        "alertmethod_$(alert_method)", "alertthreshold_$(alert_threshold)"
    )
    testingdirpath = joinpath(
        "perc_visit_clinic_$(percent_visit_clinic)",
        "perc_clinic_tested_$(percent_clinic_tested)",
    )

    dirpath = if alert_method == "dailythreshold"
        joinpath(
            alertdirpath,
            testingdirpath,
        )
    else
        joinpath(
            alertdirpath,
            "moveavglag_$(moving_average_lag)",
            testingdirpath,
        )
    end

    return OutbreakDetectionSpecification(
        alert_threshold,
        moving_average_lag,
        percent_visit_clinic,
        percent_clinic_tested,
        percent_visit_clinic * percent_clinic_tested,
        AlertMethod(alert_method),
        dirpath,
    )
end

struct IndividualTestSpecification
    sensitivity::Float64
    specificity::Float64
    test_result_lag::Int64
end

function get_test_description(test_specification::IndividualTestSpecification)
    test_specification == IndividualTestSpecification(1.0, 0.0, 0) && return "Clinical case definition"
    test_specification.sensitivity == test_specification.specificity < 1.0 && return "Imperfect Test ($(Int64(test_specification.sensitivity * 100))% Sensitive & Specific)"
    (test_specification.sensitivity < 1.0 || test_specification.specificity < 1.0) && return "Imperfect Test ($(Int64(test_specification.sensitivity * 100))% Sensitive & $(Int64(test_specification.specificity * 100))% Specific)"
    test_specification.sensitivity == test_specification.specificity == 1.0 && return "Perfect Test"
    @error "Don't have a description matching the test specification"
    return
end

abstract type AbstractNoiseSpecification end

struct PoissonNoise
    noise_mean_scaling::Float64
end

struct DynamicalNoise
    R_0::Float64
    latent_period::Int64
    duration_infection::Int64
    correlation::String
    noise_mean_scaling::Float64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
end

@sumtype NoiseSpecification(PoissonNoise, DynamicalNoise) <: AbstractNoiseSpecification


function calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range = 0.2,
    )
    @assert mean_vaccination_coverage <= 1.0

    min_vaccination_range = minimum(
        [
            max_vaccination_range,
            1.0 - mean_vaccination_coverage,
            mean_vaccination_coverage,
        ]
    )

    min_vaccination_coverage = round(
        mean_vaccination_coverage - min_vaccination_range; digits = 4
    )
    max_vaccination_coverage = round(
        mean_vaccination_coverage + min_vaccination_range; digits = 4
    )
    return min_vaccination_coverage, max_vaccination_coverage
end

struct DynamicalNoiseSpecification
    R0::Float64
    latent_period::Int64
    duration_infection::Int64
    correlation::String
    poisson_component::Float64
    vaccination_bounds::Vector{Float64}
    susceptible_bounds::Vector{Float64}
    max_vaccination_range::Float64
    function DynamicalNoiseSpecification(
            R0::Float64,
            latent_period::Int64,
            duration_infection::Int64,
            correlation::String,
            poisson_component::Float64,
            vaccination_bounds::Vector{Float64},
            susceptible_bounds::Vector{Float64},
            max_vaccination_range::Float64,
        )

        @assert length(vaccination_bounds) == 2
        @assert vaccination_bounds[1] < vaccination_bounds[2]
        @assert length(susceptible_bounds) == 2
        @assert susceptible_bounds[1] < susceptible_bounds[2]
        return new(
            R0,
            latent_period,
            duration_infection,
            correlation,
            poisson_component,
            vaccination_bounds,
            susceptible_bounds,
            max_vaccination_range

        )
    end
end

function DynamicalNoiseSpecification(;
        R0::Float64,
        latent_period::Int64,
        duration_infection::Int64,
        correlation::String,
        poisson_component::Float64,
        vaccination_bounds::Vector{Float64} = [0.0, 1.0],
        susceptible_bounds::Vector{Float64} = [0.01, 0.99],
        max_vaccination_range::Float64 = 0.2
    )
    return DynamicalNoiseSpecification(
        R0,
        latent_period,
        duration_infection,
        correlation,
        poisson_component,
        vaccination_bounds,
        susceptible_bounds,
        max_vaccination_range
    )
end

struct NoiseVaccinationOptimizationParameters
    n_sobol_points::Int64
    local_algorithm
    maxeval::Int64
    xtol_rel::Float64
    xtol_abs::Float64
    atol::Float64
end

"""
    NoiseVaccinationOptimizationParameters(;
        n_sobol_points::Int64 = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval::Int64 = 1000,
        xtol_rel::Float64 = 1.0e-3,
        xtol_abs::Float64 = 1.0e-3,
        atol::Float64 = 1.0e-4
    )

Parameters for optimizing noise vaccination levels.

# Arguments
- `n_sobol_points::Int64`: Number of Sobol sequence points for global search
- `local_algorithm`: NLopt algorithm for local optimization
- `maxeval::Int64`: Maximum number of function evaluations for local optimization
- `xtol_rel::Float64`: Relative tolerance on parameter changes for local optimization
- `xtol_abs::Float64`: Absolute tolerance on parameter changes for local optimization
- `atol::Float64`: Absolute difference tolerance threshold; program errors if not met
"""
function NoiseVaccinationOptimizationParameters(;
        n_sobol_points::Int64 = 100,
        local_algorithm = NLopt.LN_BOBYQA,
        maxeval::Int64 = 1000,
        xtol_rel::Float64 = 1.0e-3,
        xtol_abs::Float64 = 1.0e-3,
        atol::Float64 = 1.0e-3
    )
    return NoiseVaccinationOptimizationParameters(
        n_sobol_points,
        local_algorithm,
        maxeval,
        xtol_rel,
        xtol_abs,
        atol
    )
end

get_noise_description(noise_specification::NoiseSpecification) = get_noise_description(variant(noise_specification))

function get_noise_description(noise_specification::PoissonNoise)
    return "poisson"
end

function get_noise_description(noise_specification::DynamicalNoise)
    return string("dynamical, ", noise_specification.correlation)
end

get_noise_magnitude_description(noise_specification::NoiseSpecification) = get_noise_magnitude_description(variant(noise_specification))

function get_noise_magnitude_description(noise_specification::Union{PoissonNoise, DynamicalNoise})
    return string("Poisson scaling: ", noise_specification.noise_mean_scaling)
end

get_noise_magnitude(noise_specification::NoiseSpecification) = get_noise_magnitude(variant(noise_specification))

function get_noise_magnitude(noise_specification::Union{PoissonNoise, DynamicalNoise})
    return noise_specification.noise_mean_scaling
end

noise_table_description(noise_specification::NoiseSpecification) = noise_table_description(variant(noise_specification))

function noise_table_description(noise_specification::PoissonNoise)
    noise_scaling = if noise_specification.noise_mean_scaling == 7
        "High"
    elseif noise_specification.noise_mean_scaling == 1
        "Low"
    else
        "Uncharacterized"
    end
    return "$(noise_scaling) Static Noise"
end

function noise_table_description(noise_specification::DynamicalNoise)
    avg_vaccination = round(
        mean(
            [
                noise_specification.min_vaccination_coverage,
                noise_specification.max_vaccination_coverage,
            ]
        );
        digits = 4,
    )
    noise_scaling = if avg_vaccination == 0.102
        "High"
    elseif avg_vaccination == 0.8734
        "Low"
    else
        "Uncharacterized"
    end
    return "$(noise_scaling) Dynamical Noise"
end

# function get_noise_magnitude(
#     noise_specification::DynamicalNoiseSpecification
# )
#     return string("Rubella vax: ", noise_specification.vaccination_coverage)
# end

function getdirpath(spec::NoiseSpecification)
    return getdirpath(variant(spec))
end

function getdirpath(spec::Union{PoissonNoise, DynamicalNoise})
    return reduce(
        joinpath,
        map(
            p -> "$(p)_$(getproperty(spec, p))",
            propertynames(spec),
        ),
    )
end

struct SEIRRun
    states::Vector{SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

struct NoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
    mean_poisson_noise::Float64
    mean_dynamic_noise::Float64
end

abstract type AbstractThresholds end

struct Thresholds <: AbstractThresholds
    lower_bounds::Vector{Int64}
    upper_bounds::Vector{Int64}
    duration::Vector{Int64}
end

struct OutbreakThresholds <: AbstractThresholds
    lower_bounds::Vector{Int64}
    upper_bounds::Vector{Int64}
    duration::Vector{Int64}
    num_infections_during_bounds::Vector{Int64}
end


abstract type AbstractEWSMethod end
struct Backward end
struct Centered end
@sumtype EWSMethod(Backward, Centered) <: AbstractEWSMethod

struct EWSMetricSpecification
    method::EWSMethod
    aggregation::Dates.Day
    bandwidth::Dates.Day
    lag::Int64
    dirpath::String
end

function EWSMetricSpecification(
        method::EWSMethod,
        aggregation::Dates.Day,
        bandwidth::Dates.Day,
        lag::Int64,
    )
    aggregation_days_val = Dates.value(aggregation)
    bandwidth_days_val = Dates.value(bandwidth)

    return EWSMetricSpecification(
        method,
        aggregation,
        bandwidth,
        lag,
        _EWSMetricSpecification_path(
            method,
            aggregation_days_val,
            bandwidth_days_val,
            lag,
        ),
    )
end

function EWSMetricSpecification(
        method::EWSMethod, aggregation::Dates.DatePeriod, bandwidth::Dates.DatePeriod, lag::Int64
    )
    aggregation_days_val = Dates.days(aggregation)
    bandwidth_days_val = Dates.days(bandwidth)
    aggregation_days = Dates.Day(aggregation_days_val)
    bandwidth_days = Dates.Day(bandwidth_days_val)

    return EWSMetricSpecification(
        method,
        aggregation_days,
        bandwidth_days,
        lag,
        _EWSMetricSpecification_path(
            method,
            aggregation_days_val,
            bandwidth_days_val,
            lag,
        ),
    )
end

function _EWSMetricSpecification_path(
        method::EWSMethod,
        aggregation::Int64,
        bandwidth::Int64,
        lag::Int64,
    )
    return joinpath(
        "ews-method_$(method_string(method))",
        "ews-aggregation-days_$(aggregation)",
        "ews-bandwidth-days_$(bandwidth)",
        "ews-lag_$(lag)",
    )
end

function get_ews_metric_specification_description(ews_metric_specification)
    return "Method: $(method_string(ews_metric_specification.method)), Aggregation: $(ews_metric_specification.aggregation), Bandwidth: $(ews_metric_specification.bandwidth), Lag: $(ews_metric_specification.lag)"
end

method_string(method::EWSMethod) = lowercase(split(string(method), "::")[1])

struct EWSMetrics
    ews_specification::EWSMetricSpecification
    mean::Vector{Float64}
    variance::Vector{Float64}
    coefficient_of_variation::Vector{Float64}
    index_of_dispersion::Vector{Float64}
    skewness::Vector{Float64}
    kurtosis::Vector{Float64}
    autocovariance::Vector{Float64}
    autocorrelation::Vector{Float64}
    mean_tau::Float64
    variance_tau::Float64
    coefficient_of_variation_tau::Float64
    index_of_dispersion_tau::Float64
    skewness_tau::Float64
    kurtosis_tau::Float64
    autocovariance_tau::Float64
    autocorrelation_tau::Float64
end

abstract type AbstractEWSThresholdWindowType end
struct ExpandingThresholdWindow end
struct RollingThresholdWindow end
@sumtype EWSThresholdWindowType(ExpandingThresholdWindow, RollingThresholdWindow) <: AbstractEWSThresholdWindowType

struct ScenarioSpecification
    ensemble_specification::EnsembleSpecification
    outbreak_specification::OutbreakSpecification
    noise_specification::NoiseSpecification
    outbreak_detection_specification::OutbreakDetectionSpecification
    individual_test_specification::IndividualTestSpecification
    ewsmetric_specification::EWSMetricSpecification
    dirpath::String
end

function ScenarioSpecification(
        ensemble_specification::EnsembleSpecification,
        outbreak_specification::OutbreakSpecification,
        noise_specification::NoiseSpecification,
        outbreak_detection_specification::OutbreakDetectionSpecification,
        individual_test_specification::IndividualTestSpecification,
        ewsmetric_specification::EWSMetricSpecification,
    )
    dirpath = joinpath(
        ensemble_specification.dirpath,
        outbreak_specification.dirpath,
        getdirpath(noise_specification),
        outbreak_detection_specification.dirpath,
        "testsens_$(individual_test_specification.sensitivity)",
        "testspec_$(individual_test_specification.specificity)",
        "testlag_$(individual_test_specification.test_result_lag)",
        ewsmetric_specification.dirpath,
    )

    return ScenarioSpecification(
        ensemble_specification,
        outbreak_specification,
        noise_specification,
        outbreak_detection_specification,
        individual_test_specification,
        ewsmetric_specification,
        dirpath,
    )
end

abstract type AbstractEWSEndDateType end

struct Reff_start end
struct Reff_end end
struct Outbreak_start end
struct Outbreak_middle end
struct Outbreak_end end

@sumtype EWSEndDateType(Reff_start, Reff_end, Outbreak_start, Outbreak_middle, Outbreak_end) <: AbstractEWSEndDateType

"""
    CachedSimulationData

Pre-computed simulation data that can be reused across parameter evaluations.
This avoids expensive recomputation of noise arrays and test arrays.
"""
struct CachedSimulationData
    testarr::Array{Int64, 3}
    null_testarr::Array{Int64, 3}
    thresholds::Vector{Matrix{Int64}}
    ews_metrics::Vector{EWSMetrics}
    null_ews_metrics::Vector{EWSMetrics}
end

"""
    OptimizationScenario

Struct representing a single optimization scenario with all necessary parameters
for EWS hyperparameter optimization.
"""
struct OptimizationScenario
    ensemble_specification::EnsembleSpecification
    null_specification::EnsembleSpecification
    noise_specification::NoiseSpecification
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_threshold_burnin::Dates.Day
    ews_metric::String
end

function OptimizationScenario(
        ensemble_specification::EnsembleSpecification,
        null_specification::EnsembleSpecification,
        noise_specification::NoiseSpecification,
        test_specification::IndividualTestSpecification,
        percent_tested::Float64,
        ews_metric_specification::EWSMetricSpecification,
        ews_enddate_type::EWSEndDateType,
        ews_threshold_window::EWSThresholdWindowType,
        threshold_burnin::P,
        ews_metric::String,
    ) where {P <: Dates.Period}

    ews_threshold_burnin = Dates.Day(round(Dates.days(threshold_burnin)))

    return OptimizationScenario(
        ensemble_specification,
        null_specification,
        noise_specification,
        test_specification,
        percent_tested,
        ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window,
        ews_threshold_burnin,
        ews_metric,
    )
end

"""
    GridSearchScenario

Scenario for grid search including both base scenario and grid parameters.
"""
struct GridSearchScenario
    # Base scenario fields (from OptimizationScenario)
    ensemble_specification::EnsembleSpecification
    null_specification::EnsembleSpecification
    noise_specification::NoiseSpecification
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_threshold_burnin::Dates.Day
    ews_metric::String
    # Grid search parameters
    threshold_quantile::Float64
    consecutive_thresholds::Int64
end

function GridSearchScenario(
        ensemble_specification::EnsembleSpecification,
        null_specification::EnsembleSpecification,
        noise_specification::NoiseSpecification,
        test_specification::IndividualTestSpecification,
        percent_tested::Float64,
        ews_metric_specification::EWSMetricSpecification,
        ews_enddate_type::EWSEndDateType,
        ews_threshold_window::EWSThresholdWindowType,
        threshold_burnin::P,
        ews_metric::String,
        threshold_quantile::Float64,
        consecutive_thresholds::Int64
    ) where {P <: Dates.Period}

    ews_threshold_burnin = Dates.Day(round(Dates.days(threshold_burnin)))

    return GridSearchScenario(
        ensemble_specification,
        null_specification,
        noise_specification,
        test_specification,
        percent_tested,
        ews_metric_specification,
        ews_enddate_type,
        ews_threshold_window,
        ews_threshold_burnin,
        ews_metric,
        threshold_quantile,
        consecutive_thresholds
    )
end

struct OptimizationResult
    # Scenario fields (from OptimizationScenario)
    ensemble_specification::EnsembleSpecification
    null_specification::EnsembleSpecification
    noise_specification::NoiseSpecification
    test_specification::IndividualTestSpecification
    percent_tested::Float64
    ews_metric_specification::EWSMetricSpecification
    ews_enddate_type::EWSEndDateType
    ews_threshold_window::EWSThresholdWindowType
    ews_threshold_burnin::Dates.Day
    ews_metric::String
    # Result fields (from OptimizedValues)
    threshold_quantile::Float64
    consecutive_thresholds::Int64
    accuracy::Float64
    sensitivity::Float64
    specificity::Float64
end

struct EnsembleSpecsParameters
    burnin_years::Int64
    nyears::Int64
    annual_births_per_k::Int64
    ensemble_state_specification::StateParameters
    R_0::Float64
    gamma::Float64
    sigma::Float64
    target_Reff::Float64
    target_years::Int64
    min_vaccination_coverage::Float64
    max_vaccination_coverage::Float64
    nsims::Int64
end

function EnsembleSpecsParameters(;
        burnin_years::Int,
        nyears::Int,
        annual_births_per_k::Int64 = ANNUAL_BIRTHS_PER_K,
        ensemble_state_specification::StateParameters = StateParameters(
            500_000,
            Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95)
        ),
        R_0::Float64 = R0,
        gamma::Float64 = GAMMA,
        sigma::Float64 = SIGMA,
        target_Reff::Float64 = 0.9,
        target_years::Int = 2 * burnin_years,
        min_vaccination_coverage::Float64 = 0.6,
        max_vaccination_coverage::Float64 = 0.8,
        nsims::Int = 1000
    )
    @assert nyears >= target_years
    @assert min_vaccination_coverage < max_vaccination_coverage

    return EnsembleSpecsParameters(
        burnin_years,
        nyears,
        annual_births_per_k,
        ensemble_state_specification,
        R_0,
        gamma,
        sigma,
        target_Reff,
        target_years,
        min_vaccination_coverage,
        max_vaccination_coverage,
        nsims
    )
end


# end
