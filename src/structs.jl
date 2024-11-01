# module ODStructs
#
# export SimTimeParameters, EnsembleSpecification, DynamicsParameters,
#     StateParameters, OutbreakThresholdChars, OutbreakDetectionSpecification,
#     OutbreakSpecification, IndividualTestSpecification, NoiseSpecification

using StaticArrays
using LabelledArrays
using StructArrays
using Match
using SumTypes
using Distributions: Distributions
using Random: Random
using UnPack: @unpack
using Dates: Dates

# include("transmission-functions.jl")
# using .TransmissionFunctions

struct SimTimeParameters{
    T1<:AbstractFloat,
    T2<:StepRangeLen,
    T3<:Tuple{T1,T1},
    T4<:Int,
}
    burnin::T1
    tmin::T1
    tmax::T1
    tstep::T1
    trange::T2
    tspan::T3
    tlength::T4
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

struct DynamicsParameterSpecification{
    T1<:AbstractFloat,T2<:Union{<:Integer,T1},T3<:Function
}
    beta_mean::T1
    beta_force::T1
    seasonality::T3
    sigma::T1
    gamma::T1
    mu::T1
    annual_births_per_k::T2
    epsilon::T1
    R_0::T1
    min_burnin_vaccination_coverage::T1
    max_burnin_vaccination_coverage::T1
    min_vaccination_coverage::T1
    max_vaccination_coverage::T1
end

function DynamicsParameterSpecification(
    beta_mean,
    beta_force,
    seasonality,
    sigma,
    gamma,
    mu,
    annual_births_per_k,
    epsilon,
    R_0,
    burnin_vaccination_params::Tuple{<:Nothing,<:AbstractFloat},
    vaccination_params::Tuple{<:AbstractFloat,<:AbstractFloat},
)
    return DynamicsParameterSpecification(
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        calculate_min_burnin_vaccination_coverage(R_0),
        burnin_vaccination_params[2],
        vaccination_params[1],
        vaccination_params[2],
    )
end

function calculate_min_burnin_vaccination_coverage(R_0; adjustment = 0.00)
    return min(calculate_min_burnin_vaccination_coverage(R_0) + adjustment, 1.0)
end

function calculate_herd_immunity_threshold(R_0)
    return (1 - 1 / R_0)
end

function DynamicsParameterSpecification(
    beta_mean,
    beta_force,
    seasonality,
    sigma,
    gamma,
    mu,
    annual_births_per_k,
    epsilon,
    R_0,
    burnin_vaccination_params::Tuple{
        <:Nothing,<:AbstractFloat,<:AbstractFloat,<:Integer
    },
    vaccination_params::Union{
        Tuple{<:AbstractFloat,<:AbstractFloat},Tuple{<:Nothing,<:Nothing}
    },
    initial_states,
)
    return DynamicsParameterSpecification(
        beta_mean,
        beta_force,
        seasonality,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        calculate_vaccination_rate_to_achieve_Reff(
            burnin_vaccination_params[3],
            burnin_vaccination_params[4],
            initial_states.S,
            initial_states.N,
            R_0,
            mu,
        ),
        burnin_vaccination_params[2],
        vaccination_params[1],
        vaccination_params[2],
    )
end

function DynamicsParameterSpecification(
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
    min_vaccination_coverage::T1,
    max_vaccination_coverage::T1,
) where {T1<:Nothing}
    return DynamicsParameterSpecification(
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
            "Target Reff cannot be reached in the burn-in period. Initial Reff = $(R_0 * S/N). Try a longer burn-in period or a smaller target Reff."
        )
    end
    return round(vaccination_coverage; digits = 4)
end

struct DynamicsParameters{
    T1<:AbstractFloat,
    T2<:Real,
    T3<:Function,
}
    beta_mean::T1
    beta_force::T1
    seasonality::T3
    sigma::T1
    gamma::T1
    mu::T1
    annual_births_per_k::T2
    epsilon::T1
    R_0::T1
    min_burnin_vaccination_coverage::T1
    max_burnin_vaccination_coverage::T1
    min_vaccination_coverage::T1
    max_vaccination_coverage::T1
    burnin_vaccination_coverage::T1
    vaccination_coverage::T1
end

function DynamicsParameters(
    dynamic_parameter_specification::T1; seed = 1234
) where {T1<:DynamicsParameterSpecification}
    Random.seed!(seed)

    burnin_vaccination_coverage = round(
        rand(
            Distributions.Uniform(
                dynamic_parameter_specification.min_burnin_vaccination_coverage,
                dynamic_parameter_specification.max_burnin_vaccination_coverage,
            ),
        ); digits = 4)

    vaccination_coverage =
        if dynamic_parameter_specification.min_burnin_vaccination_coverage ==
           dynamic_parameter_specification.min_vaccination_coverage &&
            dynamic_parameter_specification.max_burnin_vaccination_coverage ==
           dynamic_parameter_specification.max_vaccination_coverage
            burnin_vaccination_coverage
        else
            round(
                rand(
                    Distributions.Uniform(
                        dynamic_parameter_specification.min_vaccination_coverage,
                        dynamic_parameter_specification.max_vaccination_coverage,
                    ),
                ); digits = 4)
        end

    return DynamicsParameters(
        [
            getfield(dynamic_parameter_specification, f) for
            f in fieldnames(DynamicsParameterSpecification)
        ]
        ...,
        burnin_vaccination_coverage,
        vaccination_coverage,
    )
end

struct StateParameters{T1<:SLArray,T2<:SLArray}
    init_states::T1
    init_state_props::T2
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

struct EnsembleSpecification{
    T1<:Tuple,
    T2<:StateParameters,
    T3<:DynamicsParameterSpecification,
    T4<:SimTimeParameters,
    T5<:Integer,
    T6<:AbstractString,
}
    modeltypes::T1
    state_parameters::T2
    dynamics_parameter_specification::T3
    time_parameters::T4
    nsims::T5
    dirpath::T6
end

function EnsembleSpecification(
    modeltypes::Tuple,
    state_parameters::StateParameters,
    dynamics_parameter_specification::DynamicsParameterSpecification,
    time_parameters::SimTimeParameters,
    nsims::Int64,
)
    dirpath = outdir(
        "ensemble",
        modeltypes...,
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
        modeltypes,
        state_parameters,
        dynamics_parameter_specification,
        time_parameters,
        nsims,
        dirpath,
    )
end

struct OutbreakSpecification{T1<:Integer,T2<:AbstractString}
    outbreak_threshold::T1
    minimum_outbreak_duration::T1
    minimum_outbreak_size::T1
    dirpath::T2
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

struct AlertMethod{T1<:AbstractString}
    method_name::T1
    function AlertMethod(method_name::T1) where {T1<:AbstractString}
        available_test_methods = [
            "dailythreshold", "movingavg", "dailythreshold_movingavg",
            "inferred_movingavg",
        ]
        if !in(method_name, available_test_methods)
            error(
                "$(method_name) is not a valid test method. It must be one of $(available_test_methods)"
            )
        end
        return new{T1}(method_name)
    end
end

struct OutbreakDetectionSpecification{
    T1<:Integer,T2<:AbstractFloat,T3<:AlertMethod,T4<:AbstractString
}
    alert_threshold::T1
    moving_average_lag::T1
    percent_visit_clinic::T2
    percent_clinic_tested::T2
    percent_tested::T2
    alert_method::T3
    dirpath::T4
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

    dirpath = @match alert_method begin
        "dailythreshold" => joinpath(
            alertdirpath,
            testingdirpath,
        )
        _ => joinpath(
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

struct IndividualTestSpecification{T1<:AbstractFloat,T2<:Integer}
    sensitivity::T1
    specificity::T1
    test_result_lag::T2
end

function get_test_description(test_specification::IndividualTestSpecification)
    description = Match.@match test_specification begin
        IndividualTestSpecification(1.0, 0.0, 0) => "Clinical case definition"
        IndividualTestSpecification(x::AbstractFloat, x::AbstractFloat, 0) where {x<1.0} => "RDT-like ($(test_specification.sensitivity * 100)% sens/spec)"
        IndividualTestSpecification(1.0, 1.0, x::Int) => "ELISA-like ($(x) day lag)"
    end
    return description
end

abstract type NoiseSpecification end

struct PoissonNoiseSpecification{
    T1<:AbstractString,T2<:AbstractFloat
} <: NoiseSpecification
    noise_type::T1
    noise_mean_scaling::T2
end

function PoissonNoiseSpecification(
    noise_mean_scaling::T
) where {T<:AbstractFloat}
    return PoissonNoiseSpecification("poisson", noise_mean_scaling)
end

struct DynamicalNoiseSpecification{
    T1<:AbstractString,T2<:AbstractFloat,T3<:Integer
} <: NoiseSpecification
    noise_type::T1
    R_0::T2
    latent_period::T3
    duration_infection::T3
    correlation::T1
    noise_mean_scaling::T2
end

function DynamicalNoiseSpecification(
    R_0::T2,
    latent_period::T3,
    duration_infection::T3,
    correlation::T1,
    noise_mean_scaling::T2,
) where {T1<:AbstractString,T2<:AbstractFloat,T3<:Integer}
    return DynamicalNoiseSpecification(
        "dynamical",
        R_0,
        latent_period,
        duration_infection,
        correlation,
        noise_mean_scaling,
    )
end

function get_noise_description(
    noise_specification::T
) where {T<:NoiseSpecification}
    return noise_specification.noise_type
end

function get_noise_description(noise_specification::DynamicalNoiseSpecification)
    return string(
        noise_specification.noise_type, ", ", noise_specification.correlation
    )
end

function get_noise_magnitude_description(
    noise_specification::T
) where {T<:NoiseSpecification}
    return string("Poisson scaling: ", noise_specification.noise_mean_scaling)
end

function get_noise_magnitude(
    noise_specification::T
) where {T<:NoiseSpecification}
    return noise_specification.noise_mean_scaling
end

# function get_noise_magnitude(
#     noise_specification::DynamicalNoiseSpecification
# )
#     return string("Rubella vax: ", noise_specification.vaccination_coverage)
# end

function getdirpath(spec::NoiseSpecification)
    return reduce(
        joinpath,
        map(
            p -> "$(p)_$(getproperty(spec, p))",
            propertynames(spec),
        ),
    )
end

@sum_type EWSMethod begin
    Backward
    Centered
end

struct EWSMetricSpecification{T1<:Dates.Day,T2<:Integer,T3<:AbstractString}
    method::EWSMethod
    aggregation::T1
    bandwidth::T1
    lag::T2
    dirpath::T3
end

function EWSMetricSpecification(
    method::EWSMethod,
    aggregation::Dates.Day,
    bandwidth::Dates.Day,
    lag::T1,
) where {T1<:Integer}
    aggregation_days_val = Dates.value(aggregation)
    bandwidth_days_val = Dates.value(bandwidth)

    return EWSMetricSpecification(
        method,
        aggregation,
        bandwidth,
        lag,
        joinpath(
            "ewsmethod_$(method_string(method))",
            "ewsaggregationdays_$(aggregation_days_val)",
            "ewsbandwidth_$(bandwidth_days_val)",
            "ewslag_$(lag)",
        ),
    )
end

function EWSMetricSpecification(
    method::EWSMethod, aggregation::T1, bandwidth::T2, lag::T3
) where {
    T1<:Dates.DatePeriod,
    T2<:Dates.DatePeriod,
    T3<:Integer,
}
    aggregation_days_val = Dates.days(aggregation)
    bandwidth_days_val = Dates.days(bandwidth)
    aggregation_days = Dates.Day(aggregation_days_val)
    bandwidth_days = Dates.Day(bandwidth_days_val)

    return EWSMetricSpecification(
        method,
        aggregation_days,
        bandwidth_days,
        lag,
        joinpath(
            "ewsmethod_$(method_string(method))",
            "ewsaggregationdays_$(aggregation_days_val)",
            "ewsbandwidth_$(bandwidth_days_val)",
            "ewslag_$(lag)",
        ),
    )
end

method_string(method::EWSMethod) = lowercase(split(string(method), "::")[1])

struct EWSMetrics{
    T1<:EWSMetricSpecification,
    T2<:AbstractFloat,
    T3<:AbstractArray{T2},
}
    ews_specification::T1
    mean::T3
    variance::T3
    coefficient_of_variation::T3
    index_of_dispersion::T3
    skewness::T3
    kurtosis::T3
    autocovariance::T3
    autocorrelation::T3
    mean_tau::T2
    variance_tau::T2
    coefficient_of_variation_tau::T2
    index_of_dispersion_tau::T2
    skewness_tau::T2
    kurtosis_tau::T2
    autocovariance_tau::T2
    autocorrelation_tau::T2
end

struct ScenarioSpecification{
    T1<:EnsembleSpecification,
    T2<:OutbreakSpecification,
    T3<:NoiseSpecification,
    T4<:OutbreakDetectionSpecification,
    T5<:IndividualTestSpecification,
    T6<:EWSMetricSpecification,
    T7<:AbstractString,
}
    ensemble_specification::T1
    outbreak_specification::T2
    noise_specification::T3
    outbreak_detection_specification::T4
    individual_test_specification::T5
    ewsmetric_specification::T6
    dirpath::T7
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

# end
