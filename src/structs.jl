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
    min_burnin_vaccination_coverage::Nothing,
    max_burnin_vaccination_coverage,
    min_vaccination_coverage,
    max_vaccination_coverage,
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
        min_vaccination_coverage,
        max_vaccination_coverage;
        max_burnin_vaccination_coverage = max_burnin_vaccination_coverage,
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
    min_vaccination_coverage,
    max_vaccination_coverage;
    max_burnin_vaccination_coverage = 1.0,
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
        max_burnin_vaccination_coverage,
        min_vaccination_coverage,
        max_vaccination_coverage,
    )
end

function calculate_min_burnin_vaccination_coverage(R_0; adjustment = 0.00)
    return min(1 - 1 / R_0 + adjustment, 1.0)
end

function DynamicsParameterSpecification(
    sigma::Float64,
    gamma::Float64,
    R_0::Float64;
    max_burnin_vaccination_coverage::Float64 = 1.0,
    max_vaccination_coverage::Float64 = 0.8,
    min_vaccination_coverage::Float64 = 0.0,
    kwargs...,
)
    kwargs_dict = Dict(kwargs)

    if haskey(kwargs_dict, :min_burnin_vaccination_coverage)
        min_burnin_vaccination_coverage = kwargs_dict[:min_burnin_vaccination_coverage]
    else
        if haskey(kwargs_dict, :burnin_adjustment)
            burnin_adjustment = kwargs_dict[:burnin_adjustment]
            min_burnin_vaccination_coverage = calculate_min_burnin_vaccination_coverage(
                R_0; adjustment = burnin_adjustment
            )
        else
            min_burnin_vaccination_coverage = calculate_min_burnin_vaccination_coverage(
                R_0
            )
        end
    end

    return DynamicsParameterSpecification(
        BETA_MEAN,
        BETA_FORCE,
        cos,
        sigma,
        gamma,
        MU,
        ANNUAL_BIRTHS_PER_K,
        EPSILON,
        R_0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_vaccination_coverage,
        max_vaccination_coverage,
    )
end

function DynamicsParameterSpecification(
    N::Int64,
    annual_births_per_k::Int64,
    beta_force::Float64,
    sigma::Float64,
    gamma::Float64,
    R_0::Float64,
    min_vaccination_coverage::Float64,
    max_vaccination_coverage::Float64;
    max_burnin_vaccination_coverage::Float64 = 1.0,
    kwargs...,
)
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
    epsilon = calculate_import_rate(mu, R_0, N)

    kwargs_dict = Dict(kwargs)

    if haskey(kwargs_dict, :min_burnin_vaccination_coverage)
        min_burnin_vaccination_coverage = kwargs_dict[:min_burnin_vaccination_coverage]
    else
        if haskey(kwargs_dict, :burnin_adjustment)
            burnin_adjustment = kwargs_dict[:burnin_adjustment]
            min_burnin_vaccination_coverage = calculate_min_burnin_vaccination_coverage(
                R_0; adjustment = burnin_adjustment
            )
        else
            min_burnin_vaccination_coverage = calculate_min_burnin_vaccination_coverage(
                R_0
            )
        end
    end

    return DynamicsParameterSpecification(
        beta_mean,
        beta_force,
        cos,
        sigma,
        gamma,
        mu,
        annual_births_per_k,
        epsilon,
        R_0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_vaccination_coverage,
        max_vaccination_coverage,
    )
end

function DynamicsParameterSpecification(
    N::Int64, annual_births_per_k::Int64, beta_force::Float64;
    max_burnin_vaccination_coverage::Float64 = 1.0,
    min_vaccination_coverage::Float64 = 0.0,
    max_vaccination_coverage::Float64 = 0.8,
    kwargs...,
)
    mu = calculate_mu(annual_births_per_k)
    beta_mean = calculate_beta(R0, GAMMA, mu, 1, N)
    epsilon = calculate_import_rate(mu, R0, N)

    kwargs_dict = Dict(kwargs)
    if haskey(kwargs_dict, :min_burnin_vaccination_coverage)
        min_burnin_vaccination_coverage = kwargs_dict[:min_burnin_vaccination_coverage]
    else
        if haskey(kwargs_dict, :burnin_adjustment)
            burnin_adjustment = kwargs_dict[:burnin_adjustment]
            min_burnin_vaccination_coverage = calculate_min_burnin_vaccination_coverage(
                R0; adjustment = burnin_adjustment
            )
        else
            min_burnin_vaccination_coverage = calculate_min_burnin_vaccination_coverage(
                R0
            )
        end
    end

    return DynamicsParameterSpecification(
        beta_mean,
        beta_force,
        cos,
        SIGMA,
        GAMMA,
        mu,
        annual_births_per_k,
        epsilon,
        R0,
        min_burnin_vaccination_coverage,
        max_burnin_vaccination_coverage,
        min_vaccination_coverage,
        max_vaccination_coverage,
    )
end

struct DynamicsParameters{
    T1<:AbstractFloat,
    T2<:Union{<:Integer,T1},
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
        ); digits = 2)

    vaccination_coverage = round(
        rand(
            Distributions.Uniform(
                dynamic_parameter_specification.min_vaccination_coverage,
                dynamic_parameter_specification.max_vaccination_coverage),
        ); digits = 2)

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

struct OutbreakThresholdChars{
    T1<:AbstractFloat,
    T2<:Integer,
    T3<:Vector{<:AbstractFloat},
    T4<:Vector{<:Integer},
    T5<:AbstractMatrix{<:Integer},
}
    daily_sensitivity::T1
    daily_specificity::T1
    daily_ppv::T1
    daily_npv::T1
    accuracy::T1
    matchedoutbreakbounds::T5
    noutbreaks::T2
    nalerts::T2
    detected_outbreak_size::T4
    missed_outbreak_size::T4
    n_true_outbreaks_detected::T2
    n_missed_outbreaks::T2
    n_correct_alerts::T2
    n_false_alerts::T2
    n_alerts_per_outbreak::T4
    period_sum_per_outbreak::T4
    perc_true_outbreaks_detected::T1
    perc_true_outbreaks_missed::T1
    falsealert_trueoutbreak_prop::T1
    correctalert_trueoutbreak_prop::T1
    trueoutbreak_alerts_prop::T1
    outbreaksmissed_alerts_prop::T1
    perc_alerts_false::T1
    perc_alerts_correct::T1
    detectiondelays::T4
    cases_before_alerts::T4
    cases_perc_before_alerts::T3
    cases_after_alerts::T4
    cases_perc_after_alerts::T3
    unavoidable_cases::T2
    avoidable_cases::T2
    n_outbreak_cases::T2
    n_tests::T2
    noise_rubella_prop::T1
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

struct EWSMetricSpecification{T1<:Integer,T2<:AbstractString}
    method::EWSMethod
    aggregation::T1
    bandwidth::T1
    lag::T1
    dirpath::T2
end

function EWSMetricSpecification(
    method::EWSMethod, aggregation::T1, bandwidth::T1, lag::T1
) where {T1<:Integer}
    return EWSMetricSpecification(
        method,
        aggregation,
        bandwidth,
        lag,
        joinpath(
            "ewsmethod_$(method_string(method))",
            "ewsaggregationdays_$(aggregation)",
            "ewsbandwidth_$(bandwidth)",
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

struct TestPositivity{T1<:AbstractArray{<:AbstractFloat}}
    one_day::T1
    seven_day::T1
    fourteen_day::T1
    thirty_day::T1
end

function TestPositivity(true_positive_vec, total_test_vec, alert_vec)
    return TestPositivity(
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 1
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 7
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 14
        ),
        calculate_test_positivity(
            true_positive_vec, total_test_vec, alert_vec, 30
        ),
    )
end

struct OptimalThresholdCharacteristics{
    T1<:StructVector{<:OutbreakThresholdChars},
    T2<:IndividualTestSpecification,
    T3<:NoiseSpecification,
    T4<:AbstractFloat,
    T5<:Integer,
}
    outbreak_threshold_chars::T1
    individual_test_specification::T2
    noise_specification::T3
    percent_clinic_tested::T4
    alert_threshold::T5
    accuracy::T4
end

# end
