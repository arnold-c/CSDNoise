export OutbreakSpecification,
    AlertMethod,
    OutbreakDetectionSpecification,
    AbstractThresholds,
    Thresholds,
    OutbreakThresholds

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

# TODO: Update to use sumtype and add implementation to testing vecs creation
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
