export OutbreakSpecification,
    AlertMethod,
    OutbreakDetectionSpecification,
    AbstractThresholds,
    Thresholds,
    OutbreakThresholds

"""
    OutbreakSpecification

Configuration for defining outbreak detection criteria and thresholds.

This struct specifies the parameters used to identify and classify disease outbreaks
in epidemiological simulations. It defines the minimum conditions that must be met
for a period to be considered an outbreak, including threshold values, duration
requirements, and size constraints.

The `dirpath` field is automatically constructed from the parameter values to create
a unique filesystem path for storing outbreak-related results and analyses.

# Fields
- `outbreak_threshold::Int64`: Minimum number of cases required to trigger outbreak detection
- `minimum_outbreak_duration::Int64`: Minimum number of time periods an outbreak must persist
- `minimum_outbreak_size::Int64`: Minimum total number of cases required for outbreak classification
- `dirpath::String`: Automatically generated directory path for storing outbreak analysis results

# Constructor
    OutbreakSpecification(; outbreak_threshold, minimum_outbreak_duration, minimum_outbreak_size, dirpath)

# Example
```julia
# Create outbreak specification with specific thresholds
outbreak_spec = OutbreakSpecification(
    outbreak_threshold = 10,
    minimum_outbreak_duration = 5,
    minimum_outbreak_size = 50,
    dirpath = "min_outbreak_dur_5/min_outbreak_size_50/outbreak_threshold_10"
)

# Check outbreak criteria
println("Outbreak threshold: \$(outbreak_spec.outbreak_threshold) cases")
println("Minimum duration: \$(outbreak_spec.minimum_outbreak_duration) periods")
```

# See Also
- [`OutbreakDetectionSpecification`](@ref): Alert system configuration for outbreak detection
- [`OutbreakThresholds`](@ref): Threshold specifications for outbreak analysis
"""
Base.@kwdef struct OutbreakSpecification
    outbreak_threshold::Int64
    minimum_outbreak_duration::Int64
    minimum_outbreak_size::Int64
    dirpath::String
end

"""
    OutbreakSpecification(outbreak_threshold, minimum_outbreak_duration, minimum_outbreak_size)

Construct an OutbreakSpecification with automatic directory path generation.

This constructor creates an outbreak specification by taking the core parameters and
automatically generating a hierarchical directory path that uniquely identifies this
outbreak configuration. The directory path incorporates all parameter values to ensure
unique storage locations for different outbreak detection scenarios.

The generated directory structure follows this hierarchy:
- `min_outbreak_dur_\$(minimum_outbreak_duration)`
- `min_outbreak_size_\$(minimum_outbreak_size)`
- `outbreak_threshold_\$(outbreak_threshold)`

# Arguments
- `outbreak_threshold`: Minimum number of cases required to trigger outbreak detection
- `minimum_outbreak_duration`: Minimum number of time periods an outbreak must persist
- `minimum_outbreak_size`: Minimum total number of cases required for outbreak classification

# Returns
- `OutbreakSpecification`: Complete outbreak specification with auto-generated directory path

# Example
```julia
# Create outbreak specification with automatic path generation
outbreak_spec = OutbreakSpecification(10, 5, 50)

# The dirpath will be: "min_outbreak_dur_5/min_outbreak_size_50/outbreak_threshold_10"
println(outbreak_spec.dirpath)
```
"""
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
"""
    AlertMethod

Specification for outbreak alert detection methods.

This struct defines the method used for detecting outbreak alerts in disease surveillance
systems. It validates that the specified method is one of the supported alert detection
algorithms and provides a type-safe way to specify detection strategies.

Currently supported methods include daily threshold-based detection, moving average
approaches, and hybrid methods that combine multiple detection strategies.

# Fields
- `method_name::String`: Name of the alert detection method (validated against available methods)

# Constructor
    AlertMethod(method_name::String)

The constructor validates that the method name is one of the supported options:
- `"dailythreshold"`: Simple daily case threshold detection
- `"movingavg"`: Moving average-based detection
- `"dailythreshold_movingavg"`: Combined daily threshold and moving average
- `"inferred_movingavg"`: Inferred moving average detection

# Example
```julia
# Create alert method specification
alert_method = AlertMethod("dailythreshold")

# This will throw an error for invalid methods
try
    invalid_method = AlertMethod("invalid_method")
catch e
    println("Error: \$e")
end
```

# See Also
- [`OutbreakDetectionSpecification`](@ref): Uses AlertMethod for outbreak detection configuration
"""
Base.@kwdef struct AlertMethod
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

"""
    OutbreakDetectionSpecification

Configuration for outbreak detection and alert systems in disease surveillance.

This struct defines all parameters needed for outbreak detection including alert thresholds,
testing protocols, and detection methods. It combines epidemiological surveillance parameters
with diagnostic testing specifications to create a complete outbreak detection system.

The `dirpath` field is automatically constructed from the parameter values to create a unique
filesystem path for storing detection-related results and analyses.

# Fields
- `alert_threshold::Int64`: Number of cases required to trigger an outbreak alert
- `moving_average_lag::Int64`: Time lag (in days) for moving average calculations
- `percent_visit_clinic::Float64`: Proportion of infected individuals who visit a clinic (0.0-1.0)
- `percent_clinic_tested::Float64`: Proportion of clinic visitors who receive diagnostic testing (0.0-1.0)
- `percent_tested::Float64`: Overall proportion of infected individuals tested (calculated as percent_visit_clinic × percent_clinic_tested)
- `alert_method::AlertMethod`: Method used for outbreak alert detection
- `dirpath::String`: Automatically generated directory path for storing detection analysis results

# Constructor
    OutbreakDetectionSpecification(; alert_threshold, moving_average_lag, percent_visit_clinic,
                                   percent_clinic_tested, percent_tested, alert_method, dirpath)

# Example
```julia
# Create outbreak detection specification
detection_spec = OutbreakDetectionSpecification(
    alert_threshold = 5,
    moving_average_lag = 7,
    percent_visit_clinic = 0.8,
    percent_clinic_tested = 0.9,
    percent_tested = 0.72,  # 0.8 × 0.9
    alert_method = AlertMethod("movingavg"),
    dirpath = "alertmethod_movingavg/alertthreshold_5/moveavglag_7/perc_visit_clinic_0.8/perc_clinic_tested_0.9"
)
```

# See Also
- [`AlertMethod`](@ref): Alert detection method specification
- [`OutbreakSpecification`](@ref): Outbreak definition parameters
"""
Base.@kwdef struct OutbreakDetectionSpecification
    alert_threshold::Int64
    moving_average_lag::Int64
    percent_visit_clinic::Float64
    percent_clinic_tested::Float64
    percent_tested::Float64
    alert_method::AlertMethod
    dirpath::String
end

"""
    OutbreakDetectionSpecification(alert_threshold, moving_average_lag, percent_visit_clinic,
                                  percent_clinic_tested, alert_method)

Construct an OutbreakDetectionSpecification with automatic directory path generation.

This constructor creates a complete outbreak detection specification by taking the core
parameters and automatically generating a hierarchical directory path and calculating
derived values. The directory path incorporates all relevant parameter values to ensure
unique storage locations for different detection scenarios.

The constructor automatically:
- Calculates `percent_tested` as the product of `percent_visit_clinic` and `percent_clinic_tested`
- Generates a hierarchical directory path based on the alert method and parameters
- Creates an `AlertMethod` object from the method string

The generated directory structure varies by alert method:
- For "dailythreshold": `alertmethod_\$(alert_method)/alertthreshold_\$(alert_threshold)/\$(testingdirpath)`
- For other methods: `alertmethod_\$(alert_method)/alertthreshold_\$(alert_threshold)/moveavglag_\$(moving_average_lag)/\$(testingdirpath)`

Where `testingdirpath` is: `perc_visit_clinic_\$(percent_visit_clinic)/perc_clinic_tested_\$(percent_clinic_tested)`

# Arguments
- `alert_threshold`: Number of cases required to trigger an outbreak alert
- `moving_average_lag`: Time lag (in days) for moving average calculations
- `percent_visit_clinic`: Proportion of infected individuals who visit a clinic (0.0-1.0)
- `percent_clinic_tested`: Proportion of clinic visitors who receive diagnostic testing (0.0-1.0)
- `alert_method`: String name of the alert detection method

# Returns
- `OutbreakDetectionSpecification`: Complete detection specification with auto-generated values

# Example
```julia
# Create detection specification with automatic calculations
detection_spec = OutbreakDetectionSpecification(
    5,           # alert_threshold
    7,           # moving_average_lag
    0.8,         # percent_visit_clinic
    0.9,         # percent_clinic_tested
    "movingavg"  # alert_method
)

# percent_tested will be automatically calculated as 0.72 (0.8 × 0.9)
println("Overall testing rate: \$(detection_spec.percent_tested)")
```
"""
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

"""
    AbstractThresholds

Abstract base type for threshold specifications in epidemiological surveillance systems.

This abstract type serves as the parent for all threshold-related specifications
used in outbreak detection and effective reproduction number (Reff) analysis. It provides
a common interface for different types of threshold configurations while allowing for
specialized implementations with different field requirements.

Thresholds can represent:
- Outbreak detection: Periods when incidence exceeds a specified threshold
- Reff analysis: Periods when the effective reproduction number (Reff) is >= 1

# Subtypes
- [`Thresholds`](@ref): Basic threshold specification with bounds and duration
- [`OutbreakThresholds`](@ref): Extended threshold specification including infection counts

# See Also
- [`Thresholds`](@ref): Basic threshold implementation
- [`OutbreakThresholds`](@ref): Extended threshold implementation with infection tracking
- [`Reff_ge_than_one`](@ref): Function for calculating Reff thresholds
"""
abstract type AbstractThresholds end

"""
    Thresholds <: AbstractThresholds

Basic threshold specification for epidemiological surveillance analysis.

This struct defines threshold parameters used in outbreak detection and effective reproduction
number (Reff) analysis, specifying the bounds and duration criteria for identifying significant
epidemiological events. It provides the fundamental threshold configuration for time-series
analysis of disease surveillance data.

The thresholds can represent:
- Outbreak periods: Time intervals when incidence exceeds a specified threshold
- Reff periods: Time intervals when the effective reproduction number (Reff) is >= 1

# Fields
- `lower_bounds::Vector{Int64}`: Starting time indices of threshold exceedance periods
- `upper_bounds::Vector{Int64}`: Ending time indices of threshold exceedance periods
- `duration::Vector{Int64}`: Duration (in time periods) of each threshold exceedance

# Constructor
    Thresholds(; lower_bounds, upper_bounds, duration)

# Example
```julia
# Create threshold specification for outbreak detection
outbreak_thresholds = Thresholds(
    lower_bounds = [5, 10, 15],
    upper_bounds = [20, 30, 40],
    duration = [16, 21, 26]  # upper - lower + 1
)

# Create threshold specification for Reff >= 1 periods
reff_thresholds = Reff_ge_than_one(Reff_vec)

# Access threshold parameters
println("Lower bounds: \$(outbreak_thresholds.lower_bounds)")
println("Duration requirements: \$(outbreak_thresholds.duration)")
```

# See Also
- [`AbstractThresholds`](@ref): Parent abstract type
- [`OutbreakThresholds`](@ref): Extended threshold specification with infection counts
- [`Reff_ge_than_one`](@ref): Function for calculating Reff >= 1 thresholds
"""
Base.@kwdef struct Thresholds <: AbstractThresholds
    lower_bounds::Vector{Int64}
    upper_bounds::Vector{Int64}
    duration::Vector{Int64}
end

"""
    OutbreakThresholds <: AbstractThresholds

Extended threshold specification for outbreak detection with infection count tracking.

This struct extends the basic threshold specification to include tracking of infection
counts during threshold exceedance periods. It provides comprehensive threshold
configuration for detailed outbreak analysis that requires monitoring both detection
criteria and the epidemiological burden during outbreak periods.

Unlike the basic `Thresholds` type, this struct is specifically designed for outbreak
analysis (incidence-based thresholds) and includes the total number of infections
that occurred during each threshold exceedance period. This is used to classify
whether a threshold exceedance qualifies as a true outbreak based on minimum
duration and size criteria.

# Fields
- `lower_bounds::Vector{Int64}`: Starting time indices of outbreak periods
- `upper_bounds::Vector{Int64}`: Ending time indices of outbreak periods
- `duration::Vector{Int64}`: Duration (in time periods) of each outbreak
- `num_infections_during_bounds::Vector{Int64}`: Total number of infections during each outbreak period

# Constructor
    OutbreakThresholds(; lower_bounds, upper_bounds, duration, num_infections_during_bounds)

# Example
```julia
# Create extended threshold specification with infection tracking
outbreak_thresholds = OutbreakThresholds(
    lower_bounds = [5, 10, 15],
    upper_bounds = [20, 30, 40],
    duration = [16, 21, 26],
    num_infections_during_bounds = [50, 150, 300]
)

# Access infection count data
println("Infections during bounds: \$(outbreak_thresholds.num_infections_during_bounds)")
println("Total outbreaks detected: \$(length(outbreak_thresholds.lower_bounds))")
```

# See Also
- [`AbstractThresholds`](@ref): Parent abstract type
- [`Thresholds`](@ref): Basic threshold specification without infection tracking
- [`OutbreakSpecification`](@ref): Outbreak definition parameters
- [`calculate_outbreak_thresholds`](@ref): Function for calculating outbreak thresholds
"""
Base.@kwdef struct OutbreakThresholds <: AbstractThresholds
    lower_bounds::Vector{Int64}
    upper_bounds::Vector{Int64}
    duration::Vector{Int64}
    num_infections_during_bounds::Vector{Int64}
end
