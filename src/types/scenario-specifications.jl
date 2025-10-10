export ScenarioSpecification

"""
    ScenarioSpecification

Complete specification for a disease surveillance scenario combining all modeling components.

This struct aggregates all the individual specifications needed to define a complete
epidemiological surveillance scenario, including disease dynamics, outbreak detection,
diagnostic testing, early warning systems, and noise characteristics. It serves as the
top-level configuration object that encapsulates all parameters needed for simulation
and analysis.

The `dirpath` field is automatically constructed from the component specifications to
create a unique filesystem path for storing results associated with this scenario.

# Fields
- `ensemble_specification::EnsembleSpecification`: Disease dynamics and simulation ensemble parameters
- `outbreak_specification::OutbreakSpecification`: Outbreak definition and detection thresholds
- `noise_specification::NoiseSpecification`: Noise model specification (Poisson or dynamical)
- `outbreak_detection_specification::OutbreakDetectionSpecification`: Alert system configuration
- `individual_test_specification::IndividualTestSpecification`: Diagnostic test characteristics
- `ewsmetric_specification::EWSMetricSpecification`: Early warning system metric parameters
- `dirpath::String`: Automatically generated directory path for scenario results

# Constructor
    ScenarioSpecification(; ensemble_specification, outbreak_specification, noise_specification,
                         outbreak_detection_specification, individual_test_specification,
                         ewsmetric_specification)

# Example
```julia
# Create a complete scenario specification
scenario = ScenarioSpecification(
    ensemble_specification = ensemble_spec,
    outbreak_specification = outbreak_spec,
    noise_specification = PoissonNoise(noise_mean_scaling = 0.1),
    outbreak_detection_specification = detection_spec,
    individual_test_specification = test_spec,
    ewsmetric_specification = ews_spec
)

# Access the auto-generated directory path
println(scenario.dirpath)  # Shows hierarchical path based on all component parameters
```

# See Also
- [`EnsembleSpecification`](@ref): Disease dynamics and simulation parameters
- [`OutbreakSpecification`](@ref): Outbreak definition parameters
- [`NoiseSpecification`](@ref): Noise model specifications
- [`OutbreakDetectionSpecification`](@ref): Alert system configuration
- [`IndividualTestSpecification`](@ref): Diagnostic test parameters
- [`EWSMetricSpecification`](@ref): Early warning system parameters
"""
# TODO: Doesn't seem to be used
Base.@kwdef struct ScenarioSpecification
    ensemble_specification::EnsembleSpecification
    outbreak_specification::OutbreakSpecification
    noise_specification::NoiseSpecification
    outbreak_detection_specification::OutbreakDetectionSpecification
    individual_test_specification::IndividualTestSpecification
    ewsmetric_specification::EWSMetricSpecification
    dirpath::String
end

"""
    ScenarioSpecification(ensemble_specification, outbreak_specification, noise_specification,
                         outbreak_detection_specification, individual_test_specification,
                         ewsmetric_specification)

Construct a ScenarioSpecification with automatic directory path generation.

This constructor creates a complete scenario specification by combining all the individual
component specifications and automatically generating a hierarchical directory path that
uniquely identifies this scenario configuration. The directory path is constructed by
joining the directory paths from each component specification along with test-specific
parameters.

The generated directory structure follows this hierarchy:
- Ensemble specification path (disease dynamics, population, simulation parameters)
- Outbreak specification path (outbreak thresholds and detection criteria)
- Noise specification path (noise model type and parameters)
- Outbreak detection specification path (alert methods and testing coverage)
- Test sensitivity, specificity, and lag parameters
- EWS metric specification path (early warning system configuration)

# Arguments
- `ensemble_specification::EnsembleSpecification`: Disease dynamics and simulation ensemble parameters
- `outbreak_specification::OutbreakSpecification`: Outbreak definition and detection thresholds
- `noise_specification::NoiseSpecification`: Noise model specification (Poisson or dynamical)
- `outbreak_detection_specification::OutbreakDetectionSpecification`: Alert system configuration
- `individual_test_specification::IndividualTestSpecification`: Diagnostic test characteristics
- `ewsmetric_specification::EWSMetricSpecification`: Early warning system metric parameters

# Returns
- `ScenarioSpecification`: Complete scenario specification with auto-generated directory path

# Example
```julia
# Create component specifications
ensemble_spec = EnsembleSpecification(...)
outbreak_spec = OutbreakSpecification(...)
noise_spec = PoissonNoise(noise_mean_scaling = 0.1)
detection_spec = OutbreakDetectionSpecification(...)
test_spec = IndividualTestSpecification(sensitivity = 0.95, specificity = 0.99, test_result_lag = 2)
ews_spec = EWSMetricSpecification(...)

# Create complete scenario
scenario = ScenarioSpecification(
    ensemble_spec, outbreak_spec, noise_spec,
    detection_spec, test_spec, ews_spec
)

# The dirpath will be something like:
# "ensemble/.../outbreak/.../noise/.../detection/.../testsens_0.95/testspec_0.99/testlag_2/ews/..."
```

# Notes
- The directory path is deterministic and uniquely identifies the scenario configuration
- Test parameters are explicitly included in the path for easy identification
- The path structure facilitates organized storage and retrieval of simulation results
"""
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
