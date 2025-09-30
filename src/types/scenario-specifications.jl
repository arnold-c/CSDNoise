export ScenarioSpecification

# TODO: Doesn't seem to be used
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
