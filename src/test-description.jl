export get_test_description

function get_test_description(test_specification::IndividualTestSpecification)
    test_specification == IndividualTestSpecification(1.0, 0.0, 0) &&
        return "Clinical case definition"

    test_specification.sensitivity == test_specification.specificity < 1.0 &&
        return "Imperfect Test ($(Int64(test_specification.sensitivity * 100))% Sensitive & Specific)"

    (
        test_specification.sensitivity < 1.0 ||
            test_specification.specificity < 1.0
    ) &&
        return "Imperfect Test ($(Int64(test_specification.sensitivity * 100))% Sensitive & $(Int64(test_specification.specificity * 100))% Specific)"

    test_specification.sensitivity == test_specification.specificity == 1.0 &&
        return "Perfect Test"

    return error("Don't have a description matching the test specification")
end
