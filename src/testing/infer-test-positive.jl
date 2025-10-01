function infer_true_positives(
        test_positivity_rate_vec,
        test_positive_vec,
        total_test_vec,
        total_clinic_visit_vec,
        test_result_lag,
    )
    inferred_positives = zeros(Float64, length(total_clinic_visit_vec))

    if test_result_lag == 0
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            test_positive_vec,
            total_test_vec,
            total_clinic_visit_vec,
        )
    else
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            total_clinic_visit_vec,
        )
    end

    return inferred_positives
end

function infer_true_positives!(
        inferred_positives,
        test_positivity_rate_vec,
        test_positive_vec,
        total_test_vec,
        total_clinic_visit_vec,
        test_result_lag,
    )
    if test_result_lag == 0
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            test_positive_vec,
            total_test_vec,
            total_clinic_visit_vec,
        )
    else
        infer_true_positives!(
            inferred_positives,
            test_positivity_rate_vec,
            total_clinic_visit_vec,
        )
    end

    return nothing
end
function infer_true_positives!(
        inferred_positives,
        test_positivity_rate_vec,
        test_positive_vec,
        total_test_vec,
        total_clinic_visit_vec,
    )
    @. inferred_positives =
        test_positive_vec +
        (total_clinic_visit_vec - total_test_vec) * test_positivity_rate_vec

    return nothing
end

function infer_true_positives!(
        inferred_positives,
        test_positivity_rate_vec,
        total_clinic_visit_vec,
    )
    return @. inferred_positives = total_clinic_visit_vec * test_positivity_rate_vec
end
