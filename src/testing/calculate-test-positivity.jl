function calculate_test_positivity_rate(
        test_positive_vec, total_test_vec, agg_days
    )
    test_positivity_rate_vec = zeros(Float64, length(total_test_vec))

    calculate_test_positivity_rate!(
        test_positivity_rate_vec, test_positive_vec, total_test_vec, agg_days
    )

    return test_positivity_rate_vec
end

function calculate_test_positivity_rate!(
        test_positivity_rate_vec, test_positive_vec, total_test_vec, agg_days
    )
    for i in axes(test_positivity_rate_vec, 1)
        if i < agg_days
            test_positivity_rate_vec[i] =
                sum(test_positive_vec[begin:i]) / sum(total_test_vec[begin:i])
            continue
        end

        test_positivity_rate_vec[i] =
            sum(test_positive_vec[(i - agg_days + 1):i]) /
            sum(total_test_vec[(i - agg_days + 1):i])
    end
    return
end

function calculate_test_positivity(
        true_positive_vec, total_test_vec, alert_vec, agg_days
    )
    @views outvec = zeros(Float64, length(true_positive_vec) ÷ agg_days, 2)
    @inbounds for i in axes(outvec, 1)
        start_ind = 1 + (i - 1) * agg_days
        end_ind = start_ind + (agg_days - 1)

        @views total_test_sum = sum(total_test_vec[start_ind:end_ind])
        @views true_positive_sum = sum(true_positive_vec[start_ind:end_ind])
        @views num_outbreak_days = sum(alert_vec[start_ind:end_ind])
        agg_outbreak_status = num_outbreak_days >= agg_days / 2 ? 1 : 0

        outvec[i, 1] = true_positive_sum / total_test_sum
        outvec[i, 2] = agg_outbreak_status
    end
    return outvec
end
