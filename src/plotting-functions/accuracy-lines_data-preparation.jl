export prepare_line_plot_data

function prepare_line_plot_data(
        results::StructVector{OptimizationResult},
        ewsmetric = "mean",
        tests = [
            IndividualTestSpecification(0.8, 0.8, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
            IndividualTestSpecification(0.95, 0.95, 0),
            IndividualTestSpecification(0.96, 0.96, 0),
            IndividualTestSpecification(0.97, 0.97, 0),
            IndividualTestSpecification(0.98, 0.98, 0),
            IndividualTestSpecification(0.99, 0.99, 0),
            IndividualTestSpecification(1.0, 1.0, 0),
        ];
        tiebreaker_preference = "specificity",
    )
    tiebreaker_field = if tiebreaker_preference == "speed"
        :consecutive_thresholds
    elseif tiebreaker_preference == "specificity"
        :specificity
    else
        error(
            "Invalid preference: $tiebreaker_preference. Please choose either \"speed\" or \"specificity\"."
        )
    end
    tiebreaker_rev = tiebreaker_preference == "specificity"

    metric_mask = results.ews_metric .== ewsmetric
    test_mask = [t in tests for t in results.test_specification]
    combined_mask = metric_mask .& test_mask
    filtered_results = results[combined_mask]

    groups = group_structvector(
        filtered_results,
        :noise_level,
        :noise_type_description,
        :test_specification
    )

    selected_results_vec = OptimizationResult[]
    for (key, group) in groups
        tiebreaker_values = getproperty(group, tiebreaker_field)
        if tiebreaker_rev
            best_idx = argmax(tiebreaker_values)
        else
            best_idx = argmin(tiebreaker_values)
        end
        push!(selected_results_vec, group[best_idx])
    end

    selected_results = StructVector(selected_results_vec)

    sort_order = sortperm(
        selected_results;
        by = r -> (
            -r.test_specification.sensitivity,
            -r.test_specification.specificity,
            r.test_specification.test_result_lag,
        )
    )

    return selected_results[sort_order]
end
