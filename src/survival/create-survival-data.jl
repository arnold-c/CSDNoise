using DataFrames: nrow

export create_ews_survival_data

function create_ews_survival_data(
        ews_optimal_simulation_df
    )
    nsims = nrow(ews_optimal_simulation_df)
    unique_detection_indices = sort(
        filter(
            !isnothing, unique(ews_optimal_simulation_df.detection_index)
        ),
    )

    unique_null_detection_indices = sort(
        filter(
            !isnothing,
            unique(
                ews_optimal_simulation_df.null_detection_index
            ),
        ),
    )

    detection_survival_vec, detection_indices_vec = calculate_detection_indices_and_survival(
        ews_optimal_simulation_df.detection_index,
        unique_detection_indices,
        nsims,
    )

    null_survival_vec, null_indices_vec = calculate_detection_indices_and_survival(
        ews_optimal_simulation_df.null_detection_index,
        unique_null_detection_indices,
        nsims,
    )

    return (
        (; detection_survival_vec, detection_indices_vec),
        (; null_survival_vec, null_indices_vec),
    )
end
