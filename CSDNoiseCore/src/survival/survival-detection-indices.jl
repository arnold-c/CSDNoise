export calculate_detection_indices_and_survival

function calculate_detection_indices_and_survival(
        detection_vec,
        detection_indices,
        nsims,
    )
    if isempty(detection_indices)
        @assert nsims ==
            sum(detection_vec .== nothing)

        survival_vec = [nsims, nsims]
        indices_vec = Int64[]
    else
        detection_indices_counts = map(
            ((i, v),) -> sum(
                detection_vec .== v
            ),
            enumerate(detection_indices),
        )

        @assert nsims ==
            sum(detection_indices_counts) +
            sum(detection_vec .== nothing)

        survival_vec = nsims .- cumsum(detection_indices_counts)

        survival_vec = vcat(
            nsims, nsims, survival_vec..., survival_vec[end]
        )
        detection_indices = Int64.(detection_indices)

        indices_vec = vcat(
            detection_indices[1],
            detection_indices,
        )
    end

    return survival_vec, indices_vec
end
