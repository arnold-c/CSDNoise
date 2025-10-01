export calculate_min_max_vaccination_range

function calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range = 0.2,
    )
    @assert mean_vaccination_coverage <= 1.0

    min_vaccination_range = minimum(
        [
            max_vaccination_range,
            1.0 - mean_vaccination_coverage,
            mean_vaccination_coverage,
        ]
    )

    min_vaccination_coverage = round(
        mean_vaccination_coverage - min_vaccination_range; digits = 4
    )
    max_vaccination_coverage = round(
        mean_vaccination_coverage + min_vaccination_range; digits = 4
    )
    return min_vaccination_coverage, max_vaccination_coverage
end
