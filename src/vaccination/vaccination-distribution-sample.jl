using Distributions: Distributions

export sample_vaccination_coverage

function sample_vaccination_coverage(
        min_coverage,
        max_coverage,
        digits = 4
    )
    return round(
        rand(
            Distributions.Uniform(
                min_coverage,
                max_coverage
            )
        );
        digits = digits
    )
end
