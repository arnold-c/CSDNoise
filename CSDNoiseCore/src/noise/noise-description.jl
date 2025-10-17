export get_noise_description,
    get_noise_magnitude,
    get_noise_magnitude_description,
    getdirpath

get_noise_description(noise_specification::NoiseSpecification) = get_noise_description(LightSumTypes.variant(noise_specification))

function get_noise_description(noise_specification::PoissonNoise)
    return "poisson"
end

function get_noise_description(noise_specification::DynamicalNoise)
    return string("dynamical, ", noise_specification.correlation)
end

get_noise_magnitude_description(noise_specification::NoiseSpecification) = get_noise_magnitude_description(LightSumTypes.variant(noise_specification))

function get_noise_magnitude_description(noise_specification::Union{PoissonNoise, DynamicalNoise})
    return string("Poisson scaling: ", noise_specification.noise_mean_scaling)
end

get_noise_magnitude(noise_specification::NoiseSpecification) = get_noise_magnitude(LightSumTypes.variant(noise_specification))

function get_noise_magnitude(noise_specification::Union{PoissonNoise, DynamicalNoise})
    return noise_specification.noise_mean_scaling
end

noise_table_description(noise_specification::NoiseSpecification) = noise_table_description(LightSumTypes.variant(noise_specification))

function noise_table_description(noise_specification::PoissonNoise)
    noise_scaling = if noise_specification.noise_mean_scaling == 7
        "High"
    elseif noise_specification.noise_mean_scaling == 1
        "Low"
    else
        "Uncharacterized"
    end
    return "$(noise_scaling) Static Noise"
end

function noise_table_description(noise_specification::DynamicalNoise)
    avg_vaccination = round(
        mean(
            [
                noise_specification.min_vaccination_coverage,
                noise_specification.max_vaccination_coverage,
            ]
        );
        digits = 4,
    )
    noise_scaling = if avg_vaccination == 0.102
        "High"
    elseif avg_vaccination == 0.8734
        "Low"
    else
        "Uncharacterized"
    end
    return "$(noise_scaling) Dynamical Noise"
end

function getdirpath(spec::NoiseSpecification)
    return getdirpath(LightSumTypes.variant(spec))
end

function getdirpath(spec::Union{PoissonNoise, DynamicalNoise})
    return reduce(
        joinpath,
        map(
            p -> "$(p)_$(getproperty(spec, p))",
            propertynames(spec),
        ),
    )
end
