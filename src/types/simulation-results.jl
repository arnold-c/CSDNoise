export SEIRRun,
    EnsembleSEIRRun,
    NoiseRun


Base.@kwdef struct SEIRRun
    states::Vector{StaticArrays.SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

Base.@kwdef struct EnsembleSEIRRun
    emergent_seir_run::StructVector{SEIRRun}
    null_seir_run::StructVector{SEIRRun}
end

abstract type AbstractNoiseRun end

Base.@kwdef struct DynamicalNoiseRun
    incidence::Vector{Vector{Int64}}
    Reff::Vector{Vector{Float64}}
    mean_noise::Float64
    mean_poisson_noise::Float64
    mean_dynamic_noise::Float64
end

Base.@kwdef struct PoissonNoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
end

LightSumTypes.@sumtype NoiseRun(DynamicalNoiseRun, PoissonNoiseRun) <: AbstractNoiseRun
