using StaticArrays: SVector

export SEIRRun,
    NoiseRun

struct SEIRRun
    states::Vector{SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

struct NoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
    mean_poisson_noise::Float64
    mean_dynamic_noise::Float64
end
