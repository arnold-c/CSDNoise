export SEIRRun,
    NoiseRun

Base.@kwdef struct SEIRRun
    states::Vector{StaticArrays.SVector{5, Int64}}
    incidence::Vector{Int64}
    Reff::Vector{Float64}
end

Base.@kwdef struct NoiseRun
    incidence::Vector{Vector{Int64}}
    mean_noise::Float64
    mean_poisson_noise::Float64
    mean_dynamic_noise::Float64
end
