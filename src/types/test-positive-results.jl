export EnsembleTestResultRun

Base.@kwdef struct EnsembleTestResultRun
    emergent_test_positives::Vector{Vector{Int64}}
    null_test_positives::Vector{Vector{Int64}}
end
