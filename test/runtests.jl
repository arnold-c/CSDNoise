using CSDNoise
using Test
using Aqua
using JET

@testset "CSDNoise.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CSDNoise; ambiguities = false)
        @testset "Ambiguities" begin
            Aqua.test_ambiguities(CSDNoise)
        end
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(CSDNoise; target_defined_modules = true)
    end
    include("ews-functions.jl")
    include("ews-metrics.jl")
    include("ews-multistart-optimization.jl")
end
