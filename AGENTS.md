# AGENTS.md - Development Guidelines for CSDNoise

## Build/Test Commands
- **Run all tests**: `julia --project=. -e "using Pkg; Pkg.test()"`
- **Run single test file**: `julia --startup-file=no -e 'using Test; @testset "filename" include("test/filename.jl")'`
- **Run optimizations**: `just optimizations` (runs gridsearch with auto-threading)
- **Build manuscript**: `just manuscript` (requires Just and Typst)
- **Format code**: `runic -i path/to/file` (ALWAYS format after editing/creating files)

## Code Style Guidelines
- **Language**: Julia with functional programming style preferred
- **Imports**: Use explicit imports at module level, group by package (see CSDNoiseCore/src/CSDNoiseCore.jl:1-35)
- **Types**: Use `Base.@kwdef` for structs with defaults, `LightSumTypes.@sumtype` for ADTs, multiple dispatch over nested if statements
- **Naming**: snake_case for functions/variables/files, PascalCase for types/modules
- **Error handling**: Use Try.jl and TryExperimental.jl for operations with potential failures
- **Performance**: Follow Julia performance tips, ensure type stability with JET.jl, use `@unstable` from DispatchDoctor for unstable functions
- **Testing**: Include Aqua.jl and JET.jl tests for code quality and type stability (see CSDNoiseCore/test/runtests.jl)
- **Documentation**: Use docstrings for exported functions with Args/Returns/Throws/Examples sections (see CSDNoiseCore/src/ews-metrics/ews-metrics.jl:3-30)
- **Exports**: Export at top of each file, not in module definition

## Project Structure
- `CSDNoiseCore/`: Core package with all simulation/optimization logic
- `src/`: Top-level plotting and analysis functions
- `test/`: Test files (mirror src/ structure)
- `scripts/`: Analysis scripts using DrWatson @quickactivate
- `manuscript/`: Typst manuscript files and generation scripts

## Key Dependencies
DrWatson, LightSumTypes, Try/TryExperimental, DispatchDoctor, JET, Aqua, MultistartOptimization, OhMyThreads, StructArrays, Distributions