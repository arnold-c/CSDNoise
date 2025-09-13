# AGENTS.md - Development Guidelines for CSDNoise

## Build/Test Commands
- **Run all tests**: `julia --project=. -e "using Pkg; Pkg.test()"`
- **Run single test file**: `julia --startup-file=no -e 'using Test; @testset "filename" include("test/filename.jl")'`
- **Build manuscript**: `just manuscript` (requires Just and Typst)
- **Run tests + manuscript**: `just` (default target)
- **Format code**: `runic -i path/to/file` (use Runic.jl for formatting)

## Code Style Guidelines
- **Language**: Julia with functional programming style preferred
- **Imports**: Use explicit imports, group by package (see src/ews-functions.jl:1-10)
- **Types**: Use multiple dispatch over nested if statements, leverage SumTypes for ADTs
- **Naming**: snake_case for functions/variables, PascalCase for types/modules
- **Error handling**: Use Try.jl package for operations with potential failures
- **Performance**: Follow Julia performance tips, ensure type stability with JET.jl
- **Testing**: Include Aqua.jl and JET.jl tests for code quality and type stability
- **Documentation**: Use docstrings for exported functions, follow existing patterns

## Project Structure
- `src/`: Core Julia modules and functions
- `test/`: Test files (mirror src/ structure)
- `scripts/`: Analysis and exploration scripts
- `manuscript/`: Typst manuscript files and generation scripts

## Dependencies
Uses DrWatson.jl for project management. Key packages: DataFrames, StatsBase, Distributions, CairoMakie, JET, Try, SumTypes