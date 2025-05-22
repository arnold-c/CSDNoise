# CSDNoise

This code base is using the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named
> CSDNoise

## Getting Started
### Installing dependencies

To (locally) reproduce this project, do the following:

0. Download this code base.
You can do this by cloning the repository or by downloading the zip file.
Notice that raw data are not included in the git-history and may need to be downloaded independently.
1. Install Julia. I would recommend you use the [`juliaup` installer](https://github.com/JuliaLang/juliaup) as it makes it much easier to deal with multiple versions of Julia, as well as keep them up to date.
2. Open a Julia console and do:

   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> using DrWatson
   julia> @quickactivate "CSDNoise"
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:

```julia
using DrWatson
@quickactivate "CSDNoise"
```
which auto-activate the project and enable local path handling from DrWatson.

### Running the project

To run the project, you can use the Justfile.
First, you should make sure you have both Just and Typst installed on your computer.
You can do so by following the instructions at [the Just manual](https://just.systems/man/en/introduction.html) and [Typst's GitHub Page for the CLI](https://github.com/typst/typst).

Once you have a version of Just and Typst installed you can use the following terminal command to run the Justfile:

```bash
just manuscript
```


## Project Structure

```bash
❯ tree -L 2
.
├── Justfile
├── Manifest.toml
├── manuscript
│   ├── combined-manuscript.pdf
│   ├── combined-manuscript.typ
│   ├── CSD.bib
│   ├── manuscript_files
│   ├── manuscript.pdf
│   ├── manuscript.typ
│   ├── scripts
│   │   ├── optimal-thresholds.jl
│   │   └── plotting-setup.jl
│   ├── supplemental_files
│   ├── supplemental-appendix.pdf
│   ├── supplemental-appendix.typ
│   └── template.typ
├── out
│   └── ensemble
│       ├── ews-hyperparam-optimization
│       └── seasonal-infectivity-import
├── plots
├── Project.toml
├── README.md
├── src
│   ├── ARCHIVE
│   ├── cairomakie-plotting-setup.jl
│   ├── CSDNoise.jl
│   ├── detection-thresholds.jl
│   ├── diag-testing-functions.jl
│   ├── dynamics-constants.jl
│   ├── ensemble-functions.jl
│   ├── ensemble-sim_single-scenario_plots.jl
│   ├── ews-functions.jl
│   ├── ews-hyperparam-optimization.jl
│   ├── ews-metrics.jl
│   ├── ews-survival.jl
│   ├── helpers.jl
│   ├── makie-plotting-setup.jl
│   ├── noise-functions.jl
│   ├── plotting-functions
│   ├── SEIR-model.jl
│   ├── structs.jl
│   ├── test-constants.jl
│   └── transmission-functions.jl
├── test
│   ├── ews-functions.jl
│   ├── ews-metrics.jl
│   └── runtests.jl
└── workflows
    └── CI.yml

25 directories, 86 files
```

- `plots` contains all output plots
- `scripts` contains the Julia scripts used to examine single and ensemble simulations, using plotting and other functions defined in `src/*.jl` files
- `src` contains all Julia source files and functions used in the analysis pipeline and exploration scripts. These files are separated by purpose.
- `manuscript/` contains all files relevant to the manuscript
    - `combined-manuscript.typ` is an outer file that includes both the `manuscript.typ` and `supplemental-appendix.typ` files for cross-referencing.
    - `scripts/optimal-thresholds.jl` contains all function calls to generate the underlying simulation results and figures and tables for the manuscript
- `test` contains all test scripts
- `workflows` contains the CI workflow using GitHub Actions. Currently it only contains a file that can run tests on push to the `main` branch, but it is not active
