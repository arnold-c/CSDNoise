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

To run the project, you can use the Makefile.
First, you should make sure you have Make installed on your computer.
If you have macOS, it should be pre-installed, but if you'd like to use the most up-to-date version, you can  install it using homebrew.

```bash
brew install make
```

> If you use the homebrew version of Make then you will need to use `gmake` instead of `make` to run the Makefile.

If you are on Linux, you can install Make with your package manager.
If you are on Windows, you can use the [Chocolatey](https://chocolatey.org/) package manager, download it [directly](https://gnuwin32.sourceforge.net/packages/make.html), or use the ubuntu package manager on [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) (recommended as it will give you access to a linux terminal).

Once you have a version of Make installed you can use the following terminal command to run the Makefile:

```bash
make
```

If you would like to run a specific target, you can specify it after the `make` command e.g. `make ensemble-sim`

Make will track when each file was last run and save a temporary file to the `tmp/` directory.
If you want to force a re-run then you can delete the associated temporary file, either manually, or by running the associated Make clean command e.g. `make clean-ensemble-sims` to delete all ensemble simulation files.

## Project Structure

```bash
.
├── Makefile
├── Manifest.toml
├── Project.toml
├── README.md
├── _research
├── data
│   ├── CFR_2022.csv
│   ├── input-populations.csv
│   ├── optimal-threshold-results
│   │   ├── R0_12.0
│   │   │   ├── noise_type_dynamical
│   │   │   └── noise_type_poisson
│   │   ├── R0_16.0
│   │   │   ├── noise_type_dynamical
│   │   │   └── noise_type_poisson
│   │   ├── R0_20.0
│   │   │   ├── noise_type_dynamical
│   │   │   └── noise_type_poisson
│   │   └── R0_8.0
│   │       ├── noise_type_dynamical
│   │       └── noise_type_poisson
│   ├── seasonal-infectivity-import
│   │   └── tau-leaping
│   └── singlesim
├── notebooks
│   ├── Julia
│   └── R
├── plots
│   ├── ensemble
│   │   ├── optimal-thresholds
│   │   │   ├── R0_12.0
│   │   │   │   ├── noise_type_dynamical
│   │   │   │   └── noise_type_poisson
│   │   │   │       └── ...
│   │   │   │           ├── clinic-tested
│   │   │   │           └── tests
│   │   │   ├── R0_16.0
│   │   │   │   ├── noise_type_dynamical
│   │   │   │   └── noise_type_poisson
│   │   │   │       └── ...
│   │   │   │           ├── clinic-tested
│   │   │   │           └── tests
│   │   │   ├── R0_20.0
│   │   │   │   ├── noise_type_dynamical
│   │   │   │   └── noise_type_poisson
│   │   │   │       └── ...
│   │   │   │           ├── clinic-tested
│   │   │   │           └── tests
│   │   │   └── R0_8.0
│   │   │   │   ├── noise_type_dynamical
│   │   │   │   └── noise_type_poisson
│   │   │   │       └── ...
│   │   │   │           ├── clinic-tested
│   │   │   │           └── tests
│   │   ├── single-scenario
│   │   └── testing-comparison
│   │       ├── noise_type_dynamical
│   │       └── noise_type_poisson
│   └── singlesim
├── renv.lock
├── scripts
│   ├── debugging.jl
│   ├── ensemble-diag-testing_optimal-thresholds.jl
│   ├── ensemble-diag-testing_scenarios_plots.jl
│   ├── ensemble-sim.jl
│   ├── ensemble-sim_noise-visualizations.jl
│   ├── ensemble-sim_single-scenario.jl
│   ├── single-sim.jl
│   └── single-sim_plots.jl
├── src
│   ├── DrWatson-helpers.jl
│   ├── OutbreakDetection.jl
│   ├── R
│   │   └── app.R
│   ├── SEIR-model.jl
│   ├── cleaning-functions.jl
│   ├── detection-thresholds.jl
│   ├── diag-testing-functions.jl
│   ├── dynamics-constants.jl
│   ├── ensemble-functions.jl
│   ├── ensemble-parameters.jl
│   ├── ensemble-sim_single-scenario_plots.jl
│   ├── makie-plotting-setup.jl
│   ├── noise-functions.jl
│   ├── optimal-threshold-functions.jl
│   ├── plotting-functions.jl
│   ├── single-sim_setup.jl
│   ├── structs.jl
│   ├── test-constants.jl
│   ├── threshold_comparison_plots.jl
│   └── transmission-functions.jl
├── test
│   └── runtests.jl
├── tmp
└── workflows
    └── CI.yml
```

- `_research`
- `data/` contains input and output data files
    - `CFR_2022.csv` contains CFR rates for representative countries in 2022
    - `input-populations.csv` contains population sizes for representative
    - `optimal-threshold-results` contains output excel tables of the optimal threshold results, separated into subdirectories by `R_0` of the simulation and the noise type
    - `seasonal-infectivity-import` contains the output data files of the outbreak detection characteristics for the ensemble simulations, separated into subdirectories by model specification. Files are saved in the Julia's HDF5-compliant `.jld2` format
    - `singlesim` contains data file for a single simulation (setup files and the output arrays)
- `notebooks` contains short notebooks to perform temporary analyses using Quarto and Rmarkdown documents
- `plots` contains all output plots
    - `ensemble` contains all plots related to the ensemble simulation
        - `optimal-thresholds` contains plots related to the optimal alert thresholds for each simulation type, separated by simulation `R_0` and noise type
            - `clinic-tested` contains the optimal threshold plots where each plot refers to a different level of the % of clinic visits that are tested, and the rows refer to the test type
            - `tests` contains the optimal threshold plots where each plot refers to a different test type and the rows refer to a different % of clinic visits that are tested
        - `single-scenario` contains the plots for a single scenario of the ensemble simulations, with subdirectories for noise type where appropriate (i.e., for alert-related metrics)
        - `testing-comparison` contains plots for alert metrics compared across test type and testing rate, separated into subdirectories by noise type. These figures are computationally expensive to compute so only produced for `R_0` = 16
    - `singlesim` contains plots for the single simulation
- `renv` contains the R package versions used for the prototype app examining the trade-off between test sensitivity and specificity and the detections ability when a number of true positives are tested
- `scripts` contains the Julia scripts used to examine single and ensemble simulations, using plotting and other functions defined in `src/*.jl` files
- `src` contains all Julia source files and functions used in the analysis pipeline and exploration scripts. These files are separated by purpose e.g., `cleaning-functions.jl` contains functions for cleaning the simulation arrays into dataframes for simpler plotting and manipulation, and `ensemble-functions.jl` contains all functions related to running the ensemble simulations. There is also an `R/` directory that runs the prototype R shiny app
- `test` contains all test scripts
- `tmp` contains all temporary files created to track dependencies for Make. This will likely be removed when the shifting to use [Just](https://github.com/casey/just) for the pipeline
- `workflows` contains the CI workflow using GitHub Actions. Currently it only contains a file that can run tests on push to the `main` branch, but it is not active
