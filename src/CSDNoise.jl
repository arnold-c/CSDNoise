module CSDNoise
using CSDNoiseCore
using GLMakie

# Plotting Functions
include("./plotting-functions/plotting-setup_makie.jl")
include("./plotting-functions/helpers_plots.jl")
include("./plotting-functions/tau-heatmap_plots.jl")
include("./plotting-functions/lead-time_plots.jl")
include("./plotting-functions/single-simulation_plots.jl")
include("./plotting-functions/simulation-ews_plots.jl")
include("./plotting-functions/single-scenario_plots.jl")
include("./plotting-functions/simulation-optimal-ews_plots.jl")
include("./plotting-functions/hyperparam-debugging_plots.jl")
include("./plotting-functions/accuracy-lines_data-preparation.jl")
include("./plotting-functions/accuracy-lines_plots.jl")
include("./plotting-functions/ews-heatmap_plots.jl")
include("./plotting-functions/survival/ews-survival_facet-parameter-preparation.jl")
include("./plotting-functions/survival/ews-survival_facet.jl")
include("./plotting-functions/survival/ews-survival_plot-wrapper.jl")
include("./plotting-functions/survival/ews-survival_Reff-histogram.jl")
include("./plotting-functions/survival/ews-survival_simulate-and-plot.jl")
include("./plotting-functions/survival/ews-survival_survival-lines.jl")
include("./plotting-functions/survival/Reff_histogram_plot.jl")

end
