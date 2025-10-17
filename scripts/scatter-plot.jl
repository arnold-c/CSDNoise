#%%
using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using CSDNoise: CSDNoiseCore
using CSDNoiseCore: IndividualTestSpecification, load_previous_gridsearch_results_structvector, OptimizationResult
using GLMakie
using StructArrays: StructVector

# Setup plotting theme
include(scriptsdir("plotting-setup.jl"));

#%%
# Load gridsearch optimization results
results_dir = projectdir("out", "ensemble", "ews-hyperparam-gridsearch")
filename_base = "ews-hyperparam-gridsearch-structvector"

println("Loading gridsearch results from: $results_dir")
all_results = load_previous_gridsearch_results_structvector(results_dir, filename_base);
println("Loaded $(length(all_results)) optimization results")

#%%
# Get unique ensemble names
unique_ensembles = unique([r.ensemble_specification.label for r in all_results])
println("Available ensembles: ", unique_ensembles)

# Filter for a specific ensemble (e.g., measles)
# Change this to the ensemble you want to plot
ensemble_name = "measles"
ensemble_mask = [r.ensemble_specification.label == ensemble_name for r in all_results]
ensemble_results = all_results[ensemble_mask]

println("Filtered to $(length(ensemble_results)) results for ensemble: $ensemble_name")

#%%
# Get unique EWS metrics
unique_metrics = unique(ensemble_results.ews_metric)
println("Available EWS metrics: ", unique_metrics)

#%%
# Select metric to plot
metric = "autocovariance"

#%%
# Create scatter plot
fig = accuracy_scatter_plot(
    ensemble_results,
    metric;
    title = "$(ensemble_name) - $(metric)",
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    xticklabelsize = xticklabelsize,
    yticklabelsize = yticklabelsize,
    legendsize = legendsize,
);

# Save the plot
plots_dir = plotsdir("accuracy-scatter-plots")
mkpath(plots_dir)

filename = "$(ensemble_name)_$(metric)_scatter"
for ext in [".png", ".svg"]
    filepath = joinpath(plots_dir, filename * ext)
    save(filepath, fig)
    println("Saved: $filepath")
end

println("Scatter plot created successfully!")
