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
# Define test specifications to include in the plot
test_specs = [
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(0.95, 0.95, 0),
    IndividualTestSpecification(0.96, 0.96, 0),
    IndividualTestSpecification(0.97, 0.97, 0),
    IndividualTestSpecification(0.98, 0.98, 0),
    IndividualTestSpecification(0.99, 0.99, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

#%%
autocovariance_lineplot_data = prepare_line_plot_data(
    ensemble_results,
    "autocovariance",
    test_specs;
    tiebreaker_preference = "specificity"
)

#%%
# Create line plots for each EWS metric
plots_dir = plotsdir("accuracy-line-plots")
mkpath(plots_dir)

for metric in unique_metrics
    println("Creating line plot for metric: $metric")

    # Prepare data for this metric
    prepared_data = prepare_line_plot_data(
        ensemble_results,
        metric,
        test_specs;
        tiebreaker_preference = "specificity"
    )

    if length(prepared_data) == 0
        @warn "No data available for metric: $metric"
        continue
    end

    println("  Prepared $(length(prepared_data)) data points")

    # Create the plot
    local fig = line_plot(
        prepared_data;
        xlabel = "Test Sensitivity & Specificity",
        ylabel = "Alert Accuracy",
        nbanks = 2,
        facet_fontsize = facet_fontsize,
        legendsize = legendsize,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
    )

    # Save the plot
    local filename = "$(ensemble_name)_$(metric)_accuracy_line_plot"
    for ext in [".png", ".svg"]
        filepath = joinpath(plots_dir, filename * ext)
        save(filepath, fig)
        println("  Saved: $filepath")
    end
end

println("All line plots created successfully!")

#%%
# Optional: Create a combined plot with multiple metrics
# Select up to 4 metrics to plot together
selected_metrics = unique_metrics[1:min(4, length(unique_metrics))]

println("Creating combined plot with metrics: ", selected_metrics)

# Filter results to only include selected metrics
combined_mask = [r.ews_metric in selected_metrics for r in ensemble_results]
combined_results = ensemble_results[combined_mask]

# Prepare data for all selected metrics
combined_prepared = StructVector(OptimizationResult[])
for metric in selected_metrics
    metric_data = prepare_line_plot_data(
        combined_results,
        metric,
        test_specs;
        tiebreaker_preference = "specificity"
    )
    combined_prepared = vcat(combined_prepared, metric_data)
end

# Create combined plot
combined_fig = line_plot(
    combined_prepared;
    xlabel = "Test Sensitivity & Specificity",
    ylabel = "Alert Accuracy",
    nbanks = 2,
    facet_fontsize = facet_fontsize,
    legendsize = legendsize,
    xlabelsize = xlabelsize,
    ylabelsize = ylabelsize,
    xticklabelsize = xticklabelsize,
    yticklabelsize = yticklabelsize,
)

# Save combined plot
combined_filename = "$(ensemble_name)_combined_accuracy_line_plot"
for ext in [".png", ".svg"]
    filepath = joinpath(plots_dir, combined_filename * ext)
    save(filepath, combined_fig)
    println("Saved combined plot: $filepath")
end

println("Done!")
