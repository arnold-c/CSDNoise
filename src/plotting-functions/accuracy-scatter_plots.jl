export accuracy_scatter_plot

function accuracy_scatter_plot(
        results::StructVector{OptimizationResult},
        metric::String;
        xlabel = "Test Sensitivity",
        ylabel = "Alert Accuracy",
        title = nothing,
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 22,
        yticklabelsize = 22,
        legendsize = 22,
        markersize = 10,
        jitter_amount = 0.005,
        kwargs...,
    )
    metric_mask = results.ews_metric .== metric
    filtered_results = results[metric_mask]

    if length(filtered_results) == 0
        error("No results found for metric $(metric)")
    end

    unique_tests = unique(filtered_results.test_specification)
    test_labels = map(
        test -> "$(Int64(round(test.sensitivity * 100; digits = 0)))%",
        unique_tests,
    )
    test_positions = [1.0 - t.sensitivity for t in unique_tests]

    fig = Figure()

    plot_title = if isnothing(title)
        "$(metric)"
    else
        title
    end

    ax = Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = plot_title,
        xticks = (test_positions, test_labels),
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
    )

    groups = group_structvector(
        filtered_results,
        :noise_level,
        :noise_type_description
    )

    num_groups = length(groups)
    jitter_offsets = range(-jitter_amount, jitter_amount, length = num_groups)

    for (idx, (key, group)) in enumerate(groups)
        noise_level = key.noise_level
        noise_type = key.noise_type_description

        noise_scaling = if noise_level == 7.0
            "High"
        elseif noise_level == 1.0
            "Low"
        else
            "Uncharacterized"
        end

        noise_description = if noise_type == :static
            "$(noise_scaling) Static Noise"
        else
            "$(noise_scaling) Dynamical Noise"
        end

        test_sensitivities = [1.0 - t.sensitivity for t in group.test_specification]
        jittered_sensitivities = test_sensitivities .+ jitter_offsets[idx]

        scatter!(
            ax,
            jittered_sensitivities,
            group.accuracy;
            color = Makie.wong_colors()[idx],
            label = noise_description,
            markersize = markersize,
        )
    end

    Legend(
        fig[1, 2],
        ax;
        labelsize = legendsize,
    )

    return fig
end
