export line_plot

function line_plot(
        results::StructVector{OptimizationResult};
        xlabel = "Test Sensitivity & Specificity",
        ylabel = "Alert Accuracy",
        facet_fontsize = 20,
        legendsize = 22,
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 22,
        yticklabelsize = 22,
        legend_rowsize = Relative(0.05),
        xlabel_rowsize = Relative(0.03),
        ylabel_rowsize = Relative(0.02),
        kwargs...,
    )
    kwargs_dict = Dict{Symbol, Any}(kwargs)

    if !haskey(kwargs_dict, :ylims)
        ylims = (
            maximum([0, minimum(results.accuracy) - 0.05]),
            minimum([1, maximum(results.accuracy) + 0.05]),
        )
    else
        ylims = kwargs_dict[:ylims]
    end

    unique_noise_combos = unique(
        collect(zip(results.noise_level, results.noise_type_description))
    )
    noise_descriptions = map(unique_noise_combos) do (level, type_desc)
        if type_desc == :static
            noise_scaling = if level == 7.0
                "High"
            elseif level == 1.0
                "Low"
            else
                "Uncharacterized"
            end
            "$(noise_scaling) Static Noise"
        else
            noise_scaling = if level ≈ 0.102
                "High"
            elseif level ≈ 0.8734
                "Low"
            else
                "Uncharacterized"
            end
            "$(noise_scaling) Dynamical Noise"
        end
    end
    num_noise_descriptions = length(noise_descriptions)

    if !haskey(kwargs_dict, :nbanks)
        nbanks = num_noise_descriptions
    else
        nbanks = kwargs_dict[:nbanks]
    end

    unique_metrics = unique(results.ews_metric)
    num_metrics = length(unique_metrics)
    if num_metrics > 4
        error(
            "Trying to plot too many metric facets $(num_metrics). Max allowed is 4"
        )
    end

    fig = Figure()

    for (metric_num, metric) in enumerate(unique_metrics)
        ax_position = if metric_num == 1
            (1, 1)
        elseif metric_num == 2
            (1, 2)
        elseif metric_num == 3
            (2, 1)
        elseif metric_num == 4
            (2, 2)
        else
            @error "Too many metric types"
        end

        metric_mask = results.ews_metric .== metric
        metric_results = results[metric_mask]

        gl = fig[ax_position...] = GridLayout()
        line_plot_facet!(
            gl,
            metric_results,
            unique_noise_combos;
            ax_position = ax_position,
            num_metrics = num_metrics,
            facet_title = sentencecase(replace(metric, "_" => " ")),
            ylims = ylims,
            facet_fontsize = facet_fontsize,
            xticklabelsize = xticklabelsize,
            yticklabelsize = yticklabelsize,
        )
    end

    Legend(
        fig[0, :],
        [
            PolyElement(; color = col) for
                col in Makie.wong_colors()[1:num_noise_descriptions]
        ],
        noise_descriptions,
        "";
        nbanks = nbanks,
        labelsize = legendsize,
    )

    rowsize!(fig.layout, 0, legend_rowsize)

    Label(
        fig[:, 0],
        ylabel;
        fontsize = ylabelsize,
        font = :bold,
        rotation = pi / 2,
    )
    colsize!(fig.layout, 0, ylabel_rowsize)

    xlabel_position = num_metrics > 2 ? 3 : 2
    Label(
        fig[xlabel_position, :],
        xlabel;
        fontsize = xlabelsize,
        font = :bold,
    )
    rowsize!(fig.layout, xlabel_position, xlabel_rowsize)

    return fig
end

function line_plot_facet!(
        gl,
        results::StructVector{OptimizationResult},
        unique_noise_combos;
        ax_position = (1, 1),
        num_metrics = 1,
        facet_title = "",
        ylims = nothing,
        facet_fontsize = 20,
        legendsize = 22,
        xlabel = "",
        ylabel = "",
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 22,
        yticklabelsize = 22,
        kwargs...,
    )
    unique_tests = unique(results.test_specification)
    test_labels = map(
        test -> "$(Int64(round(test.sensitivity * 100; digits = 0)))%",
        unique_tests,
    )

    ax = Axis(
        gl[2, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        xticks = (collect(1:length(test_labels)), test_labels),
        titlesize = facet_fontsize,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
    )

    Box(gl[1, 1]; color = :lightgray, strokevisible = false)
    Label(
        gl[1, 1],
        facet_title;
        fontsize = facet_fontsize,
        padding = (0, 0, 0, 0),
        valign = :bottom,
        tellwidth = false,
    )
    rowsize!(gl, 2, Relative(0.9))

    for (noise_level, noise_type) in unique_noise_combos
        noise_mask = (results.noise_level .== noise_level) .&
            (results.noise_type_description .== noise_type)
        noise_results = results[noise_mask]

        if length(noise_results) > 0
            lines!(ax, 1:length(noise_results), noise_results.accuracy)
        end
    end

    if !(isnothing(ylims))
        ylims!(ax, ylims)
    end

    if ax_position[2] != 1
        hideydecorations!(ax)
    end

    nrows = ceil(num_metrics / 2)
    if ax_position[1] != nrows && ax_position[2] + 2 <= num_metrics
        hidexdecorations!(ax)
    end
    return nothing
end
