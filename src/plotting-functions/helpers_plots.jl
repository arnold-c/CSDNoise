using GLMakie

export BASE_COLOR,
    OUTBREAK_COLOR,
    REFF_GT_ONE_COLOR

const BASE_COLOR = "#5E5C6C"
const OUTBREAK_COLOR = "#F0780F"
const REFF_GT_ONE_COLOR = "#00857E"

function clean_ews_metric_names(metric_vector)
    return replace.(
        titlecase.(
            replace.(
                metric_vector,
                "_" => " ",
            )
        ), "Of" => "of"
    )
end

function test_axis_label(test)
    val = Int64(round(test.sensitivity * 100; digits = 0))
    return "$val%"
end

function line_and_hline!(
        ax,
        times,
        yvec,
        hlinevec;
        linewidth = 3,
        hlinecolor = :black,
        hlinewidth = 2,
        hlinestyle = :dash,
        kwargs...,
    )
    kwargs_dict = Dict(kwargs)

    color_keys = [:color, :colormap]

    if !(
            sum(k -> haskey(kwargs_dict, k), color_keys) in (0, length(color_keys))
        )
        error(
            "Must provide both $color_keys for the `lines!` function, or none. You have specified: $(keys(kwargs_dict))"
        )
    end

    lines!(
        ax,
        times,
        yvec;
        linewidth = linewidth,
        kwargs...,
    )
    hlines!(
        ax,
        hlinevec;
        color = hlinecolor,
        linewidth = hlinewidth,
        linestyle = hlinestyle,
    )
    return nothing
end
