using CairoMakie
# using GLMakie

function theme_adjustments()
    return Theme(;
        fontsize = 16,
        Axis = (;
            xlabelsize = 28,
            ylabelsize = 28,
            xticklabelsize = 24,
            yticklabelsize = 24,
            xlabelfont = :bold,
            ylabelfont = :bold,
        ),
        Colorbar = (;
            labelsize = 20,
            labelfont = :bold,
        ),
    )
end
custom_theme = merge(theme_adjustments(), theme_minimal())

set_theme!(
    custom_theme;
    fontsize = 16,
    linewidth = 4,
)

update_theme!(; size = (1300, 800))
CairoMakie.activate!(; type = "svg", pt_per_unit = 1.5)
# GLMakie.activate!()

facet_fontsize = 24
legendsize = 24
legendticklabelsize = 22
legendwidth = 20
xlabelsize = 28
ylabelsize = 28
xticklabelsize = 24
yticklabelsize = 24
