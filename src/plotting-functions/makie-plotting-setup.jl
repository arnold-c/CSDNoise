using GLMakie
function theme_adjustments()
    return Theme(;
        fontsize = 16,
        Axis = (;
            xlabelsize = 20,
            ylabelsize = 20,
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

# Laptop monitor
# update_theme!(; size = (1220, 640))
update_theme!(; size = (2100, 1200)) # Home monitor
# update_theme!(; size = (950, 550)) # Work monitor
# CairoMakie.activate!(; type = "png", px_per_unit = 1.0)
GLMakie.activate!(; float = true)
