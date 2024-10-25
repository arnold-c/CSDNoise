# %%
library(graphics) # standard R package
library(grDevices) # standard R package
library(zoo)
library(colorspace) # for constructing color
library(latex2exp) # used to include LaTeX expressions in plot text
library(scales)

## Additional packages for managing real data
library(lubridate)
library(padr)
library(imputeTS)


## Set UTC time
localTZ <- Sys.timezone()
Sys.setenv(TZ = "UTC")

# Load data
load(here::here("out", "tycho", "R-scripts", "tycho-params.RData"))

tychoL1.measles.CA.cases.imp.zoo.plotwindow <- readRDS(here::here("out", "tycho", "R-scripts", "tycho_CA_measles_plotwindow.rds"))
tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk <- readRDS(here::here("out", "tycho", "R-scripts", "tycho_CA_measles_plotwindow_2wk.rds"))
tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk <- readRDS(here::here("out", "tycho", "R-scripts", "tycho_CA_measles_plotwindow_4wk.rds"))

tychoL1.measles.CA.cases.imp.zoo.statswindow <- readRDS( here::here("out", "tycho", "R-scripts","tycho_CA_measles_statswindow_1wk.rds"))
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk <- readRDS( here::here("out", "tycho", "R-scripts","tycho_CA_measles_statswindow_2wk.rds"))
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk <- readRDS( here::here("out", "tycho", "R-scripts","tycho_CA_measles_statswindow_4wk.rds"))

tychoL1.measles.CA.cases.imp.zoo.statswindow.stats <- readRDS( here::here("out", "tycho", "R-scripts","tycho_CA_measles_stats_1wk.rds"))
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats <- readRDS( here::here("out", "tycho", "R-scripts","tycho_CA_measles_stats_2wk.rds"))
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats <- readRDS( here::here("out", "tycho", "R-scripts","tycho_CA_measles_stats_4wk.rds"))

# %%
# Color palettes ----------------------------------------------------------
unipalette.lorder <- c(
  "black" = "#252525", # Nero (grey)
  "purple" = "#5e2b7b", # Blue Diamond (violet)
  "red" = "#a11c3e", # Fire Brick (red)
  "blue" = "#226e83", # Allports (blue)
  "green" = "#319045", # Sea Green (green)
  "lightblue" = "#5798d1", # Picton Blue (blue)
  "pink" = "#e2908c" # Sea Pink (red)
)
unipalette <- unipalette.diff <- unipalette.lorder[c(3, 6, 1, 5, 2, 7, 4)]
unipalette.hybrid <- unipalette.lorder[c(1, 3, 2, 5, 4, 7, 6)]

AUC.colors <- diverge_hcl(
  n = 20, # number of colors
  h = c(45, 225), # Hues (low, hi)
  c = 100, # fixed Chroma or Chroma Range (edges, center)
  l = c(90, 10), # Lightness range (edges, center)
  power = 1 # exponent
)


# Set Y ranges ------------------------------------------------------------

# %%
cases.min <- 0
cases.max <- max(
  window(tychoL1.measles.CA.cases.imp.zoo.plotwindow, start = plotxmin, end = plotxmax),
  window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk, start = plotxmin, end = plotxmax),
  window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk, start = plotxmin, end = plotxmax)
)
cases.tick.interval <- 250
cases.ticks <- seq(from = cases.min, by = cases.tick.interval, length.out = ceiling(cases.max / cases.tick.interval) + 1)

var.min <- 0
var.max.actual <- max(
  tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$stats$variance,
  tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$stats$variance,
  tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$stats$variance
)
var.max <- 20000 # manually set plot limit

var.tick.interval <- 2000
var.ticks <- seq(from = var.min, to = var.max, by = var.tick.interval)

# Visual Variables --------------------------------------------------------

# %%
## Colors
color.1wk <- unipalette["black"]
color.2wk <- unipalette["lightblue"]
color.4wk <- unipalette["red"]
color.omitted <- rgb(0, 0, 0, .1)

## Line widths
lwd.1wk <- 0.5
lwd.2wk <- 1.0
lwd.4wk <- 2.0

## Filled line border widths and fill alpha
border.1wk <- 0.5
border.2wk <- 1.0
border.4wk <- 2.0
fill.alpha.1wk <- 0.7
fill.alpha.2wk <- 0.2
fill.alpha.4wk <- 0.0

## Typography
font.family <- "Times"
font.sizes <- seq(
  from = 8, # publisher's minimum point size (points)
  to = 12, # publisher's maximum point size (points)
  length.out = 5
)
font.size.normal <- mean(font.sizes)
font.scales <- font.sizes / mean(font.sizes)
names(font.scales) <- names(font.sizes) <- c("XS", "S", "M", "L", "XL")

## Figure dimensions
figure.widths <- c(min = 2.63, page = 7.5, column = 5.2) # in inches, as defined by publisher
figure.heights <- c(min = 1, page = 8.75) # in inches, as defined by publisher

## Margins and Figure Bounds
margins <- c(4, 5, 4, 8) + 0.1
top.pannel <- c(0, 1, .45, 1) # panel bounds: x0,x1,y0,y1 as fraction of figure region
bottom.pannel <- c(0, 1, 0, .6) # panel bounds: x0,x1,y0,y1 as fraction of figure region

# Init Figure -------------------------------------------------------------

# %%
## PDF output
plot_dir <- here::here("plots", "tycho", "R-plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
filename <- file.path(plot_dir, "fig1.pdf")

pdf(
  file = filename,
  title = "Figure 1", # displayed in title bar of PDF Reader
  width = figure.widths["page"], # inches.  Must fit publisher's min and max figure dimensions
  height = figure.heights["page"] * .7, # inches.  Must fit publisher's min and max figure dimensions # change to .7
  family = font.family,
  pointsize = font.size.normal # default size of text (points).
)

## init figure region
par(mar = margins)
par(lend = "butt", ljoin = "mitre")

# Top Panel (Variance) ----------------------------------------------------

## Initialize plot of stats
par(fig = top.pannel)
plot(0, 0,
  type = "n", axes = FALSE, ann = FALSE, yaxt = "n", xaxs = "i",
  xlim = c(plotxmin, plotxmax), ylim = range(var.ticks)
)

## Shade unanalyzed times
rect(obsdate, var.min, plotxmax + bandwidth, var.max, border = F, col = color.omitted)

## Grid
abline(h = var.ticks, col = "grey")

## weekly reports
var_1wk <- zoo(
  tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$stats$variance,
  time(tychoL1.measles.CA.cases.imp.zoo.statswindow)
)
lines(var_1wk, col = color.1wk, lwd = lwd.1wk)

## bi-weekly reports
var_2wk <- zoo(
  tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$stats$variance,
  time(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk)
)
lines(var_2wk, col = color.2wk, lwd = lwd.2wk)

## four-weekly reports
var_4wk <- zoo(
  tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$stats$variance,
  time(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk)
)
lines(var_4wk, col = color.4wk, lwd = lwd.4wk)

## Axes

axis(
  side = 1, # 1 specifies bottom axis - major ticks
  at = seq(plotxmin, plotxmax, "years"), # tick locations (in data units).
  labels = FALSE,
  tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
  line = 1, # axis position (in lines)
  lty = "solid",
  lwd = 0, # axis line weight
  lwd.ticks = 1, # tick line weight
  col = "grey" # axis line color
)
axis(
  side = 3, # 1 specifies bottom axis - major ticks
  at = seq(plotxmin, plotxmax, "years"), # tick locations (in data units).
  labels = year(seq(plotxmin, plotxmax, "years")),
  tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
  line = .5, # axis position (in lines)
  lty = "solid",
  lwd = 0, # axis line weight
  lwd.ticks = 1, # tick line weight
  col = "grey" # axis line color
)
axis(
  side = 2, # 1 specifies left axis
  at = var.ticks, # tick locations (in data units).
  labels = var.ticks, # vector of tick labels
  las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
  tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
  line = .5, # axis position, in lines
  lty = "solid",
  lwd = 0, # axis line weight
  lwd.ticks = 1, # tick line weight
  col = "grey" # axis line color
)

## Axis Labels
mtext(text = "Year", side = 3, line = 2.5)
title(ylab = "Variance", line = 4)

## Taus (Tau depends on all plotted values of the variance)
text(
  x = c(statsend - 108 * 7, statsend + 6 * 7, statsend + 6 * 7),
  y = c(9500, 9500, 4500),
  adj = 0,
  cex = font.scales["M"],
  col = c(
    color.4wk,
    color.2wk,
    color.1wk
  ),
  labels = c(
    TeX(sprintf("$\\tau = $ %.2f", tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats$taus$variance)),
    TeX(sprintf("$\\tau = $ %.2f", tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats$taus$variance)),
    TeX(sprintf("$\\tau = $ %.2f", tychoL1.measles.CA.cases.imp.zoo.statswindow.stats$taus$variance))
  )
)

## Legend
legend("topleft",
  xpd = NA, inset = c(1.01, 0), xjust = 0, yjust = 0, cex = font.scales["XS"], y.intersp = 3,
  legend = c("Four-weekly\nReports", "Bi-weekly\nReports", "Weekly\nReports"),
  col = c(color.4wk, color.2wk, color.1wk),
  lwd = c(lwd.4wk, lwd.2wk, lwd.1wk),
  lty = 1,
  bty = "n"
)

# Bottom Panel ------------------------------------------------------------

## Initialize plot of cases
par(fig = bottom.pannel, new = TRUE)
plot(0, 0,
  type = "n", axes = FALSE, ann = FALSE, yaxt = "n", xaxs = "i", yaxs = "i",
  xlim = c(plotxmin, plotxmax), ylim = range(cases.ticks)
)

## shade unanalyzed times
rect(obsdate, range(cases.ticks)[1], plotxmax + bandwidth, range(cases.ticks)[2], border = F, col = color.omitted)

## Grid
abline(h = cases.ticks, col = "grey")

## Plot time series of cases, aggreaged weekly, bi-weekly, four-weekly

filledLine <- function(data, lwd, color, fill.alpha) {
  x <- time(data)
  y <- as.numeric(data)
  data.x <- c(head(x, 1), x, tail(x, 1))
  data.x <- c(rbind(x, x)) ## double each time
  data.y <- c(
    0,
    c(rbind(y[2:length(y)], y[2:length(y)])),
    0
  )
  polygon(x = data.x, y = data.y, lwd = lwd, border = color, col = alpha(color, fill.alpha))
}

# four-weekly
filledLine(data = tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk, lwd = border.4wk, color = color.4wk, fill.alpha = fill.alpha.4wk)
# bi-weekly
filledLine(data = tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk, lwd = border.2wk, color = color.2wk, fill.alpha = fill.alpha.2wk)
# weekly
filledLine(data = tychoL1.measles.CA.cases.imp.zoo.plotwindow, lwd = border.1wk, color = color.1wk, fill.alpha = fill.alpha.1wk)


## Axes

axis(
  side = 1, # 1 specifies bottom axis - major ticks
  at = seq(plotxmin, plotxmax, "year"), # tick locations (in data units).
  labels = year(seq(plotxmin, plotxmax, "year")),
  tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
  line = .5, # axis position, in lines
  lty = "solid",
  lwd = 0, # axis line weight
  lwd.ticks = 1, # tick line weight
  col = "grey" # axis line color
)

axis(
  side = 2, # 1 specifies left axis
  at = cases.ticks, # tick locations (in data units).
  labels = cases.ticks, # vector of tick labels
  las = 1, # label orientation (0:parall., 1:horiz., 2:perp., 3:vert.)
  tcl = -.5, # tick length (fraction of plot width, neg to draw outside plot)
  line = .5, # axis position, in lines
  lty = "solid",
  lwd = 0, # axis line weight
  lwd.ticks = 1, # tick line weight
  col = "grey" # axis line color
)

## Axis Labels
title(ylab = "Number", line = 4)
title(xlab = "Year", line = 2.5)

# Close PDF ---------------------------------------------------------------

# %%
invisible(dev.off())


# FOOTER ------------------------------------------------------------------

# %%
## reset Timezone to previous stored local timezone
Sys.setenv(TZ = localTZ)
