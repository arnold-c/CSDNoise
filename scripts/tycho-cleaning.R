# %%
# all packages used in ericfigs/
library(utils) # standard R package
library(stats) # standard R package
library(scales)
library(zoo)
library(reshape)
library(magrittr)
library(pomp) # required by spaero
library(spaero) # model simulation and statistics
library(here)
library(tidyverse)

## Additional packages for managing real data
library(lubridate)
library(padr)
library(imputeTS)


# Helper functions --------------------------------------------------------
source(here::here("src", "tycho-functions.R"))

# %%
## Sample process
sample_process <- function(external_forcing = 1 / 7, host_lifetime = 70 * 365,
                           infectious_days = 7, observation_days = 20 * 365,
                           sampling_interval = 1, population_size = 1e6,
                           process_reps = 1) {
  ## Runs model and returns sample of state variables and flow from I to R, which is the cases column
  ## Default parameters are those used in the simulation study.

  times <- seq(0, observation_days, by = sampling_interval)
  params <- c(
    gamma = 1 / infectious_days, mu = 1 / host_lifetime,
    d = 1 / host_lifetime, eta = external_forcing / population_size,
    beta = 0, rho = 0.1, S_0 = 1, I_0 = 0, R_0 = 0, N_0 = population_size
  )
  beta_final <- (params["gamma"] + params["d"]) / population_size
  covalt <- data.frame(
    gamma_t = c(0, 0), mu_t = c(0, 0), d_t = c(0, 0),
    eta_t = c(0, 0), beta_t = c(0, beta_final),
    time = c(0, observation_days)
  )

  simtest <- spaero::create_simulator(times = times, params = params, covar = covalt)

  ret <- list()
  do_sim <- function(obj, nsim = process_reps) {
    cols_to_delete <- c("reports", "gamma_t", "mu_t", "d_t", "eta_t")
    ret <- pomp::simulate(obj, nsim = nsim, as.data.frame = TRUE)
    ret[, !colnames(ret) %in% cols_to_delete]
  }
  do_sim(simtest)
}

# %%
## Sample observations taking a time series as input
sample_observation <- function(cases, sampling_interval = 1, tau = 1, reporting_prob = 1, dispersion_parameter = 100) {
  tots <- aggregate.ts(cases, nfrequency = 1 / (tau / sampling_interval))
  mu <- tots * reporting_prob
  n <- length(mu)
  sampled <- rnbinom(n = n, mu = mu, size = dispersion_parameter)
  ts(sampled, start = start(tots), end = end(tots), frequency = frequency(tots))
}

# %%
## Statistical Analysis
analysis <- function(data, params) {
  get_stats(
    data,
    center_trend = params$center_trend,
    stat_trend = params$stat_trend,
    center_kernel = params$center_kernel,
    stat_kernel = params$stat_kernel,
    center_bandwidth = params$center_bandwidth,
    stat_bandwidth = params$stat_bandwidth,
    lag = params$lag
  )
}

# %%
## CDC Epiweek to Dates
## returns the start date of the CDC epiweek
cdcweekToDate <- function(epiweek, weekday = 0, year = NULL, week = NULL) {
  if (missing(epiweek)) {
    year <- year
    week <- week
  } else {
    year <- epiweek %/% 1e2
    week <- epiweek %% 1e2
  }
  jan1 <- as.Date(ISOdate(year, 1, 1))
  jan1.wday <- as.POSIXlt(jan1)$wday
  if (jan1.wday < 4) {
    origin <- jan1 - jan1.wday
  } else {
    origin <- jan1 + (7 - jan1.wday)
  }
  date <- origin + (week - 1) * 7 + weekday
  return(date)
}

# %%
tmpfcn <- function(x, agg.interval = 1) {
  index <- index(x)
  interval <- mean(diff(index))
  agg.groups <- rep(seq(index[1], index[length(index)], interval * agg.interval), each = agg.interval)
  agg.groups <- agg.groups[1:length(x)]
  out <- aggregate(x, by = agg.groups, sum)
  return(out)
}


# %%
## Set UTC time
localTZ <- Sys.timezone()
Sys.setenv(TZ = "UTC")

# DATA --------------------------------------------------------------------

# Project Tycho level 1 data

## local copy of Project Tycho level 1 data from:
# browseURL:("https://www.healthdata.gov/sites/default/files/ProjectTycho_Level1_v1.0.0.csv")
# Metadata:
# browseURL('https://catalog.data.gov/harvest/object/cb73ca20-e127-4b96-8a48-646d0d4a606f')

## tychoL1 <- read.csv('data/ProjectTycho_Level1_v1.0.0.csv', stringsAsFactors = F, na.strings = "\\N")
# %%
tychoL1 <- read.csv(here::here("data", "Project_Tycho___Level_1_Data_20240813.csv"))

# The original copy of Project Tycho level 1 data may be loaded directly using the line of code below:
# tychoL1 <- read.csv(url("https://www.healthdata.gov/sites/default/files/ProjectTycho_Level1_v1.0.0.csv", stringsAsFactors = F, na.strings = "\\N"))

# %%
## Add column for epiweek enddate; calculate epiweek enddate from epiweek
# Note: cdcweekToDate() calculates the start date from the 6 digit epiweek code (yyyyww)
# Note: cdcweekToDate() is not a vectorized function and may take several minutes to run on a large dataset.
tychoL1 <- cbind(enddate = NA, tychoL1)
tychoL1$enddate <- do.call(
  "c",
  lapply(tychoL1$epi_week, cdcweekToDate, weekday = 6)
)

# %%
## Extract CA Measles cases and pad missing weeks with NA's
tychoL1.measles.CA <- pad(
  tychoL1[(tychoL1$disease == "MEASLES") & (tychoL1$state == "CA"), c("enddate", "epi_week", "cases")],
  interval = "week"
)

# %%
## Impute cases with Kalman filter
tychoL1.measles.CA$cases.imp <- round(na.kalman(as.numeric(tychoL1.measles.CA$cases)))

# %%
## Convert to zoo time series
tychoL1.measles.CA.cases.imp.zoo <- zoo(tychoL1.measles.CA$cases.imp, tychoL1.measles.CA$enddate)

# Set windows -----------------------------------------------------

# %%
## Statistics bandwidth
bandwidth_weeks <- 52
bandwidth <- bandwidth_weeks * 7 # days

# %%
## date of peak of epidemic
peakdate <- cdcweekToDate(year = 1990, week = 34, weekday = 6)

# %%
## date prior to peak of outbreak from which past data is used to calculate statistics.
obsdate <- cdcweekToDate(year = 1990, week = 34 - 31, weekday = 6)

# %%
## end date for plot of stats
statsend <- obsdate

# %%
## set range of years for overall plat
plotxmin.year <- 1984
plotxmax.year <- 1992

# %%
## cdcweekToDate() calculates the date from the epiweek (year and week)
plotxmin <- cdcweekToDate(year = plotxmin.year, week = 1, weekday = 6)
plotxmax <- cdcweekToDate(year = plotxmax.year, week = 1, weekday = 6)

save(plotxmin, plotxmax, obsdate, statsend, bandwidth_weeks, bandwidth, file = here::here("out", "tycho", "R-scripts", "tycho-params.RData"))

# Window data -------------------------------------------------------------

# %%
## window data to plot region
tychoL1.measles.CA.cases.imp.zoo.plotwindow <- window(tychoL1.measles.CA.cases.imp.zoo,
  start = plotxmin,
  end = plotxmax
)

# Aggregate biweekly and 4-weekly -----------------------------------------

# %%
# aggregate windowed dataset
tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo.plotwindow, 2)
tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk <- tmpfcn(tychoL1.measles.CA.cases.imp.zoo.plotwindow, 4)

saveRDS(tychoL1.measles.CA.cases.imp.zoo.plotwindow, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_plotwindow.rds"))
saveRDS(tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_plotwindow_2wk.rds"))
saveRDS(tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_plotwindow_4wk.rds"))

# %%
zoo_to_tibble <- function(zoo_vec, agg_weeks = 1) {
  tibble(
    date = index(zoo_vec),
    aggregation_weeks = agg_weeks,
    cases = as.numeric(zoo_vec)
  )
}

tycho_CA_measles_long_plotdata <- bind_rows(
  zoo_to_tibble(tychoL1.measles.CA.cases.imp.zoo.plotwindow),
  zoo_to_tibble(tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk, agg_weeks = 2),
  zoo_to_tibble(tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk, agg_weeks = 4)
)

tycho_CA_measles_wide_plotdata <- tycho_CA_measles_long_plotdata %>%
  pivot_wider(
    names_from = aggregation_weeks,
    values_from = cases,
    names_glue = "{.value}_{aggregation_weeks}wk"
  )

write_csv(tycho_CA_measles_long_plotdata, here::here("out", "tycho", "tycho_CA_measles_long_plotdata.csv"))
write_csv(tycho_CA_measles_wide_plotdata, here::here("out", "tycho", "tycho_CA_measles_wide_plotdata.csv"))

tycho_CA_measles_long_statsdata <- tycho_CA_measles_long_plotdata %>%
  filter(date <= obsdate)

tycho_CA_measles_wide_statsdata <- tycho_CA_measles_wide_plotdata %>%
  filter(date <= obsdate)

write_csv(tycho_CA_measles_long_statsdata, here::here("out", "tycho", "tycho_CA_measles_long_statsdata.csv"))
write_csv(tycho_CA_measles_wide_statsdata, here::here("out", "tycho", "tycho_CA_measles_wide_statsdata.csv"))

# Window aggregated data for analysis -------------------------------------

## window to end at a chosen "observation date."
# This windowed data is used below to calculate statistics

# %%
tychoL1.measles.CA.cases.imp.zoo.statswindow <- window(tychoL1.measles.CA.cases.imp.zoo.plotwindow,
  end = obsdate
)
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk <- window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.2wk,
  end = obsdate
)
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk <- window(tychoL1.measles.CA.cases.imp.zoo.plotwindow.4wk,
  end = obsdate
)

saveRDS(tychoL1.measles.CA.cases.imp.zoo.statswindow, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_statswindow_1wk.rds"))
saveRDS(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_statswindow_2wk.rds"))
saveRDS(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_statswindow_4wk.rds"))


# Get Stats ---------------------------------------------------------------

# %%
analysis_params_1wk <- list(
  center_trend = "local_constant",
  stat_trend = "local_constant",
  center_kernel = "uniform",
  stat_kernel = "uniform",
  center_bandwidth = bandwidth_weeks,
  stat_bandwidth = bandwidth_weeks,
  lag = 1
)

# %%
# weekly reports
analysis_params_1wk$center_bandwidth <- analysis_params_1wk$stat_bandwidth <- bandwidth_weeks
tychoL1.measles.CA.cases.imp.zoo.statswindow.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow, analysis_params_1wk)

# %%
# bi-weekly reports
analysis_params_2wk <- analysis_params_1wk
analysis_params_2wk$center_bandwidth <- analysis_params_2wk$stat_bandwidth <- bandwidth_weeks / 2
tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk, analysis_params_2wk)

# %%
# four-weekly reports
analysis_params_4wk <- analysis_params_1wk
analysis_params_4wk$center_bandwidth <- analysis_params_4wk$stat_bandwidth <- bandwidth_weeks / 4
tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats <- analysis(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk, analysis_params_4wk)

# %%
saveRDS(tychoL1.measles.CA.cases.imp.zoo.statswindow.stats, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_stats_1wk.rds"))
saveRDS(tychoL1.measles.CA.cases.imp.zoo.statswindow.2wk.stats, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_stats_2wk.rds"))
saveRDS(tychoL1.measles.CA.cases.imp.zoo.statswindow.4wk.stats, here::here("out", "tycho", "R-scripts", "tycho_CA_measles_stats_4wk.rds"))


