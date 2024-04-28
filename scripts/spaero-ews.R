library(spaero)
library(tidyverse)

inc_array <- read_csv(here::here("out", "incidence-array_sim-1.csv"))

get_stats_wrapper <- function(incvec, bandwidth = 30, backward_only = FALSE) {
  get_stats(
    incvec,
    center_bandwidth = bandwidth,
    stat_bandwidth = bandwidth,
    center_trend = "local_constant",
    stat_trend = "local_constant",
    center_kernel = "uniform",
    stat_kernel = "uniform",
    lag = 1,
    backward_only = backward_only
  )$stats
}

spaero_ews_centered <- get_stats_wrapper(
  inc_array[, "incidence"],
  backward_only = FALSE
) %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())

spaero_ews_backward <- get_stats_wrapper(
  inc_array[, "incidence"],
  backward_only = TRUE
) %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())
