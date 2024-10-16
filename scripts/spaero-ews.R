library(spaero)
library(tidyverse)

incvec_1d <- read_csv(here::here("out", "tycho", "tests", "incidence-array_1d_sim-1.csv"), col_names = FALSE)
incvec_30d <- read_csv(here::here("out", "tycho", "tests", "incidence-array_30d_sim-1.csv"), col_names = FALSE)

get_stats_wrapper <- function(incvec, bandwidth = 35, backward_only = FALSE, lag = 1) {
  get_stats(
    incvec,
    center_bandwidth = bandwidth,
    stat_bandwidth = bandwidth,
    center_trend = "local_constant",
    stat_trend = "local_constant",
    center_kernel = "uniform",
    stat_kernel = "uniform",
    lag = lag,
    backward_only = backward_only
  )$stats
}

spaero_ews_centered_1d <- get_stats_wrapper(
	incvec_1d,
  backward_only = FALSE,
) %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())

spaero_ews_centered_30d <- get_stats_wrapper(
	incvec_30d,
  backward_only = FALSE,
) %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())

spaero_ews_backward_1d <- get_stats_wrapper(
	incvec_1d,
  backward_only = TRUE,
) %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())

spaero_ews_backward_30d <- get_stats_wrapper(
	incvec_30d,
  backward_only = TRUE,
) %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())

