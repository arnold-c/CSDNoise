library(spaero)
library(tidyverse)

inc_array <- read_csv(here::here("out", "incidence-array_sim-1.csv"))

spaero_inc_ews_stats <- get_stats(
  inc_array[, "incidence"],
  center_bandwidth = 30,
  stat_bandwidth = 30,
  center_trend = "local_constant",
  stat_trend = "local_constant",
  center_kernel = "uniform",
  stat_kernel = "uniform",
  lag = 1,
  # backward_only = FALSE
  backward_only = TRUE
)$stats

spaero_ensemble_single_inc_ews_arr <- spaero_inc_ews_stats %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything())
