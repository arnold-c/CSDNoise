library(spaero)
library(tidyverse)

inc_array <- read_csv(here::here("out", "incidence-array_sim-1.csv"))

ews_stats <- get_stats(inc_array[, "incidence"],
  center_kernel = "uniform",
  center_trend = "local_constant", center_bandwidth = 360,
  stat_bandwidth = 360
)

ewsmetrics <-
  ews_stats$stats %>%
  as_tibble() %>%
  mutate(time = row_number()) %>%
  select(time, everything()) %>%
  pivot_longer(cols = c(everything(), -time))

ggplot(ewsmetrics, aes(x = time, y = value, color = name)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y")
