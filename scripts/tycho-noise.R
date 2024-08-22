# %%
library(spaero)
library(tidyverse)
library(here)
library(patchwork)

theme_set(
  theme_minimal() +
    theme(
      title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )
)

# %%
source(here::here("src", "tycho-functions.R"))

bandwidth_weeks <- 52
bandwidth <- bandwidth_weeks * 7 # days

obsdate <- cdcweekToDate(year = 1990, week = 34 - 31, weekday = 6)
plotxmin.year <- 1984
plotxmin <- cdcweekToDate(year = plotxmin.year, week = 1, weekday = 6)
plotxmax.year <- 1992
plotxmax <- cdcweekToDate(year = plotxmax.year, week = 1, weekday = 6) + bandwidth

# %%
tychoL1 <- read.csv(here::here("data", "Project_Tycho___Level_1_Data_20240813.csv"))

tycho_CA_measles_long_plotdata <- read_csv(here::here("out", "tycho_CA_measles_long_plotdata.csv"))
tycho_CA_measles_wide_plotdata <- read_csv(here::here("out", "tycho_CA_measles_wide_plotdata.csv"))

tycho_CA_measles_long_statsdata <- read_csv(here::here("out", "tycho_CA_measles_long_statsdata.csv"))
tycho_CA_measles_wide_statsdata <- read_csv(here::here("out", "tycho_CA_measles_wide_statsdata.csv"))

# %%
# confirm all data is weekly aggregated
unique(diff(tycho_CA_measles_wide_plotdata$date)) == 7

# %%
average_incidence <- tycho_CA_measles_wide_plotdata %>%
  summarize(across(.cols = contains("cases"), .fns = ~ mean(.x, na.rm = TRUE)))

# %%
mean_inc_noise <- tycho_CA_measles_long_plotdata %>%
  group_by(aggregation_weeks) %>%
  summarize(mean_cases = mean(cases, na.rm = TRUE)) %>%
  mutate(
    noise_10pc = round(mean_cases * 0.1, 0),
    noise_25pc = round(mean_cases * 0.25, 0),
    noise_50pc = round(mean_cases * 0.5, 0)
  )

# %%
fill_aggregate_cases <- function(data) {
  start_date <- min(data$date)
  end_date <- max(data$date)

  tibble(
    date = seq.Date(from = start_date, to = end_date, by = 1)
  ) %>%
    left_join(data, by = "date") %>%
    fill(c(-date), .direction = "down")
}

filled_aggregate_plotdata <- fill_aggregate_cases(tycho_CA_measles_wide_plotdata)
filled_aggregate_statsdata <- fill_aggregate_cases(tycho_CA_measles_wide_statsdata)

# %%
# Recreate incidence plot with ggplot2
plot_colors <- c("1wk" = "gray20", "2wk" = "blue", "4wk" = "darkred")

incidence_plot <- filled_aggregate_plotdata %>%
  pivot_longer(
    cols = contains("cases"),
    names_to = "aggregation",
    names_pattern = "cases_(.+)",
    values_to = "cases",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk"))) %>%
  ggplot(
    aes(x = date, y = cases, color = aggregation, fill = aggregation)
  ) +
  annotate("rect", xmin = obsdate, xmax = plotxmax, ymin = 0, ymax = Inf, alpha = 0.1, fill = "red") +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y", limits = c(as.Date(plotxmin), as.Date(plotxmax))) +
  labs(title = "Actual Incidence", x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation")

incidence_plot

# %%
add_noise <- function(incidence_df, mean_incidence_df, mean_incidence = NULL, aggregation_weeks = 1, noise_pc = 100) {
  cases_col <- paste0("cases_", aggregation_weeks, "wk")
  noise_col <- paste0("noise_", aggregation_weeks, "wk")
  obs_col <- paste0("obs_", aggregation_weeks, "wk")

  if (is.null(mean_incidence)) {
    mean_incidence <- mean(incidence_df[[cases_col]], na.rm = TRUE)
  }

  noise_prop <- noise_pc / 100
  incidence_df[[noise_col]] <- rpois(n = length(incidence_df[[cases_col]]), lambda = round(noise_prop * mean_incidence_df[[cases_col]], 0))
  incidence_df[[noise_col]] <- if_else(is.na(incidence_df[[cases_col]]), NA, incidence_df[[noise_col]])
  incidence_df[[obs_col]] <- incidence_df[[cases_col]] + incidence_df[[noise_col]]

  return(incidence_df)
}

# %%
set.seed(1234)

noise_df <- add_noise(tycho_CA_measles_wide_plotdata, average_incidence, aggregation_weeks = 1, noise_pc = 100) %>%
  add_noise(average_incidence, aggregation_weeks = 2, noise_pc = 100) %>%
  add_noise(average_incidence, aggregation_weeks = 4, noise_pc = 100)


filled_noise_df <- fill_aggregate_cases(noise_df)

# %%
filled_noise_long_df <- filled_noise_df %>%
  pivot_longer(
    cols = c(-date),
    names_to = c("type", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "values",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk")))

filled_noise_long_df

# %%
filled_noise_long_df %>%
  mutate(type = factor(type, levels = c("obs", "noise", "cases"))) %>%
  ggplot(
    aes(x = date, y = values, color = aggregation, fill = aggregation)
  ) +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  labs(x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation")

# %%
calculate_test_characteristics <- function(inc_noise_df, aggregation_weeks = 1, test_chars = list(prop_tested = 1.0, sensitivity = 1.0, specificity = 1.0)) {
  cases_col <- paste0("cases_", aggregation_weeks, "wk")
  noise_col <- paste0("noise_", aggregation_weeks, "wk")
  total_tested_col <- paste0("total_tested_", aggregation_weeks, "wk")
  test_pos_col <- paste0("test_pos_", aggregation_weeks, "wk")
  test_neg_col <- paste0("test_neg_", aggregation_weeks, "wk")

  cases_tested <- round(inc_noise_df[[cases_col]] * test_chars$prop_tested, 0)
  noise_tested <- round(inc_noise_df[[noise_col]] * test_chars$prop_tested, 0)

  true_positives <- round(cases_tested * test_chars$sensitivity, 0)
  false_positive <- round(noise_tested * (1 - test_chars$specificity), 0)

  false_negatives <- cases_tested - true_positives
  true_negatives <- noise_tested - false_positive

  inc_noise_df[[total_tested_col]] <- cases_tested + noise_tested
  inc_noise_df[[test_pos_col]] <- true_positives + false_positive
  inc_noise_df[[test_neg_col]] <- true_negatives + false_negatives

  return(inc_noise_df)
}

# %%
perfect_test_chars <- list(prop_tested = 1.0, sensitivity = 1.0, specificity = 1.0)
perfect_test_df <- calculate_test_characteristics(noise_df, aggregation_weeks = 1, perfect_test_chars) %>%
  calculate_test_characteristics(aggregation_weeks = 2, perfect_test_chars) %>%
  calculate_test_characteristics(aggregation_weeks = 4, perfect_test_chars)

filled_perfect_test_df <- fill_aggregate_cases(perfect_test_df)

# Quickly check things are being calculated as expected
sum(filled_perfect_test_df$test_pos_1wk == filled_perfect_test_df$cases_1wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$test_neg_1wk == filled_perfect_test_df$noise_1wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$total_tested_1wk == filled_perfect_test_df$obs_1wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$test_pos_2wk == filled_perfect_test_df$cases_2wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$test_neg_2wk == filled_perfect_test_df$noise_2wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$total_tested_2wk == filled_perfect_test_df$obs_2wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$test_pos_4wk == filled_perfect_test_df$cases_4wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$test_neg_4wk == filled_perfect_test_df$noise_4wk) == nrow(filled_perfect_test_df)
sum(filled_perfect_test_df$total_tested_4wk == filled_perfect_test_df$obs_4wk) == nrow(filled_perfect_test_df)

# %%
filled_perfect_test_long_df <- filled_perfect_test_df %>%
  pivot_longer(
    cols = c(-date),
    names_to = c("type", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "values",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk")))

# %%
filled_perfect_test_long_df %>%
  mutate(
    type = factor(
      type,
      levels = c("test_pos", "test_neg", "total_tested", "obs", "noise", "cases")
    )
  ) %>%
  ggplot(
    aes(x = date, y = values, color = aggregation, fill = aggregation)
  ) +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y", limits = c(as.Date(plotxmin), as.Date(plotxmax))) +
  labs(title = "Perfect Test Characteristics", x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation")

# %%
rdt_8080_test_chars <- list(prop_tested = 1.0, sensitivity = 0.8, specificity = 0.8)
rdt_8080_test_df <- calculate_test_characteristics(noise_df, aggregation_weeks = 1, rdt_8080_test_chars) %>%
  calculate_test_characteristics(aggregation_weeks = 2, rdt_8080_test_chars) %>%
  calculate_test_characteristics(aggregation_weeks = 4, rdt_8080_test_chars)

filled_rdt_8080_test_df <- fill_aggregate_cases(rdt_8080_test_df)

# %%
filled_rdt_8080_test_long_df <- filled_rdt_8080_test_df %>%
  pivot_longer(
    cols = c(-date),
    names_to = c("type", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "values",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk")))

# %%
filled_rdt_8080_test_long_df %>%
  filter(type %in% c("total_tested", "test_pos", "test_neg", "cases")) %>%
  mutate(type = factor(type, levels = c("test_pos", "test_neg", "total_tested", "cases"))) %>%
  ggplot(
    aes(x = date, y = values, color = aggregation, fill = aggregation)
  ) +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  labs(title = "RDT 80/80, 100% Testing", x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation")

# %%
analysis_params_1w <- list(
  center_trend = "local_constant",
  stat_trend = "local_constant",
  center_kernel = "uniform",
  stat_kernel = "uniform",
  center_bandwidth = bandwidth_weeks,
  stat_bandwidth = bandwidth_weeks,
  lag = 1
)

analysis_params_4w <- analysis_params_2w <- analysis_params_1w
analysis_params_2w$center_bandwidth <- analysis_params_2w$stat_bandwidth <- bandwidth_weeks / 2
analysis_params_4w$center_bandwidth <- analysis_params_4w$stat_bandwidth <- bandwidth_weeks / 4

# %%
filter_statdata <- function(data, aggregation_wks = 1, inc_var = "test_pos", obsdate) {
  analysis_params <- list(
    center_trend = "local_constant",
    stat_trend = "local_constant",
    center_kernel = "uniform",
    stat_kernel = "uniform",
    center_bandwidth = bandwidth_weeks / aggregation_wks,
    stat_bandwidth = bandwidth_weeks / aggregation_wks,
    lag = 1
  )

  inc_col <- paste0(inc_var, "_", aggregation_wks, "wk")
  statdata <- filter(data, date <= obsdate) %>%
    select(date, inc_col) %>%
    drop_na()

  return(list(dates = statdata$date, analysis_params = analysis_params, ews = analysis(statdata[[inc_col]], analysis_params)))
}

rdt_8080_test_pos_stats <- c(1, 2, 4) %>%
  rlang::set_names(., paste0, "wk") %>%
  purrr::map(
    .,
    ~ filter_statdata(rdt_8080_test_df, aggregation_wks = .x, obsdate = obsdate)
  )

cases_stats <- c(1, 2, 4) %>%
  rlang::set_names(., paste0, "wk") %>%
  purrr::map(
    .,
    ~ filter_statdata(rdt_8080_test_df, aggregation_wks = .x, inc_var = "cases", obsdate = obsdate)
  )

# %%
# Check stats
tycho_CA_measles_stats_1wk <- readRDS(here::here("out", "tycho_CA_measles_stats_1wk.rds" ))
tycho_CA_measles_stats_2wk <- readRDS(here::here("out", "tycho_CA_measles_stats_2wk.rds" ))
tycho_CA_measles_stats_4wk <- readRDS(here::here("out", "tycho_CA_measles_stats_4wk.rds" ))

sum(tycho_CA_measles_stats_1wk$stats$variance == cases_stats[["1wk"]]$ews$stats$variance) == length(tycho_CA_measles_stats_1wk$stats$variance)

# %%
rdt_8080_test_pos_plot <- filled_rdt_8080_test_long_df %>%
  filter(type == "test_pos") %>%
  ggplot(
    aes(x = date, y = values, color = aggregation, fill = aggregation)
  ) +
  annotate("rect", xmin = obsdate, xmax = plotxmax, ymin = 0, ymax = Inf, alpha = 0.1, fill = "red") +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y", limits = c(as.Date(plotxmin), as.Date(plotxmax))) +
  labs(title = "Test Positive: RDT 80/80, 100% Testing", x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation")

rdt_8080_test_pos_plot

# %%
rdt_8080_test_pos_plot / incidence_plot

# %%
extract_ews_to_df <- function(ews_list, aggregation_weeks = 1) {
  aggregation <- paste0(aggregation_weeks, "wk")
  bind_cols(
    date = ews_list[[aggregation]]$dates,
    as_tibble(ews_list[[aggregation]]$ews$stats)
  ) %>%
     rename_with(
       .cols = -date,
       .fn = ~ paste0(.x, "_", aggregation)
     )
}


# %%
rdt_8080_ews_df <- purrr::map(
  c(1, 2, 4),
  ~ extract_ews_to_df(rdt_8080_test_pos_stats, aggregation_weeks = .x)
) %>%
  reduce(left_join, by = "date")

rdt_8080_ews_long_df <- rdt_8080_ews_df %>%
  pivot_longer(
    cols = -date,
    names_to = c("statistic", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "value"
  )

# %%
cases_ews_df <- purrr::map(
  c(1, 2, 4),
  ~ extract_ews_to_df(cases_stats, aggregation_weeks = .x)
) %>%
  reduce(left_join, by = "date")

cases_ews_long_df <- cases_ews_df %>%
  pivot_longer(
    cols = -date,
    names_to = c("statistic", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "value"
  )

# %%
# check stats df
sum(tycho_CA_measles_stats_1wk$stats$variance == cases_ews_df$variance_1wk) == length(tycho_CA_measles_stats_1wk$stats$variance)
sum(cases_stats[["1wk"]]$ews$stats$variance == cases_ews_df$variance_1wk) == length(cases_stats[["1wk"]]$ews$stats$variance)

sum(tycho_CA_measles_stats_2wk$stats$variance == cases_ews_df$variance_2wk[!is.na(cases_ews_df$variance_2wk)]) == length(tycho_CA_measles_stats_2wk$stats$variance)
sum(cases_stats[["2wk"]]$ews$stats$variance == cases_ews_df$variance_2wk[!is.na(cases_ews_df$variance_2wk)]) == length(cases_stats[["2wk"]]$ews$stats$variance)

sum(tycho_CA_measles_stats_4wk$stats$variance == cases_ews_df$variance_4wk[!is.na(cases_ews_df$variance_4wk)]) == length(tycho_CA_measles_stats_4wk$stats$variance)
sum(cases_stats[["4wk"]]$ews$stats$variance == cases_ews_df$variance_4wk[!is.na(cases_ews_df$variance_4wk)]) == length(cases_stats[["4wk"]]$ews$stats$variance)

# %%
extract_tau <- function(ews_list, aggregation = "1wk") {
  tibble(
    aggregation = aggregation,
    statistic = names(ews_list$ews$taus),
    value = as.numeric(ews_list$ews$taus)
  )
}

rdt_8080_ews_tau_df <- map2_dfr(
  names(rdt_8080_test_pos_stats),
  rdt_8080_test_pos_stats,
  function(.x, .y) {
    extract_tau(.y, aggregation = .x)
  }
)

cases_ews_tau_df <- map2_dfr(
  names(cases_stats),
  cases_stats,
  function(.x, .y) {
    extract_tau(.y, aggregation = .x)
  }
)

# %%
# check tau values
tycho_CA_measles_stats_1wk$taus$variance == filter(cases_ews_tau_df, aggregation == "1wk", statistic == "variance")$value
tycho_CA_measles_stats_2wk$taus$variance == filter(cases_ews_tau_df, aggregation == "2wk", statistic == "variance")$value
tycho_CA_measles_stats_4wk$taus$variance == filter(cases_ews_tau_df, aggregation == "4wk", statistic == "variance")$value

# %%
ews_plot <- function(ews_df, tau_df, ews = variance, ews_type = "Test Positive", test_type = "RDT 80/80, 100% Testing") {
  test_title <- ifelse(is.null(test_type), "", paste0(": ", str_to_title(test_type)))
  plot <- ews_df %>%
    filter(statistic == {{ ews }}) %>%
    drop_na() %>%
    ggplot(
      aes(x = date, y = value, color = aggregation)
    ) +
    geom_line(aes(group = aggregation), na.rm = TRUE) +
    scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
    scale_x_date(date_breaks = "1 years", date_labels = "%Y", limits = c(as.Date(plotxmin), as.Date(plotxmax))) +
    labs(
      title = paste0(str_to_title(ews_type), " ", str_to_title(ews), test_title),
      x = "Date",
      y = str_to_title(ews),
      color = "Aggregation"
    )

  tau_label_df <- ews_df %>%
    filter(statistic == {{ ews }}) %>%
    drop_na() %>%
    group_by(aggregation) %>%
    filter(date == last(date)) %>%
    transmute(
      date = date + weeks(4),
      y = value,
      aggregation = aggregation
    )

  taus <- tau_df %>%
    filter(statistic == {{ ews }}) %>%
    left_join(
      .,
      tau_label_df,
      by = "aggregation"
    )

  plot <- plot +
    geom_text(
      data = taus,
      aes(
        x = date,
        y = y,
        label = TeX(sprintf("$\\tau = $ %.2f", value), output = "character")
      ),
      parse = TRUE,
      hjust = "left",
      show.legend = FALSE
    )

  return(plot)
}

# %%
rdt_8080_var_plot <- ews_plot(rdt_8080_ews_long_df, rdt_8080_ews_tau_df, ews = "variance")
rdt_8080_var_plot

# %%
rdt_8080_var_test_pos_plot <- rdt_8080_var_plot / rdt_8080_test_pos_plot + plot_layout(axes = "collect")

rdt_8080_var_test_pos_plot

ggsave(rdt_8080_var_test_pos_plot, filename = here::here("plots", "tycho", "rdt_8080_var_test_pos_plot.png"), width = 10, height = 5, dpi = 600)

# %%
cases_var_plot <- ews_plot(cases_ews_long_df, cases_ews_tau_df, ews = "variance", ews_type = "Actual Incidence", test_type = NULL)
cases_var_plot

# %%
cases_var_incidence_plot <- cases_var_plot / incidence_plot + plot_layout(axes = "collect")

cases_var_incidence_plot

ggsave(cases_var_incidence_plot, filename = here::here("plots", "tycho", "cases_var_incidence_plot.png"), width = 10, height = 5, dpi = 600)

ggsave(cases_var_incidence_plot, filename = here::here("plots", "tycho", "cases_var_incidence_plot.png"), width = 10, height = 5, dpi = 600)
