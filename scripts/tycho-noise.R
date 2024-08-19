# %%
library(spaero)
library(tidyverse)
library(here)

tychoL1 <- read.csv(here::here("data", "Project_Tycho___Level_1_Data_20240813.csv"))

tycho_CA_measles_long_plotdata <- read_csv(here::here("out", "tycho_CA_measles_long_plotdata.csv"))
tycho_CA_measles_wide_plotdata <- read_csv(here::here("out", "tycho_CA_measles_wide_plotdata.csv"))

tycho_CA_measles_long_statsdata <- read_csv(here::here("out", "tycho_CA_measles_long_statsdata.csv"))
tycho_CA_measles_wide_statsdata <- read_csv(here::here("out", "tycho_CA_measles_wide_statsdata.csv"))

# %%
average_incidence <- tycho_CA_measles_wide_plotdata %>%
    summarize(across(.cols = contains("cases"), .fns = ~mean(.x, na.rm = TRUE)))

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

filled_aggregate_plotdata %>%
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
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  labs(x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation") +
  theme_minimal()

# %%
add_noise <- function(incidence_df, mean_incidence_df = NULL, aggregation_weeks = 1, noise_pc = 100) {
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

set.seed(1234)

noise_df <- add_noise(tycho_CA_measles_wide_plotdata, average_incidence, aggregation_weeks = 1, noise_pc = 100) %>%
  add_noise(average_incidence, aggregation_weeks = 2, noise_pc = 100) %>%
  add_noise(average_incidence, aggregation_weeks = 4, noise_pc = 100) %>%
  fill_aggregate_cases()

# %%
noise_long_df <- noise_df %>%
  pivot_longer(
    cols = c(-date),
    names_to = c("type", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "values",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk")))

noise_long_df

# %%
noise_long_df %>%
  ggplot(
    aes(x = date, y = values, color = aggregation, fill = aggregation)
  ) +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  labs(x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation") +
  theme_minimal()

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
  false_positive <- round(noise_tested * (1-test_chars$specificity), 0)

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
  calculate_test_characteristics(aggregation_weeks = 4, perfect_test_chars) %>%
  fill_aggregate_cases()

# Quickly check things are being calculated as expected
sum(perfect_test_df$test_pos_1wk == perfect_test_df$cases_1wk) == nrow(perfect_test_df)
sum(perfect_test_df$test_neg_1wk == perfect_test_df$noise_1wk) == nrow(perfect_test_df)
sum(perfect_test_df$total_tested_1wk == perfect_test_df$obs_1wk) == nrow(perfect_test_df)
sum(perfect_test_df$test_pos_2wk == perfect_test_df$cases_2wk) == nrow(perfect_test_df)
sum(perfect_test_df$test_neg_2wk == perfect_test_df$noise_2wk) == nrow(perfect_test_df)
sum(perfect_test_df$total_tested_2wk == perfect_test_df$obs_2wk) == nrow(perfect_test_df)
sum(perfect_test_df$test_pos_4wk == perfect_test_df$cases_4wk) == nrow(perfect_test_df)
sum(perfect_test_df$test_neg_4wk == perfect_test_df$noise_4wk) == nrow(perfect_test_df)
sum(perfect_test_df$total_tested_4wk == perfect_test_df$obs_4wk) == nrow(perfect_test_df)

# %%
perfect_test_long_df <- perfect_test_df %>%
  pivot_longer(
    cols = c(-date),
    names_to = c("type", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "values",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk")))

# %%
perfect_test_long_df %>%
  ggplot(
    aes(x = date, y = values, color = aggregation, fill = aggregation)
  ) +
  geom_line() +
  scale_color_manual(values = plot_colors, aesthetics = c("color", "fill")) +
  geom_area(aes(alpha = aggregation), position = "identity") +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_alpha_manual(values = c(0.0, 0.1, 1.0)) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  labs(x = "Date", y = "Incidence", color = "Aggregation", fill = "Aggregation", alpha = "Aggregation") +
  theme_minimal()


# %%
rdt_8080_test_chars <- list(prop_tested = 1.0, sensitivity = 0.8, specificity = 0.8)
rdt_8080_test_df <- calculate_test_characteristics(noise_df, aggregation_weeks = 1, rdt_8080_test_chars) %>%
  calculate_test_characteristics(aggregation_weeks = 2, rdt_8080_test_chars) %>%
  calculate_test_characteristics(aggregation_weeks = 4, rdt_8080_test_chars) %>%
  fill_aggregate_cases()

# %%
rdt_8080_test_long_df <- rdt_8080_test_df %>%
  pivot_longer(
    cols = c(-date),
    names_to = c("type", "aggregation"),
    names_pattern = "(.+)_(.+)",
    values_to = "values",
  ) %>%
  mutate(aggregation = factor(aggregation, levels = c("4wk", "2wk", "1wk")))

# %%
rdt_8080_test_long_df %>%
  filter(type %in% c("total_tested", "test_pos", "test_neg", "cases")) %>%
  ggplot(
    aes(x = date, y = values, color = type)
  ) +
  geom_line() +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  facet_wrap(~aggregation, scales = "free_y", ncol = 1) +
  labs(x = "Date", y = "Incidence", color = "Type") +
  theme_minimal()

