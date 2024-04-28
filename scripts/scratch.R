# %%
library(spaero)
library(testthat)

# %%
sinx <- sin(1:10)
# sinx <- 1:10
bw2 <- 2
get_stats_wrapper <- function(x, bw) {
    get_stats(
        x,
        center_bandwidth = bw,
        stat_bandwidth = bw,
        center_trend = "local_constant",
        stat_trend = "local_constant",
        center_kernel = "uniform",
        stat_kernel = "uniform",
        lag = 1
    )$stats
}
stats <- get_stats_wrapper(sinx, bw2)

# %%
central_mean <- function(x, bw, na.rm = TRUE) {
    purrr::map_dbl(seq_along(x), function(ind) {
        if (ind < bw) {
            start_ind <- 1
            end_ind <- ind + bw - 1
        } else if (ind + bw > length(x)) {
            start_ind <- ind - bw + 1
            end_ind <- length(x)
        } else {
            start_ind <- ind - bw + 1
            end_ind <- ind + bw - 1
        }
        mean(x[start_ind:end_ind], na.rm = na.rm)
    })
}

stats$mean
central_mean(sinx, bw2)
expect_equal(stats$mean, central_mean(sinx, bw2))

get_stats_wrapper(sinx, 3)$mean
central_mean(sinx, 3)
expect_equal(get_stats_wrapper(sinx, 3)$mean, central_mean(sinx, 3))

# %%
central_moment <- function(x, bw, moment = 2, na.rm = TRUE) {
    central_mean((x - central_mean(x, bw, na.rm = na.rm))^moment, bw, na.rm = na.rm)
}

# %%
stats$variance
central_moment(sinx, bw2, moment = 2)
expect_equal(stats$variance, central_moment(sinx, bw2, moment = 2))

# %%
central_cov <- function(x, bw, na.rm = TRUE) {
    sqrt(central_moment(x, bw, moment = 2, na.rm = na.rm)) / central_mean(x, bw, na.rm = na.rm)
}
stats$coefficient_of_variation
central_cov(sinx, bw2)
expect_equal(stats$coefficient_of_variation, central_cov(sinx, bw2))

# %%
central_iod <- function(x, bw, na.rm = TRUE) {
    central_moment(x, bw, moment = 2, na.rm = na.rm) / central_mean(x, bw, na.rm = na.rm)
}
stats$index_of_dispersion
central_iod(sinx, bw2)
expect_equal(stats$index_of_dispersion, central_iod(sinx, bw2))

# %%
central_skewness <- function(x, bw, na.rm = TRUE) {
    central_moment(x, bw, moment = 3, na.rm = na.rm) / (central_moment(x, bw, moment = 2, na.rm = na.rm))^1.5
}
stats$skewness
central_skewness(sinx, bw2)
expect_equal(stats$skewness, central_skewness(sinx, bw2))

# %%
central_kurtosis <- function(x, bw, na.rm = TRUE) {
    central_moment(x, bw, moment = 4, na.rm = na.rm) / (central_moment(x, bw, moment = 2, na.rm = na.rm))^2
}
stats$kurtosis
central_kurtosis(sinx, bw2)
expect_equal(stats$kurtosis, central_kurtosis(sinx, bw2))

# %%
central_autocov <- function(x, bw, lag = 1, na.rm = TRUE) {
    meandiff <- x - central_mean(x, bw, na.rm = na.rm)
    val <- central_mean(meandiff * dplyr::lag(meandiff, lag), bw, na.rm = na.rm)
    val[1] <- NA
    val
}
stats$autocovariance
central_autocov(sinx, bw2)
expect_equal(stats$autocovariance, central_autocov(sinx, bw2))

# %%
central_autocor <- function(x, bw, lag = 1, na.rm = TRUE) {
    sd <- sqrt(central_moment(x, bw, moment = 2, na.rm = na.rm))
    central_autocov(x, bw, na.rm = na.rm) / (sd * dplyr::lag(sd, lag))
}
stats$autocorrelation
central_autocor(sinx, bw2)
expect_equal(stats$autocorrelation, central_autocor(sinx, bw2))
