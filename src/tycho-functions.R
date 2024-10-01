library(spaero)
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
