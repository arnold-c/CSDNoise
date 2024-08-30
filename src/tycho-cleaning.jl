using Dates: Dates

function cdc_week_to_date(epiweek; weekday = 0)
    year = epiweek ÷ 100
    week = epiweek % 100

    return cdc_week_to_date(year, week; weekday = weekday)
end

function cdc_week_to_date(year, week; weekday = 0)
    jan1 = Dates.Date(year, 1, 1)
    jan1_wday = Dates.dayofweek(jan1)

    # Calculate closest Monday date of week containing January 1st
    if jan1_wday < 4
        origin = jan1 - Dates.Day(jan1_wday)
    else
        origin = jan1 + Dates.Day(7 - jan1_wday)
    end
    date = origin + Dates.Day((week - 1) * 7 + weekday)
    return date
end

function calculate_aggregation_cases(long_cases_df; week_aggregation = 1)
    subset_cases = subset(
        long_cases_df,
        :aggregation_weeks => x -> x .== week_aggregation,
    )
    return subset_cases[!, :cases],
    fill_aggregation_values(
        subset_cases[!, :cases]; week_aggregation = week_aggregation
    )
end

function fill_aggregation_values(vals; week_aggregation = 1)
    dimensions = ndims(vals)
    daily_aggregation = week_aggregation * 7
    inner = (daily_aggregation, repeat([1], dimensions - 1)...)
    return repeat(vals; inner = inner)
end
