#%%
using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using Dates
using CSV
using DataFrames
using GLMakie

include(srcdir("makie-plotting-setup.jl"))

#%%
bandwidth_weeks = 52
bandwidth = bandwidth_weeks * 7

obsdate = cdc_week_to_date(1990, 34; weekday = 6)
plotxmin_year = 1984
plotxmin = cdc_week_to_date(plotxmin_year, 1; weekday = 6)
plotxmax_year = 1992
plotxmax = cdc_week_to_date(plotxmax_year, 1; weekday = 6) + Day(bandwidth)

#%%
tycho_CA_measles_long_plotdata = CSV.read(
    outdir("tycho_CA_measles_long_plotdata.csv"), DataFrame
)

tycho_CA_measles_wide_plotdata = CSV.read(
    outdir("tycho_CA_measles_wide_plotdata.csv"), DataFrame
)

plot_dates =
    unique(tycho_CA_measles_long_plotdata[!, :date]) |>
    x -> collect(minimum(x):Day(1):(maximum(x) + Day(6)))

#%%
subset(tycho_CA_measles_long_plotdata, :aggregation_weeks => x -> x .== 1) |>
df -> diff(df[!, :date]) |>
      x -> unique(x) == [Day(7)]

#%%
weekly_cases, weekly_plot_cases = calculate_aggregation_cases(
    tycho_CA_measles_long_plotdata; week_aggregation = 1
)

biweekly_cases, biweekly_plot_cases = calculate_aggregation_cases(
    tycho_CA_measles_long_plotdata; week_aggregation = 2
)

monthly_cases, monthly_plot_cases = calculate_aggregation_cases(
    tycho_CA_measles_long_plotdata; week_aggregation = 4
)
monthly_plot_cases = monthly_plot_cases[1:length(plot_dates)]

#%%
incidence_epicurve = tycho_epicurve(
    plot_dates,
    weekly_plot_cases,
    biweekly_plot_cases,
    monthly_plot_cases;
    plottitle = "Incidence Epicurve",
)

incidence_epicurve

#%%
nsims = 100

weekly_cases_arr = zeros(Int64, length(weekly_cases), 1, nsims)
weekly_cases_arr .= weekly_cases

biweekly_cases_arr = zeros(Int64, length(biweekly_cases), 1, nsims)
biweekly_cases_arr .= biweekly_cases

monthly_cases_arr = zeros(Int64, length(monthly_cases), 1, nsims)
monthly_cases_arr .= monthly_cases

#%%
weekly_noise_100pc_arr =
    create_noise_arr(
        PoissonNoiseSpecification("poisson", 1.0),
        weekly_cases_arr;
    )[1] |>
    x -> fill_aggregation_values(x)

weekly_noise_800pc_arr =
    create_noise_arr(
        PoissonNoiseSpecification("poisson", 8.0),
        weekly_cases_arr;
    )[1] |>
    x -> fill_aggregation_values(x)

biweekly_noise_100pc_arr =
    create_noise_arr(
        PoissonNoiseSpecification("poisson", 1.0),
        biweekly_cases_arr;
    )[1] |>
    x -> fill_aggregation_values(x; week_aggregation = 2)

biweekly_noise_800pc_arr =
    create_noise_arr(
        PoissonNoiseSpecification("poisson", 8.0),
        biweekly_cases_arr;
    )[1] |>
    x -> fill_aggregation_values(x; week_aggregation = 2)

monthly_noise_100pc_arr =
    create_noise_arr(
        PoissonNoiseSpecification("poisson", 1.0),
        monthly_cases_arr;
    )[1] |>
    x ->
        fill_aggregation_values(x; week_aggregation = 4) |>
        x -> x[1:length(plot_dates)]

monthly_noise_800pc_arr =
    create_noise_arr(
        PoissonNoiseSpecification("poisson", 8.0),
        monthly_cases_arr;
    )[1] |>
    x ->
        fill_aggregation_values(x; week_aggregation = 4) |>
        x -> x[1:length(plot_dates)]

#%%
tycho_noise_epicurve(
    plot_dates,
    (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
    (weekly_noise_100pc_arr[:, 1],
        biweekly_noise_100pc_arr[:, 1],
        monthly_noise_100pc_arr[:, 1]);
    plottitle = "Noise (100% Poisson) Epicurve",
)
