#%%
using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using Dates
using CSV
using DataFrames
using GLMakie
using StatsBase

include(srcdir("makie-plotting-setup.jl"))

#%%
bandwidth_weeks = 52
bandwidth = bandwidth_weeks * 7

obsdate = cdc_week_to_date(1990, 3; weekday = 6)
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
weekly_noise_100pc_arr = create_noise_arr(
    PoissonNoiseSpecification("poisson", 1.0),
    weekly_cases_arr;
)[1]
filled_weekly_noise_100pc_arr = fill_aggregation_values(weekly_noise_100pc_arr)

weekly_noise_800pc_arr = create_noise_arr(
    PoissonNoiseSpecification("poisson", 8.0),
    weekly_cases_arr;
)[1]
filled_weekly_noise_800pc_arr = fill_aggregation_values(weekly_noise_800pc_arr)

biweekly_noise_100pc_arr = create_noise_arr(
    PoissonNoiseSpecification("poisson", 1.0),
    biweekly_cases_arr;
)[1]
filled_biweekly_noise_100pc_arr = fill_aggregation_values(
    biweekly_noise_100pc_arr; week_aggregation = 2
)

biweekly_noise_800pc_arr = create_noise_arr(
    PoissonNoiseSpecification("poisson", 8.0),
    biweekly_cases_arr;
)[1]
filled_biweekly_noise_800pc_arr = fill_aggregation_values(
    biweekly_noise_800pc_arr; week_aggregation = 2
)

monthly_noise_100pc_arr = create_noise_arr(
    PoissonNoiseSpecification("poisson", 1.0),
    monthly_cases_arr;
)[1]
filled_monthly_noise_100pc_arr =
    fill_aggregation_values(monthly_noise_100pc_arr; week_aggregation = 4) |>
    x -> x[1:length(plot_dates), :]

monthly_noise_800pc_arr = create_noise_arr(
    PoissonNoiseSpecification("poisson", 8.0),
    monthly_cases_arr;
)[1]
filled_monthly_noise_800pc_arr =
    fill_aggregation_values(monthly_noise_800pc_arr; week_aggregation = 4) |>
    x -> x[1:length(plot_dates), :]

#%%
tycho_noise_components_epicurve(
    plot_dates,
    (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
    (
        filled_weekly_noise_100pc_arr[:, 1],
        filled_biweekly_noise_100pc_arr[:, 1],
        filled_monthly_noise_100pc_arr[:, 1],
    );
    plottitle = "Noise (100% Poisson) Epicurve",
)

#%%
weekly_100pc_perfect_test_noise_100pc_test_arr = create_testing_arrs(
    weekly_cases_arr,
    weekly_noise_100pc_arr,
    1.0,
    IndividualTestSpecification(1.0, 1.0, 0),
)

filled_weekly_100pc_perfect_test_noise_100pc_test_arr = fill_aggregation_values(
    weekly_100pc_perfect_test_noise_100pc_test_arr
)

biweekly_100pc_perfect_test_noise_100pc_test_arr = create_testing_arrs(
    biweekly_cases_arr,
    biweekly_noise_100pc_arr,
    1.0,
    IndividualTestSpecification(1.0, 1.0, 0),
)

filled_biweekly_100pc_perfect_test_noise_100pc_test_arr = fill_aggregation_values(
    biweekly_100pc_perfect_test_noise_100pc_test_arr; week_aggregation = 2
)

monthly_100pc_perfect_test_noise_100pc_test_arr = create_testing_arrs(
    monthly_cases_arr,
    monthly_noise_100pc_arr,
    1.0,
    IndividualTestSpecification(1.0, 1.0, 0),
)

filled_monthly_100pc_perfect_test_noise_100pc_test_arr =
    fill_aggregation_values(
        monthly_100pc_perfect_test_noise_100pc_test_arr; week_aggregation = 4
    ) |>
    x -> x[1:length(plot_dates), :, :]

#%%
tycho_epicurve(
    plot_dates,
    filled_weekly_100pc_perfect_test_noise_100pc_test_arr[:, 5, 1],
    filled_biweekly_100pc_perfect_test_noise_100pc_test_arr[:, 5, 1],
    filled_monthly_100pc_perfect_test_noise_100pc_test_arr[:, 5, 1];
    plottitle = "Test Positives",
    subtitle = "Perfect Tests (100% Testing), Noise (100% Poisson) Epicurve",
)

#%%
weekly_100pc_rdt_8080_noise_100pc_test_arr = create_testing_arrs(
    weekly_cases_arr,
    weekly_noise_100pc_arr,
    1.0,
    IndividualTestSpecification(0.8, 0.8, 0),
)

filled_weekly_100pc_rdt_8080_noise_100pc_test_arr = fill_aggregation_values(
    weekly_100pc_rdt_8080_noise_100pc_test_arr
)

weekly_100pc_rdt_8080_noise_800pc_test_arr = create_testing_arrs(
    weekly_cases_arr,
    weekly_noise_800pc_arr,
    1.0,
    IndividualTestSpecification(0.8, 0.8, 0),
)

filled_weekly_100pc_rdt_8080_noise_800pc_test_arr = fill_aggregation_values(
    weekly_100pc_rdt_8080_noise_800pc_test_arr
)

biweekly_100pc_rdt_8080_noise_100pc_test_arr = create_testing_arrs(
    biweekly_cases_arr,
    biweekly_noise_100pc_arr,
    1.0,
    IndividualTestSpecification(0.8, 0.8, 0),
)

filled_biweekly_100pc_rdt_8080_noise_100pc_test_arr = fill_aggregation_values(
    biweekly_100pc_rdt_8080_noise_100pc_test_arr; week_aggregation = 2
)

biweekly_100pc_rdt_8080_noise_800pc_test_arr = create_testing_arrs(
    biweekly_cases_arr,
    biweekly_noise_800pc_arr,
    1.0,
    IndividualTestSpecification(0.8, 0.8, 0),
)

filled_biweekly_100pc_rdt_8080_noise_800pc_test_arr = fill_aggregation_values(
    biweekly_100pc_rdt_8080_noise_800pc_test_arr; week_aggregation = 2
)

monthly_100pc_rdt_8080_noise_100pc_test_arr = create_testing_arrs(
    monthly_cases_arr,
    monthly_noise_100pc_arr,
    1.0,
    IndividualTestSpecification(0.8, 0.8, 0),
)

filled_monthly_100pc_rdt_8080_noise_100pc_test_arr =
    fill_aggregation_values(
        monthly_100pc_rdt_8080_noise_100pc_test_arr; week_aggregation = 4
    ) |>
    x -> x[1:length(plot_dates), :, :]

monthly_100pc_rdt_8080_noise_800pc_test_arr = create_testing_arrs(
    monthly_cases_arr,
    monthly_noise_800pc_arr,
    1.0,
    IndividualTestSpecification(0.8, 0.8, 0),
)

filled_monthly_100pc_rdt_8080_noise_800pc_test_arr =
    fill_aggregation_values(
        monthly_100pc_rdt_8080_noise_800pc_test_arr; week_aggregation = 4
    ) |>
    x -> x[1:length(plot_dates), :, :]

#%%
tycho_epicurve(
    plot_dates,
    filled_weekly_100pc_rdt_8080_noise_100pc_test_arr[:, 5, 1],
    filled_biweekly_100pc_rdt_8080_noise_100pc_test_arr[:, 5, 1],
    filled_monthly_100pc_rdt_8080_noise_100pc_test_arr[:, 5, 1];
    plottitle = "Test Positives",
    subtitle = "RDT 80/80 (100% Testing), Noise (100% Poisson) Epicurve",
)

#%%
tycho_epicurve(
    plot_dates,
    filled_weekly_100pc_rdt_8080_noise_800pc_test_arr[:, 5, 1],
    filled_biweekly_100pc_rdt_8080_noise_800pc_test_arr[:, 5, 1],
    filled_monthly_100pc_rdt_8080_noise_800pc_test_arr[:, 5, 1];
    plottitle = "Test Positives",
    subtitle = "RDT 80/80 (100% Testing), Noise (800% Poisson) Epicurve",
)

#%%
tycho_test_positive_components_epicurve(
    plot_dates,
    (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
    (
        filled_weekly_100pc_rdt_8080_noise_100pc_test_arr[:, :, 1],
        filled_biweekly_100pc_rdt_8080_noise_100pc_test_arr[:, :, 1],
        filled_monthly_100pc_rdt_8080_noise_100pc_test_arr[:, :, 1],
    );
    plottitle = "RDT 80/80 (100% Testing), Noise (100% Poisson) Epicurves",
)

#%%
tycho_test_positive_components_epicurve(
    plot_dates,
    (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
    (
        filled_weekly_100pc_rdt_8080_noise_800pc_test_arr[:, :, 1],
        filled_biweekly_100pc_rdt_8080_noise_800pc_test_arr[:, :, 1],
        filled_monthly_100pc_rdt_8080_noise_800pc_test_arr[:, :, 1],
    );
    plottitle = "RDT 80/80 (100% Testing), Noise (800% Poisson) Epicurves",
)

#%%
weekly_obs_index = calculate_ews_enddate(
    tycho_CA_measles_long_plotdata;
    week_aggregation = 1,
    obsdate = obsdate,
)

biweekly_obs_index = calculate_ews_enddate(
    tycho_CA_measles_long_plotdata;
    week_aggregation = 2,
    obsdate = obsdate,
)

monthly_obs_index = calculate_ews_enddate(
    tycho_CA_measles_long_plotdata;
    week_aggregation = 4,
    obsdate = obsdate,
)

weekly_cases_ewsmetrics = EWSMetrics(
    EWSMetricSpecification(
        Centered,
        1,
        52,
        1,
    ),
    weekly_cases[1:weekly_obs_index],
)

biweekly_cases_ewsmetrics = EWSMetrics(
    EWSMetricSpecification(
        Centered,
        1,
        Int(52 / 2),
        1,
    ),
    biweekly_cases[1:biweekly_obs_index],
)

monthly_cases_ewsmetrics = EWSMetrics(
    EWSMetricSpecification(
        Centered,
        1,
        Int(52 / 4),
        1,
    ),
    monthly_cases[1:monthly_obs_index],
)

#%%
# Compare against R version
tycho_spaero_cases_ews_df = CSV.read(
    projectdir("out", "cases_ews_df.csv"), DataFrame
)
tycho_spaero_cases_ews_tau_df = CSV.read(
    projectdir("out", "cases_ews_tau_df.csv"), DataFrame
)

isapprox(
    tycho_spaero_cases_ews_df.variance_1wk,
    weekly_cases_ewsmetrics.variance;
    atol = 1e-5,
)

isapprox(
    parse.(
        Float64,
        String.(
            filter(x -> !(x .== "NA"), tycho_spaero_cases_ews_df.variance_4wk)
        ),
    ),
    monthly_cases_ewsmetrics.variance;
    atol = 1e-5,
)

isapprox(
    subset(tycho_spaero_cases_ews_tau_df, :aggregation => x -> x .== "1wk") |>
    df -> subset(df, :statistic => x -> x .== "variance")[1, :value],
    weekly_cases_ewsmetrics.variance_tau,
)

isapprox(
    subset(tycho_spaero_cases_ews_tau_df, :aggregation => x -> x .== "4wk") |>
    df -> subset(df, :statistic => x -> x .== "variance")[1, :value],
    monthly_cases_ewsmetrics.variance_tau,
)
