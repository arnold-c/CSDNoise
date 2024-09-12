#%%
using DrWatson
@quickactivate "CSDNoise"

using CSDNoise
using Dates
using CSV
using DataFrames
using GLMakie
using StatsBase
using StructArrays

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
base_plotdir = plotsdir("tycho")

incidence_epicurve = tycho_epicurve(
    plot_dates,
    weekly_plot_cases,
    biweekly_plot_cases,
    monthly_plot_cases;
    plottitle = "Incidence Epicurve",
)

incidence_epicurve

save(
    joinpath(base_plotdir, "incidence-epicurve.png"),
    incidence_epicurve;
    size = (2200, 1600),
)

#%%
nsims = 100

weekly_cases_arr = zeros(Int64, length(weekly_cases), 1, nsims)
weekly_cases_arr .= weekly_cases

biweekly_cases_arr = zeros(Int64, length(biweekly_cases), 1, nsims)
biweekly_cases_arr .= biweekly_cases

monthly_cases_arr = zeros(Int64, length(monthly_cases), 1, nsims)
monthly_cases_arr .= monthly_cases

#%%
sim = 1

for noise_specification in (
    PoissonNoiseSpecification(1.0),
    PoissonNoiseSpecification(8.0),
)
    weekly_noise_arr = create_noise_arr(noise_specification, weekly_cases_arr;)[1]
    filled_weekly_noise_arr = fill_aggregation_values(weekly_noise_arr)

    biweekly_noise_arr = create_noise_arr(
        noise_specification, biweekly_cases_arr;
    )[1]
    filled_biweekly_noise_arr = fill_aggregation_values(
        biweekly_noise_arr; week_aggregation = 2
    )

    monthly_noise_arr = create_noise_arr(
        noise_specification, monthly_cases_arr;
    )[1]
    filled_monthly_noise_arr =
        fill_aggregation_values(
            monthly_noise_arr; week_aggregation = 4
        ) |>
        x -> x[1:length(plot_dates), :]

    noise_epicurve = tycho_noise_components_epicurve(
        plot_dates,
        (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
        (
            filled_weekly_noise_arr[:, sim],
            filled_biweekly_noise_arr[:, sim],
            filled_monthly_noise_arr[:, sim],
        );
        plottitle = "Noise (100% Poisson) Epicurve",
    )

    noise_path = joinpath(
        base_plotdir,
        "testing-plots",
        getdirpath(noise_specification),
        "sim_$(sim)",
    )

    mkpath(noise_path)

    save(
        joinpath(noise_path, "noise-epicurve.png"),
        noise_epicurve;
        size = (2200, 1600),
    )

    for (test_specification, ews_method) in
        Iterators.product(
        (
            IndividualTestSpecification(0.5, 0.5, 0),
            IndividualTestSpecification(0.8, 0.8, 0),
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(1.0, 0.0, 0),
        ),
        (Main.Backward, Main.Centered),
    )
        tycho_testing_plots(
            (
                weekly_noise_arr,
                biweekly_noise_arr,
                monthly_noise_arr,
            ),
            (
                weekly_cases_arr,
                biweekly_cases_arr,
                monthly_cases_arr,
            ),
            (
                weekly_plot_cases,
                biweekly_plot_cases,
                monthly_plot_cases,
            ),
            tycho_CA_measles_long_plotdata;
            individual_test_specification = test_specification,
            noise_specification = noise_specification,
            ews_method = ews_method,
            sim = sim,
            plot_base_path = joinpath(base_plotdir, "testing-plots"),
            force = true,
        )
    end
end

#%%
hcat(
    1:316,
    weekly_100pc_rdt_8080_noise_800pc_backward_var_thresholds[2][:, 1],
    cumsum(weekly_100pc_rdt_8080_noise_800pc_backward_var_thresholds[2][:, 1]),
)

calculate_ews_lead_time(
    weekly_100pc_rdt_8080_noise_800pc_backward_var_thresholds[2][:, 1];
    consecutive_thresholds = 2,
    aggregation_weeks = 2,
    output_type = :years,
)
