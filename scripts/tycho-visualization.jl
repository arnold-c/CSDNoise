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
sim = 4

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
            # force = true,
            force = false,
        )
    end
end

#%%
lead_time_df = DataFrame(
    "noise_type" => String[],
    "noise_magnitude" => Float64[],
    "test_specification" => IndividualTestSpecification[],
    "week_aggregation" => Int64[],
    "ews_method" => EWSMethod[],
    "lead_time_dist" => Vector{Float64}[],
    "lead_time_median" => Float64[],
    "lead_time_lower" => Float64[],
    "lead_time_upper" => Float64[],
    "lead_time_units" => String[],
)

# push!(
#     lead_time_df,
#     (
#         noise_type = "poisson",
#         noise_magnitude = 1.0,
#         test_specification = IndividualTestSpecification(1.0, 1.0, 0),
#         week_aggregation = 1,
#         ews_method = Main.Backward,
#         lead_time = 0.0,
#         lead_time_units = "days",
#     ),
# )

res = 1

for noise_specification in (
    PoissonNoiseSpecification(1.0),
    PoissonNoiseSpecification(2.0),
    PoissonNoiseSpecification(3.0),
    PoissonNoiseSpecification(4.0),
    PoissonNoiseSpecification(5.0),
    PoissonNoiseSpecification(6.0),
    PoissonNoiseSpecification(7.0),
    PoissonNoiseSpecification(8.0),
)
    weekly_noise_arr = create_noise_arr(
        noise_specification, weekly_cases_arr; seed = 1234
    )[1]

    biweekly_noise_arr = create_noise_arr(
        noise_specification, biweekly_cases_arr;
    )[1]

    monthly_noise_arr = create_noise_arr(
        noise_specification, monthly_cases_arr;
    )[1]

    for (test_specification, ews_method) in
        Iterators.product(
        (
            IndividualTestSpecification(0.5, 0.5, 0),
            IndividualTestSpecification(0.8, 0.8, 0),
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(1.0, 0.0, 0),
        ),
        (Main.Backward,),
    )
        res = ews_lead_time_df!(
            lead_time_df,
            weekly_cases_arr,
            weekly_noise_arr,
            tycho_CA_measles_long_plotdata,
            test_specification;
            noise_type = get_noise_description(noise_specification),
            noise_magnitude = get_noise_magnitude(noise_specification),
            week_aggregation = 1,
            ews_method = ews_method,
            ews_aggregation = 1,
            ews_bandwidth = 52,
            ews_lag = 1,
            ews_metric = "variance",
            ews_threshold_window = Main.Expanding,
            ews_threshold_percentile = 0.95,
            consecutive_thresholds = 2,
            obsdate = cdc_week_to_date(1990, 3; weekday = 6),
            lead_time_units = :weeks,
            lead_time_percentile = 0.8,
            return_objects = true,
        )
    end
end

ews_lead_time_plot(
    lead_time_df;
    lead_time_units = :weeks,
)

#%%
hist(lead_time_df[1, :lead_time_dist]; bins = 140:1:350)

#%%
res_orig = deepcopy(res)

#%%
res[1][:, 5, 1] == res[1][:, 5, 100]

#%%
for sim in eachindex(res[2])
    @show res[2].variance[sim] == res_orig[2].variance[sim]
end

#%%
for sim in eachindex(res[3])
    @show res[3][sim] == res_orig[3][sim]
end

#%%
hcat(res[3][1], res_orig[3][1])

#%%
for sim in eachindex(res[2])
    @show res[2][sim].variance == res_orig[2][sim].variance
    @show res[2][1].variance == res_orig[2][sim].variance
end

#%%
for sim in eachindex(res[2])
    @assert expanding_ews_thresholds(
        res[2][sim],
        :variance,
        Main.Expanding;
        percentiles = 0.95,
    )[2] == expanding_ews_thresholds(
        res_orig[2][sim],
        :variance,
        Main.Expanding;
        percentiles = 0.95,
    )[2] "failed sim $sim"
end

#%%
hcat(
    expanding_ews_thresholds(
        res[2][1],
        :variance,
        Main.Expanding;
        percentiles = 0.95,
    )[1],
    expanding_ews_thresholds(
        res_orig[2][1],
        :variance,
        Main.Expanding;
        percentiles = 0.95,
    )[1],
)

#%%
expanding_ews_thresholds(
    res[2][1],
    :variance,
    Main.Expanding;
    percentiles = 0.95,
)[1] .==
expanding_ews_thresholds(
    res_orig[2][1],
    :variance,
    Main.Expanding;
    percentiles = 0.95,
)[1]

#%%
hcat(
    expanding_ews_thresholds(
        res[2][1],
        :variance,
        Main.Expanding;
        percentiles = 0.95,
    )...,
    expanding_ews_thresholds(
        res_orig[2][1],
        :variance,
        Main.Expanding;
        percentiles = 0.95,
    )...) |>
arr ->
    DataFrame(arr, [:res_var, :res_thresh, :orig_var, :orig_thresh]) |>
    df ->
        transform!(
            df, [:res_var, :orig_var] => ByRow((x1, x2) -> x1 - x2) => :var_diff
        ) |>
        df ->
            transform!(
                df,
                [:res_thresh, :orig_thresh] =>
                    ByRow((x1, x2) -> x1 - x2) => :thresh_diff,
            )

# df -> filter(row -> row.diff > 0, df)

#%%
@assert res[2][1].variance == res_orig[2][1].variance

#%%
for sim in eachindex(res[3])
    @assert calculate_ews_lead_time(
        res[3][sim];
        week_aggregation = 1,
        consecutive_thresholds = 2,
        output_type = :weeks,
    ) .== calculate_ews_lead_time(
        res_orig[3][sim];
        week_aggregation = 1,
        consecutive_thresholds = 2,
        output_type = :weeks,
    )
end

#%%
hcat(
    map(
        sim -> calculate_ews_lead_time(
            res[3][sim];
            week_aggregation = 1,
            consecutive_thresholds = 2,
            output_type = :weeks,
        ),
        eachindex(res[3]),
    ),
    map(
        sim -> calculate_ews_lead_time(
            res_orig[3][sim];
            week_aggregation = 1,
            consecutive_thresholds = 2,
            output_type = :weeks,
        ),
        eachindex(res_orig[3]),
    ),
)
