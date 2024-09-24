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
using ProgressMeter

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
sims = (
    1,
    # 4
)
ews_metrics = [
    "autocorrelation",
    "autocovariance",
    "coefficient_of_variation",
    "index_of_dispersion",
    "kurtosis",
    "mean",
    "skewness",
    "variance",
]

for noise_specification in (
    PoissonNoiseSpecification(1.0),
    # PoissonNoiseSpecification(8.0),
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

    for sim in sims
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
    end

    # @showprogress for (test_specification, ews_method, sim) in
    #                   Iterators.product(
    #     (
    #         IndividualTestSpecification(0.5, 0.5, 0),
    #         IndividualTestSpecification(0.8, 0.8, 0),
    #         IndividualTestSpecification(1.0, 1.0, 0),
    #         IndividualTestSpecification(1.0, 0.0, 0),
    #     ),
    #     (
    #         Main.Backward,
    #         Main.Centered,
    #     ),
    #     sims,
    # )
    #     for ews_metric in ews_metrics
    #         tycho_testing_plots(
    #             (
    #                 weekly_noise_arr,
    #                 biweekly_noise_arr,
    #                 monthly_noise_arr,
    #             ),
    #             (
    #                 weekly_cases_arr,
    #                 biweekly_cases_arr,
    #                 monthly_cases_arr,
    #             ),
    #             (
    #                 weekly_plot_cases,
    #                 biweekly_plot_cases,
    #                 monthly_plot_cases,
    #             ),
    #             tycho_CA_measles_long_plotdata;
    #             individual_test_specification = test_specification,
    #             noise_specification = noise_specification,
    #             ews_metric = ews_metric,
    #             ews_method = ews_method,
    #             sim = sim,
    #             plot_base_path = joinpath(base_plotdir, "testing-plots"),
    #             force = false,
    #         )
    #     end
    # end

    for (ews_method, (cases_arr, noise_arr, week_aggregation), sim) in
        Iterators.product(
        (Main.Backward, Main.Centered),
        zip(
            (weekly_cases_arr, biweekly_cases_arr, monthly_cases_arr),
            (weekly_noise_arr, biweekly_noise_arr, monthly_noise_arr),
            (1, 2, 4),
        ),
        sims,
    )
        ews_df = tycho_tau_heatmap_df(
            tycho_CA_measles_long_plotdata,
            cases_arr,
            noise_arr,
            (
                IndividualTestSpecification(0.5, 0.5, 0),
                IndividualTestSpecification(0.8, 0.8, 0),
                IndividualTestSpecification(1.0, 1.0, 0),
            );
            week_aggregation = week_aggregation,
            ews_metrics = [
                "autocorrelation",
                "autocovariance",
                "coefficient_of_variation",
                "index_of_dispersion",
                "kurtosis",
                "mean",
                "skewness",
                "variance",
            ],
            ews_method = ews_method,
            statistic_function = StatsBase.median,
        )
    end
end

#%%
tycho_tau_heatmap_plot(
    ews_df
)

#%%
subset(
    ews_df,
    :test_specification => ByRow(==(IndividualTestSpecification(0.5, 0.5, 0))),
) |>
df ->
    sort(df, :ews_metric) |>
    df -> df[
        indexin(["mean", "variance"], df.ews_metric),
        [:ews_metric, :ews_metric_value],
    ]

#%%
lead_time_df = DataFrame(
    "noise_type" => String[],
    "noise_magnitude" => Float64[],
    "test_specification" => IndividualTestSpecification[],
    "week_aggregation" => Int64[],
    "ews_method" => EWSMethod[],
    "lead_time_dist" => Vector{Vector{Float64}}(),
    "lead_time_median" => Float64[],
    "lead_time_lower" => Float64[],
    "lead_time_upper" => Float64[],
    "lead_time_units" => String[],
)

res = 1

for (noise_specification, (week_aggregation, cases_arr)) in Iterators.product(
    (
        PoissonNoiseSpecification(1.0),
        PoissonNoiseSpecification(2.0),
        PoissonNoiseSpecification(3.0),
        PoissonNoiseSpecification(4.0),
        PoissonNoiseSpecification(5.0),
        PoissonNoiseSpecification(6.0),
        PoissonNoiseSpecification(7.0),
        PoissonNoiseSpecification(8.0),
    ),
    zip((1, 2, 4), (weekly_cases_arr, biweekly_cases_arr, monthly_cases_arr)),
)
    noise_arr = create_noise_arr(noise_specification, cases_arr)[1]

    for (test_specification, ews_method) in
        Iterators.product(
        (
            IndividualTestSpecification(0.5, 0.5, 0),
            IndividualTestSpecification(0.6, 0.6, 0),
            IndividualTestSpecification(0.7, 0.7, 0),
            IndividualTestSpecification(0.8, 0.8, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(1.0, 0.0, 0),
        ),
        (Main.Backward, Main.Centered),
    )
        res = ews_lead_time_df!(
            lead_time_df,
            cases_arr,
            noise_arr,
            tycho_CA_measles_long_plotdata,
            test_specification;
            noise_type = get_noise_description(noise_specification),
            noise_magnitude = get_noise_magnitude(noise_specification),
            week_aggregation = week_aggregation,
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
            lead_time_percentile = 0.5,
            return_objects = true,
        )
    end
end

#%%
lead_time_plotdir = joinpath(base_plotdir, "testing-plots", "lead-time-plots")
mkpath(lead_time_plotdir)

for (week_aggregation, ews_method) in
    Iterators.product((1, 2, 4), (Main.Backward, Main.Centered))
    lead_time_plot = ews_lead_time_plot(
        lead_time_df;
        week_aggregation = week_aggregation,
        lead_time_units = :weeks,
        ews_method = ews_method,
    )

    save(
        joinpath(
            lead_time_plotdir,
            "lead-time-plot_$(week_aggregation)-aggregation_$(method_string(ews_method)).png",
        ),
        lead_time_plot;
        size = (2200, 1600),
    )

    Makie.empty!(lead_time_plot)
end
