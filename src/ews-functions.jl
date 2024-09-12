using SumTypes
using StatsBase: StatsBase
using StructArrays

@sum_type EWSThresholdWindow begin
    Expanding
    Rolling
end

function expanding_ews_thresholds(
    ewsmetrics::T1,
    metric::T2,
    window_type::T3;
    percentiles = (0.8, 0.95),
    burn_in = 10,
) where {T1<:EWSMetrics,T2<:Symbol,T3<:EWSThresholdWindow}
    ews_vec = getproperty(ewsmetrics, metric)

    ews_distributions = fill(NaN, length(ews_vec), length(percentiles))
    exceeds_thresholds = Array{Bool}(
        undef, length(ews_vec), length(percentiles)
    )

    for i in eachindex(ews_vec)
        if i <= burn_in
            continue
        end
        ews_distributions[i, :] .= map(
            p -> StatsBase.quantile(@view(ews_vec[1:i]), p),
            percentiles,
        )

        exceeds_thresholds[i, :] .= ews_vec[i] .>= ews_distributions[i, :]
    end

    return ews_distributions, exceeds_thresholds
end

function tycho_testing_plots(
    noise_arr_tuple,
    cases_arr_tuple,
    plot_cases_tuple,
    long_plotdata;
    individual_test_specification = IndividualTestSpecification(1.0, 1.0, 0),
    noise_specification = PoissonNoiseSpecification("poisson", 1.0),
    plot_base_path = DrWatson.plotsdir("tycho/testing-plots"),
    ews_method = Main.Centered,
    ews_aggregation = 1,
    ews_bandwidth = 52,
    ews_lag = 1,
    ews_metric = "variance",
    ews_threshold_window = Main.Expanding,
    ews_threshold_percentile = 0.95,
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
    sim = 1,
    force = true,
)
    @assert length(noise_arr_tuple) == length(cases_arr_tuple) == 3
    weekly_noise_arr, biweekly_noise_arr, monthly_noise_arr = noise_arr_tuple
    weekly_cases_arr, biweekly_cases_arr, monthly_cases_arr = cases_arr_tuple
    weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases =
        plot_cases_tuple

    @assert length(weekly_noise_arr) == length(weekly_cases_arr)
    @assert length(biweekly_noise_arr) == length(biweekly_cases_arr)
    @assert length(monthly_noise_arr) == length(monthly_cases_arr)
    @assert length(weekly_noise_arr) > length(biweekly_noise_arr) >
        length(monthly_noise_arr)

    plot_dates =
        unique(long_plotdata[!, :date]) |>
        x -> collect(minimum(x):Day(1):(maximum(x) + Day(6)))

    weekly_test_arr = create_testing_arrs(
        weekly_cases_arr,
        weekly_noise_arr,
        1.0,
        individual_test_specification,
    )

    filled_weekly_test_arr = fill_aggregation_values(
        weekly_test_arr
    )

    biweekly_test_arr = create_testing_arrs(
        biweekly_cases_arr,
        biweekly_noise_arr,
        1.0,
        individual_test_specification,
    )

    filled_biweekly_test_arr = fill_aggregation_values(
        biweekly_test_arr; week_aggregation = 2
    )

    monthly_test_arr = create_testing_arrs(
        monthly_cases_arr,
        monthly_noise_arr,
        1.0,
        individual_test_specification,
    )

    filled_monthly_test_arr =
        fill_aggregation_values(
            monthly_test_arr;
            week_aggregation = 4,
        ) |>
        x -> x[1:length(plot_dates), :, :]

    test_plot_description = get_test_description(individual_test_specification)
    noise_magnitude_description = get_noise_magnitude(noise_specification)

    test_path_description = "sens-$(individual_test_specification.sensitivity)_spec-$(individual_test_specification.specificity)_lag-$(individual_test_specification.test_result_lag)"
    noise_path_description = "$(
        replace(noise_magnitude_description," " => "_", ": " => "_")
    )"

    noise_path = joinpath(
        plot_base_path, getdirpath(noise_specification)
    )
    test_dir = joinpath(noise_path, test_path_description)
    mkpath(test_dir)

    plotname = "testing-epicurve_$(test_path_description)_$(noise_path_description).png"
    plotpath = joinpath(test_dir, plotname)

    if !isfile(plotpath) || force
        plot = tycho_epicurve(
            plot_dates,
            filled_weekly_test_arr[:, 5, sim],
            filled_biweekly_test_arr[:, 5, sim],
            filled_monthly_test_arr[:, 5, sim];
            plottitle = "Test Positives",
            subtitle = "$(test_plot_description), $(noise_magnitude_description) Epicurve",
        )

        save(
            plotpath,
            plot;
            size = (2200, 1600),
        )

        Makie.empty!(plot)
    end

    plotname = "testing-components-epicurve_$(test_path_description)_$(noise_path_description).png"
    plotpath = joinpath(test_dir, plotname)

    if !isfile(plotpath) || force
        plot = tycho_test_positive_components_epicurve(
            plot_dates,
            (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
            (
                filled_weekly_test_arr[:, :, sim],
                filled_biweekly_test_arr[:, :, sim],
                filled_monthly_test_arr[:, :, sim],
            );
            plottitle = "$(test_plot_description), $(noise_magnitude_description) Component Epicurves",
        )

        save(
            plotpath,
            plot;
            size = (2200, 1600),
        )

        Makie.empty!(plot)
    end

    weekly_obs_index = calculate_ews_enddate(
        long_plotdata;
        week_aggregation = 1,
        obsdate = obsdate,
    )

    biweekly_obs_index = calculate_ews_enddate(
        long_plotdata;
        week_aggregation = 2,
        obsdate = obsdate,
    )

    monthly_obs_index = calculate_ews_enddate(
        long_plotdata;
        week_aggregation = 4,
        obsdate = obsdate,
    )

    weekly_ewsmetric_specification = EWSMetricSpecification(
        ews_method,
        ews_aggregation,
        ews_bandwidth,
        ews_lag,
    )

    weekly_cases_ewsmetrics = EWSMetrics(
        weekly_ewsmetric_specification,
        weekly_cases_arr[1:weekly_obs_index, :, sim],
    )

    biweekly_cases_ewsmetrics = EWSMetrics(
        EWSMetricSpecification(
            ews_method,
            ews_aggregation,
            Int(ews_bandwidth / 2),
            ews_lag,
        ),
        biweekly_cases_arr[1:biweekly_obs_index, :, sim],
    )

    monthly_cases_ewsmetrics = EWSMetrics(
        EWSMetricSpecification(
            ews_method,
            ews_aggregation,
            Int(ews_bandwidth / 4),
            ews_lag,
        ),
        monthly_cases_arr[1:monthly_obs_index, :, sim],
    )

    weekly_test_ewsmetrics =
        map(
            k -> EWSMetrics(
                EWSMetricSpecification(
                    ews_method,
                    ews_aggregation,
                    ews_bandwidth,
                    ews_lag,
                ),
                weekly_test_arr[
                    1:weekly_obs_index, 5, k
                ],
            ),
            axes(weekly_test_arr, 3),
        ) |>
        x -> StructArray(x)

    biweekly_test_ewsmetrics =
        map(
            k -> EWSMetrics(
                EWSMetricSpecification(
                    ews_method,
                    ews_aggregation,
                    Int(ews_bandwidth / 2),
                    ews_lag,
                ),
                biweekly_test_arr[
                    1:biweekly_obs_index, 5, k
                ],
            ),
            axes(biweekly_test_arr, 3),
        ) |>
        x -> StructArray(x)

    monthly_test_ewsmetrics =
        map(
            k -> EWSMetrics(
                EWSMetricSpecification(
                    ews_method,
                    ews_aggregation,
                    Int(ews_bandwidth / 4),
                    ews_lag,
                ),
                monthly_test_arr[
                    1:monthly_obs_index, 5, k
                ],
            ),
            axes(monthly_test_arr, 3),
        ) |>
        x -> StructArray(x)

    ews_metric_tau = ews_metric * "_tau"
    ews_metric_sym = Symbol(ews_metric)
    ews_metric_tau_sym = Symbol(ews_metric_tau)

    plotname = "$(ews_metric_tau)_distribution_$(test_path_description)_$(noise_path_description).png"
    ewspath = joinpath(test_dir, weekly_ewsmetric_specification.dirpath)
    mkpath(ewspath)
    plotpath = joinpath(
        ewspath, plotname
    )

    if !isfile(plotpath) || force
        plot = tycho_tau_distribution(
            (
                weekly = weekly_test_ewsmetrics,
                biweekly = biweekly_test_ewsmetrics,
                monthly = monthly_test_ewsmetrics,
            ),
            weekly_cases_ewsmetrics,
            ews_metric_tau_sym;
            plottitle = "$(test_plot_description), $(noise_magnitude_description): $(method_string(weekly_ewsmetric_specification.method)) EWS $(ews_metric_tau) Distribution",
        )

        save(
            plotpath,
            plot;
            size = (2200, 1600),
        )

        Makie.empty!(plot)
    end

    weekly_thresholds = expanding_ews_thresholds(
        weekly_test_ewsmetrics[sim],
        ews_metric_sym,
        ews_threshold_window;
        percentiles = ews_threshold_percentile,
    )

    biweekly_thresholds = expanding_ews_thresholds(
        biweekly_test_ewsmetrics[sim],
        ews_metric_sym,
        ews_threshold_window;
        percentiles = ews_threshold_percentile,
    )

    monthly_thresholds = expanding_ews_thresholds(
        monthly_test_ewsmetrics[sim],
        ews_metric_sym,
        ews_threshold_window;
        percentiles = ews_threshold_percentile,
    )

    plotname = "$(ews_metric)_$(100*ews_threshold_percentile)-percentile_thresholds_$(test_path_description)_$(noise_path_description).png"
    plotpath = joinpath(
        ewspath, plotname
    )

    if !isfile(plotpath) || force
        plot = tycho_epicurve(
            plot_dates,
            (
                filled_weekly_test_arr[:, 5, sim],
                filled_biweekly_test_arr[:, 5, sim],
                filled_monthly_test_arr[:, 5, sim],
            ),
            (
                getproperty(weekly_test_ewsmetrics[sim], ews_metric_sym),
                getproperty(biweekly_test_ewsmetrics[sim], ews_metric_sym),
                getproperty(monthly_test_ewsmetrics[sim], ews_metric_sym),
            ),
            (
                weekly_thresholds[2],
                biweekly_thresholds[2],
                monthly_thresholds[2],
            );
            plottitle = "Test Positives",
            subtitle = "$(test_plot_description), $(noise_magnitude_description): $(method_string(weekly_ewsmetric_specification.method)) EWS $(ews_metric) Epicurve",
            ews_ylabel = ews_metric,
        )

        save(
            plotpath,
            plot;
            size = (2200, 1600),
        )

        Makie.empty!(plot)
    end

    return nothing
end
