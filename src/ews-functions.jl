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
    consecutive_thresholds = 2,
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
    sim = 1,
    return_objects = false,
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
    noise_magnitude_description = get_noise_magnitude_description(
        noise_specification
    )

    test_path_description = "sens-$(individual_test_specification.sensitivity)_spec-$(individual_test_specification.specificity)_lag-$(individual_test_specification.test_result_lag)"
    noise_path_description = "$(
        replace(noise_magnitude_description," " => "_", ": " => "_")
    )"

    noise_path = joinpath(
        plot_base_path,
        getdirpath(noise_specification),
        "sim_$(sim)",
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

    weekly_detection_index = calculate_ews_trigger_index(
        weekly_thresholds[2];
        consecutive_thresholds = consecutive_thresholds,
    )

    biweekly_detection_index = calculate_ews_trigger_index(
        biweekly_thresholds[2];
        consecutive_thresholds = consecutive_thresholds,
    )

    monthly_detection_index = calculate_ews_trigger_index(
        monthly_thresholds[2];
        consecutive_thresholds = consecutive_thresholds,
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
            ),
            (
                weekly_detection_index,
                biweekly_detection_index,
                monthly_detection_index,
            );
            plottitle = "Test Positives",
            subtitle = "$(test_plot_description), $(noise_magnitude_description): $(method_string(weekly_ewsmetric_specification.method)) EWS $(ews_metric) Epicurve",
            ews_ylabel = ews_metric,
            threshold_percentile = ews_threshold_percentile,
            consecutive_thresholds = consecutive_thresholds,
        )

        save(
            plotpath,
            plot;
            size = (2200, 1600),
        )

        Makie.empty!(plot)
    end

    if return_objects
        return (
            plot_dates,
            weekly_test_arr,
            filled_weekly_test_arr,
            biweekly_test_arr,
            filled_biweekly_test_arr,
            monthly_test_arr,
            filled_monthly_test_arr,
            weekly_obs_index,
            biweekly_obs_index,
            monthly_obs_index,
            weekly_cases_ewsmetrics,
            biweekly_cases_ewsmetrics,
            monthly_cases_ewsmetrics,
            weekly_test_ewsmetrics,
            biweekly_test_ewsmetrics,
            monthly_test_ewsmetrics,
            weekly_thresholds,
            biweekly_thresholds,
            monthly_thresholds,
        )
    end

    return nothing
end

function calculate_ews_lead_time(
    ews_thresholds;
    week_aggregation = 1,
    consecutive_thresholds = 2,
    output_type = :days,
)
    threshold_index = calculate_ews_trigger_index(
        ews_thresholds; consecutive_thresholds = consecutive_thresholds
    )

    return calculate_ews_lead_time(
        ews_thresholds,
        threshold_index;
        week_aggregation = week_aggregation,
        output_type = output_type,
    )
end

function calculate_ews_lead_time(
    ews_thresholds, threshold_index;
    week_aggregation = 1,
    output_type = :days,
)
    output_multiplier = @match output_type begin
        :days => 7
        :weeks => 1
        :months => 7 / 30.5
        :years => 7 / 365
        _ =>
            error(
                "Unknown output type: $(output_type).\nChoose between :days, :weeks, :months, or :years."
            )
    end

    if isnothing(threshold_index)
        @warn "No ews trigger index found. Returning nothing."
        return nothing
    end

    return (length(ews_thresholds) - threshold_index) * week_aggregation *
           output_multiplier
end

function calculate_ews_trigger_index(
    ews_thresholds::T1;
    consecutive_thresholds = 2,
) where {T1<:AbstractMatrix{<:Bool}}
    reshaped_ews_thresholds = reshape(ews_thresholds, :)

    @assert length(reshaped_ews_thresholds) == length(ews_thresholds[:, 1])

    return calculate_ews_trigger_index(
        reshaped_ews_thresholds;
        consecutive_thresholds = consecutive_thresholds,
    )
end

function calculate_ews_trigger_index(
    ews_thresholds::T1;
    consecutive_thresholds = 2,
) where {T1<:AbstractVector{<:Bool}}
    cumulative_thresholds = cumsum(ews_thresholds)
    for (i, v) in pairs(cumulative_thresholds)
        if i > consecutive_thresholds &&
            v >=
           cumulative_thresholds[i - consecutive_thresholds] +
           consecutive_thresholds
            return i
        end
    end
    return nothing
end

function ews_lead_time_df!(
    lead_time_df::DataFrame,
    cases_arr,
    noise_arr,
    long_plotdata,
    individual_test_specification;
    noise_type = "poisson",
    noise_magnitude = 1.0,
    week_aggregation = 1,
    ews_method = Main.Centered,
    ews_aggregation = 1,
    ews_bandwidth = 52,
    ews_lag = 1,
    ews_metric = "variance",
    ews_threshold_window = Main.Expanding,
    ews_threshold_percentile = 0.95,
    consecutive_thresholds = 2,
    obsdate = cdc_week_to_date(1990, 3; weekday = 6),
    lead_time_units = :days,
    lead_time_percentile = 0.95,
)
    test_arr = create_testing_arrs(
        cases_arr,
        noise_arr,
        1.0,
        individual_test_specification,
    )

    obs_index = calculate_ews_enddate(
        long_plotdata;
        week_aggregation = week_aggregation,
        obsdate = obsdate,
    )

    test_ewsmetrics =
        map(
            k -> EWSMetrics(
                EWSMetricSpecification(
                    ews_method,
                    ews_aggregation,
                    ews_bandwidth,
                    ews_lag,
                ),
                test_arr[
                    1:obs_index, 5, k
                ],
            ),
            axes(test_arr, 3),
        ) |>
        x -> StructArray(x)

    ews_metric_sym = Symbol(ews_metric)

    ews_lead_time = zeros(Float64, length(test_ewsmetrics))
    for sim in axes(test_ewsmetrics, 1)
        weekly_thresholds = expanding_ews_thresholds(
            test_ewsmetrics[sim],
            ews_metric_sym,
            ews_threshold_window;
            percentiles = ews_threshold_percentile,
        )[2]

        ews_lead_time[sim] = calculate_ews_lead_time(
            weekly_thresholds;
            week_aggregation = week_aggregation,
            consecutive_thresholds = consecutive_thresholds,
            output_type = lead_time_units,
        )
    end

    percentile_tail = (1 - lead_time_percentile) / 2

    push!(
        lead_time_df,
        (
            noise_type = noise_type,
            noise_magnitude = noise_magnitude,
            test_specification = individual_test_specification,
            week_aggregation = 1,
            ews_method = ews_method,
            lead_time_dist = ews_lead_time,
            lead_time_median = StatsBase.median(ews_lead_time),
            lead_time_lower = StatsBase.quantile(
                ews_lead_time, percentile_tail
            ),
            lead_time_upper = StatsBase.quantile(
                ews_lead_time, 1.0 - percentile_tail
            ),
            lead_time_units = String(lead_time_units),
        ),
    )

    return nothing
end
