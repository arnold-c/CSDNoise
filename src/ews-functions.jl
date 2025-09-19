using LightSumTypes: variant
using StatsBase: StatsBase
using StructArrays
using DataFrames
using Try: Try
using TryExperimental: trygetindex
using Dates
using UnPack: @unpack

function calculate_bandwidth_and_return_ews_metric_spec(
        ews_metric_spec_components...
    )
    @assert length(ews_metric_spec_components) == 4

    bandwidth = calculate_bandwidth(
        ews_metric_spec_components[3],
        ews_metric_spec_components[2],
    )

    return EWSMetricSpecification(
        ews_metric_spec_components[1],
        ews_metric_spec_components[2],
        bandwidth,
        ews_metric_spec_components[4],
    )
end

function calculate_bandwidth(bandwidth_days, aggregation_days)
    return bandwidth_days ÷ aggregation_days
end

function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::EWSThresholdWindowType,
        percentile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    return exceeds_ews_threshold(
        ewsmetrics,
        metric,
        variant(window_type),
        percentile,
        burn_in,
    )
end

# TODO: Implement this method. Place holder as needed for type stability
function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::RollingThresholdWindow,
        percentile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    return fill(false, 2)
end

function exceeds_ews_threshold(
        ewsmetrics::T1,
        metric::T2,
        window_type::ExpandingThresholdWindow,
        percentile::Float64 = 0.95,
        burn_in::P = Dates.Day(10),
    ) where {T1 <: EWSMetrics, T2 <: Symbol, P <: Dates.Period}
    ews_vec = get_ews_metric_vec(ewsmetrics, metric)

    @unpack aggregation = ewsmetrics.ews_specification
    burn_in_index = Int64(Dates.days(burn_in) ÷ Dates.days(aggregation))

    @assert burn_in_index >= 1 && burn_in_index <= length(ews_vec)

    return _expanding_ews_thresholds(
        ews_vec,
        percentile,
        burn_in_index,
    )
end

function get_ews_metric_vec(
        ewsmetrics::T1,
        metric::T2,
    ) where {T1 <: EWSMetrics, T2 <: Symbol}
    @assert metric in [
        :mean,
        :variance,
        :coefficient_of_variation,
        :index_of_dispersion,
        :skewness,
        :kurtosis,
        :autocovariance,
        :autocorrelation,
    ]

    ews_vec = getproperty(ewsmetrics, metric)::Vector{Float64}
    return ews_vec
end


function _expanding_ews_thresholds(
        ews_vec::Vector{Float64},
        percentile::Float64,
        burn_in_index::Int64,
    )
    # ) where {T1 <: Integer, F <: AbstractFloat}
    ews_vec_len = length(ews_vec)
    ews_distributions = fill(NaN, ews_vec_len)
    ews_worker_vec = fill(
        NaN, sum((!isnan).(ews_vec))
    )
    exceeds_thresholds = zeros(Bool, ews_vec_len)

    worker_ind = 0
    for i in eachindex(ews_vec)
        if isnan(ews_vec[i])
            if i > burn_in_index
                ews_distributions[i] = ews_distributions[(i - 1)]
            end
            continue
        end
        worker_ind += 1
        # use online stats to build up new distribution to avoid computing quantiles for vectors containing NaNs
        ews_worker_vec[worker_ind] = ews_vec[i]
        if i > burn_in_index
            ews_distributions[i] = StatsBase.quantile(
                @view(ews_worker_vec[1:worker_ind]),
                percentile
            )

            exceeds_thresholds[i] = ews_vec[i] >= ews_distributions[(i - 1)]
        end
    end

    @assert worker_ind == length(ews_worker_vec)

    return exceeds_thresholds
end

function calculate_ews_enddate(
        thresholds::T,
        enddate_type::EWSEndDateType,
    ) where {T <: AbstractArray{<:Int}}

    enddate = _calculate_ews_enddate(thresholds, enddate_type)

    if Try.isok(enddate)
        return enddate
    end

    return Try.Err(
        "Failed to calculate ews_enddate for $enddate_type"
    )
end

_calculate_ews_enddate(thresholds::T, enddate_type::EWSEndDateType) where {T <: AbstractArray{<:Int}} = _calculate_ews_enddate(thresholds, variant(enddate_type))

_calculate_ews_enddate(thresholds, enddate_type::Reff_start) = trygetindex(thresholds, 1, 1)
_calculate_ews_enddate(thresholds, enddate_type::Reff_end) = trygetindex(thresholds, 1, 2)

function _calculate_ews_enddate(thresholds, enddate_type::Outbreak_start)
    filtered_outbreak_thresholds = filter_outbreak_thresholds(thresholds)
    return trygetindex(filtered_outbreak_thresholds, 1, 1)
end

function _calculate_ews_enddate(thresholds, enddate_type::Outbreak_end)
    filtered_outbreak_thresholds = filter_outbreak_thresholds(thresholds)

    return trygetindex(filtered_outbreak_thresholds, 1, 2)
end

function _calculate_ews_enddate(thresholds, enddate_type::Outbreak_middle)
    filtered_outbreak_thresholds = filter_outbreak_thresholds(thresholds)

    if size(filtered_outbreak_thresholds, 1) == 0
        return Try.Err(BoundsError(thresholds, (1, 1)))
    end

    return Try.Ok(
        filtered_outbreak_thresholds[1, 1] +
            (
            filtered_outbreak_thresholds[1, 2] -
                filtered_outbreak_thresholds[1, 1] +
                1
        ) ÷ 2,
    )
end

# function calculate_Reff_start_ews_enddate(thresholds::T)::Union{Try.Ok{Int}, Try.Err} where {T <: AbstractArray{<:Int}}
#     return trygetindex(thresholds, 1, 1)
# end
#
# function calculate_Reff_end_ews_enddate(thresholds::T)::Union{Try.Ok{Int}, Try.Err} where {T <: AbstractArray{<:Int}}
#     return trygetindex(thresholds, 1, 2)
# end
#
# function calculate_Outbreak_start_ews_enddate(thresholds::T)::Union{Try.Ok{Int}, Try.Err} where {T <: AbstractArray{<:Int}}
#     filtered_outbreak_thresholds = filter_outbreak_thresholds(thresholds)
#
#     return trygetindex(filtered_outbreak_thresholds, 1, 1)
# end
#
# function calculate_Outbreak_end_ews_enddate(thresholds::T)::Union{Try.Ok{Int}, Try.Err} where {T <: AbstractArray{<:Int}}
#     filtered_outbreak_thresholds = filter_outbreak_thresholds(thresholds)
#
#     return trygetindex(filtered_outbreak_thresholds, 1, 2)
# end
#
# function calculate_Outbreak_ews_enddate(thresholds::T)::Union{Try.Ok{Int}, Try.Err{BoundsError}} where {T <: AbstractArray{<:Int}}
#     filtered_outbreak_thresholds = filter_outbreak_thresholds(thresholds)
#
#     if size(filtered_outbreak_thresholds, 1) == 0
#         return Try.Err(BoundsError(thresholds, (1, 1)))
#     end
#
#     return Try.Ok(
#         filtered_outbreak_thresholds[1, 1] +
#             (
#             filtered_outbreak_thresholds[1, 2] -
#                 filtered_outbreak_thresholds[1, 1] +
#                 1
#         ) ÷ 2,
#     )
# end

function filter_outbreak_thresholds(
        thresholds::T; thresholds_col = 4
    ) where {T <: AbstractArray{<:Int}}

    return thresholds[(thresholds[:, thresholds_col] .== 1), :]
end

# function tycho_testing_plots(
#         noise_arr_tuple,
#         cases_arr_tuple,
#         plot_cases_tuple,
#         long_plotdata;
#         individual_test_specification = IndividualTestSpecification(1.0, 1.0, 0),
#         noise_specification = NoiseSpecification(PoissonNoise(1.0)),
#         plot_base_path = DrWatson.plotsdir("tycho/testing-plots"),
#         ews_method = Main.Centered,
#         ews_aggregation = 1,
#         ews_bandwidth = 52,
#         ews_lag = 1,
#         ews_metric = "variance",
#         ews_threshold_window = Main.Expanding,
#         ews_threshold_percentile = 0.95,
#         consecutive_thresholds = 2,
#         obsdate = cdc_week_to_date(1990, 3; weekday = 6),
#         sim = 1,
#         return_objects = false,
#         force = true,
#     )
#     @assert length(noise_arr_tuple) == length(cases_arr_tuple) == 3
#     weekly_noise_arr, biweekly_noise_arr, monthly_noise_arr = noise_arr_tuple
#     weekly_cases_arr, biweekly_cases_arr, monthly_cases_arr = cases_arr_tuple
#     weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases =
#         plot_cases_tuple
#
#     @assert length(weekly_noise_arr) == length(weekly_cases_arr)
#     @assert length(biweekly_noise_arr) == length(biweekly_cases_arr)
#     @assert length(monthly_noise_arr) == length(monthly_cases_arr)
#     @assert length(weekly_noise_arr) > length(biweekly_noise_arr) >
#         length(monthly_noise_arr)
#
#     plot_dates =
#         unique(long_plotdata[!, :date]) |>
#         x -> collect(minimum(x):Dates.Day(1):(maximum(x) + Dates.Day(6)))
#
#     weekly_test_arr = create_testing_arrs(
#         weekly_cases_arr,
#         weekly_noise_arr,
#         1.0,
#         individual_test_specification,
#     )
#
#     filled_weekly_test_arr = fill_aggregation_values(
#         weekly_test_arr
#     )
#
#     biweekly_test_arr = create_testing_arrs(
#         biweekly_cases_arr,
#         biweekly_noise_arr,
#         1.0,
#         individual_test_specification,
#     )
#
#     filled_biweekly_test_arr = fill_aggregation_values(
#         biweekly_test_arr; week_aggregation = 2
#     )
#
#     monthly_test_arr = create_testing_arrs(
#         monthly_cases_arr,
#         monthly_noise_arr,
#         1.0,
#         individual_test_specification,
#     )
#
#     filled_monthly_test_arr =
#         fill_aggregation_values(
#         monthly_test_arr;
#         week_aggregation = 4,
#     ) |>
#         x -> x[1:length(plot_dates), :, :]
#
#     test_plot_description = get_test_description(individual_test_specification)
#     noise_magnitude_description = get_noise_magnitude_description(
#         noise_specification
#     )
#
#     test_path_description = "sens-$(individual_test_specification.sensitivity)_spec-$(individual_test_specification.specificity)_lag-$(individual_test_specification.test_result_lag)"
#     noise_path_description = "$(
#         replace(noise_magnitude_description, " " => "_", ": " => "_")
#     )"
#
#     noise_path = joinpath(
#         plot_base_path,
#         getdirpath(noise_specification),
#         "sim_$(sim)",
#     )
#     test_dir = joinpath(noise_path, test_path_description)
#     mkpath(test_dir)
#
#     plotname = "testing-epicurve_$(test_path_description)_$(noise_path_description).png"
#     plotpath = joinpath(test_dir, plotname)
#
#     if !isfile(plotpath) || force
#         plot = tycho_epicurve(
#             plot_dates,
#             filled_weekly_test_arr[:, 5, sim],
#             filled_biweekly_test_arr[:, 5, sim],
#             filled_monthly_test_arr[:, 5, sim];
#             plottitle = "Test Positives",
#             subtitle = "$(test_plot_description), $(noise_magnitude_description) Epicurve",
#         )
#
#         save(
#             plotpath,
#             plot;
#             size = (2200, 1600),
#         )
#
#         Makie.empty!(plot)
#     end
#
#     plotname = "testing-components-epicurve_$(test_path_description)_$(noise_path_description).png"
#     plotpath = joinpath(test_dir, plotname)
#
#     if !isfile(plotpath) || force
#         plot = tycho_test_positive_components_epicurve(
#             plot_dates,
#             (weekly_plot_cases, biweekly_plot_cases, monthly_plot_cases),
#             (
#                 filled_weekly_test_arr[:, :, sim],
#                 filled_biweekly_test_arr[:, :, sim],
#                 filled_monthly_test_arr[:, :, sim],
#             );
#             plottitle = "$(test_plot_description), $(noise_magnitude_description) Component Epicurves",
#         )
#
#         save(
#             plotpath,
#             plot;
#             size = (2200, 1600),
#         )
#
#         Makie.empty!(plot)
#     end
#
#     weekly_obs_index = calculate_ews_enddate(
#         long_plotdata;
#         week_aggregation = 1,
#         obsdate = obsdate,
#     )
#
#     biweekly_obs_index = calculate_ews_enddate(
#         long_plotdata;
#         week_aggregation = 2,
#         obsdate = obsdate,
#     )
#
#     monthly_obs_index = calculate_ews_enddate(
#         long_plotdata;
#         week_aggregation = 4,
#         obsdate = obsdate,
#     )
#
#     weekly_ewsmetric_specification = EWSMetricSpecification(
#         ews_method,
#         ews_aggregation,
#         ews_bandwidth,
#         ews_lag,
#     )
#
#     weekly_cases_ewsmetrics = EWSMetrics(
#         weekly_ewsmetric_specification,
#         weekly_cases_arr[1:weekly_obs_index, :, sim],
#     )
#
#     biweekly_cases_ewsmetrics = EWSMetrics(
#         EWSMetricSpecification(
#             ews_method,
#             ews_aggregation,
#             Int(ews_bandwidth / 2),
#             ews_lag,
#         ),
#         biweekly_cases_arr[1:biweekly_obs_index, :, sim],
#     )
#
#     monthly_cases_ewsmetrics = EWSMetrics(
#         EWSMetricSpecification(
#             ews_method,
#             ews_aggregation,
#             Int(ews_bandwidth / 4),
#             ews_lag,
#         ),
#         monthly_cases_arr[1:monthly_obs_index, :, sim],
#     )
#
#     weekly_test_ewsmetrics =
#         map(
#         k -> EWSMetrics(
#             EWSMetricSpecification(
#                 ews_method,
#                 ews_aggregation,
#                 ews_bandwidth,
#                 ews_lag,
#             ),
#             weekly_test_arr[
#                 1:weekly_obs_index, 5, k,
#             ],
#         ),
#         axes(weekly_test_arr, 3),
#     ) |>
#         x -> StructArray(x)
#
#     biweekly_test_ewsmetrics =
#         map(
#         k -> EWSMetrics(
#             EWSMetricSpecification(
#                 ews_method,
#                 ews_aggregation,
#                 Int(ews_bandwidth / 2),
#                 ews_lag,
#             ),
#             biweekly_test_arr[
#                 1:biweekly_obs_index, 5, k,
#             ],
#         ),
#         axes(biweekly_test_arr, 3),
#     ) |>
#         x -> StructArray(x)
#
#     monthly_test_ewsmetrics =
#         map(
#         k -> EWSMetrics(
#             EWSMetricSpecification(
#                 ews_method,
#                 ews_aggregation,
#                 Int(ews_bandwidth / 4),
#                 ews_lag,
#             ),
#             monthly_test_arr[
#                 1:monthly_obs_index, 5, k,
#             ],
#         ),
#         axes(monthly_test_arr, 3),
#     ) |>
#         x -> StructArray(x)
#
#     ews_metric_tau = ews_metric * "_tau"
#     ews_metric_sym = Symbol(ews_metric)
#     ews_metric_tau_sym = Symbol(ews_metric_tau)
#
#     plotname = "$(ews_metric_tau)_distribution_$(test_path_description)_$(noise_path_description).png"
#     ewspath = joinpath(test_dir, weekly_ewsmetric_specification.dirpath)
#     mkpath(ewspath)
#     plotpath = joinpath(
#         ewspath, plotname
#     )
#
#     if !isfile(plotpath) || force
#         plot = tycho_tau_distribution(
#             (
#                 weekly = weekly_test_ewsmetrics,
#                 biweekly = biweekly_test_ewsmetrics,
#                 monthly = monthly_test_ewsmetrics,
#             ),
#             weekly_cases_ewsmetrics,
#             ews_metric_tau_sym;
#             plottitle = "$(test_plot_description), $(noise_magnitude_description): $(method_string(weekly_ewsmetric_specification.method)) EWS $(ews_metric_tau) Distribution",
#         )
#
#         save(
#             plotpath,
#             plot;
#             size = (2200, 1600),
#         )
#
#         Makie.empty!(plot)
#     end
#
#     weekly_thresholds = expanding_ews_thresholds(
#         weekly_test_ewsmetrics[sim],
#         ews_metric_sym,
#         ews_threshold_window;
#         percentiles = ews_threshold_percentile,
#     )
#
#     biweekly_thresholds = expanding_ews_thresholds(
#         biweekly_test_ewsmetrics[sim],
#         ews_metric_sym,
#         ews_threshold_window;
#         percentiles = ews_threshold_percentile,
#     )
#
#     monthly_thresholds = expanding_ews_thresholds(
#         monthly_test_ewsmetrics[sim],
#         ews_metric_sym,
#         ews_threshold_window;
#         percentiles = ews_threshold_percentile,
#     )
#
#     weekly_detection_index = Try.@? calculate_ews_trigger_index(
#         weekly_thresholds[2],
#         consecutive_thresholds,
#     )
#
#     biweekly_detection_index = Try.@? calculate_ews_trigger_index(
#         biweekly_thresholds[2],
#         consecutive_thresholds,
#     )
#
#     monthly_detection_index = Try.@? calculate_ews_trigger_index(
#         monthly_thresholds[2],
#          consecutive_thresholds,
#     )
#
#     plotname = "$(ews_metric)_$(100 * ews_threshold_percentile)-percentile_thresholds_$(test_path_description)_$(noise_path_description).png"
#     plotpath = joinpath(
#         ewspath, plotname
#     )
#
#     # if !isfile(plotpath) || force
#     plot = tycho_epicurve(
#         plot_dates,
#         (
#             filled_weekly_test_arr[:, 5, sim],
#             filled_biweekly_test_arr[:, 5, sim],
#             filled_monthly_test_arr[:, 5, sim],
#         ),
#         (
#             getproperty(weekly_test_ewsmetrics[sim], ews_metric_sym),
#             getproperty(biweekly_test_ewsmetrics[sim], ews_metric_sym),
#             getproperty(monthly_test_ewsmetrics[sim], ews_metric_sym),
#         ),
#         (
#             weekly_thresholds[2],
#             biweekly_thresholds[2],
#             monthly_thresholds[2],
#         ),
#         (
#             weekly_detection_index,
#             biweekly_detection_index,
#             monthly_detection_index,
#         );
#         plottitle = "Test Positives",
#         subtitle = "$(test_plot_description), $(noise_magnitude_description): $(method_string(weekly_ewsmetric_specification.method)) EWS $(ews_metric) Epicurve",
#         ews_ylabel = ews_metric,
#         threshold_percentile = ews_threshold_percentile,
#         consecutive_thresholds = consecutive_thresholds,
#     )
#
#     save(
#         plotpath,
#         plot;
#         size = (2200, 1600),
#     )
#
#     Makie.empty!(plot)
#     # end
#
#     if return_objects
#         return (
#             plot_dates,
#             weekly_test_arr,
#             filled_weekly_test_arr,
#             biweekly_test_arr,
#             filled_biweekly_test_arr,
#             monthly_test_arr,
#             filled_monthly_test_arr,
#             weekly_obs_index,
#             biweekly_obs_index,
#             monthly_obs_index,
#             weekly_cases_ewsmetrics,
#             biweekly_cases_ewsmetrics,
#             monthly_cases_ewsmetrics,
#             weekly_test_ewsmetrics,
#             biweekly_test_ewsmetrics,
#             monthly_test_ewsmetrics,
#             weekly_thresholds,
#             biweekly_thresholds,
#             monthly_thresholds,
#         )
#     end
#
#     return nothing
# end

function simulation_tau_heatmap_df!(
        ews_df,
        test_ewsmetrics,
        ews_metric;
        individual_test_specification = IndividualTestSpecification(1.0, 1.0, 0),
        ews_metric_specification = EWSMetricSpecification(Centered, 7, 52, 1),
        ews_enddate_type = Reff_start::EWSEndDateType,
        statistic_function = StatsBase.mean,
    )
    @assert names(ews_df) == [
        "ews_metric",
        "test_specification",
        "ews_metric_specification",
        "ews_enddate_type",
        "ews_metric_value",
        "ews_metric_vector",
    ]

    ews_metric_tau_sym = Symbol(ews_metric, "_tau")
    ews_tau = get_tau(
        test_ewsmetrics;
        tau_metric = ews_metric_tau_sym,
        statistic_function = statistic_function,
    )

    push!(
        ews_df,
        (
            ews_metric,
            individual_test_specification,
            ews_metric_specification,
            ews_enddate_type,
            ews_tau,
            getproperty(test_ewsmetrics, ews_metric_tau_sym),
        ),
    )

    return nothing
end

# function tycho_tau_heatmap_df(
#         long_plotdata,
#         cases_arr,
#         noise_arr,
#         test_tuple;
#         week_aggregation = 1,
#         ews_metrics = [
#             :autocorrelation,
#             :autocovariance,
#             :coefficient_of_variation,
#             :index_of_dispersion,
#             :kurtosis,
#             :mean,
#             :skewness,
#             :variance,
#         ],
#         ews_method = Main.Centered,
#         ews_aggregation = 1,
#         ews_bandwidth = 52,
#         ews_lag = 1,
#         obsdate = cdc_week_to_date(1990, 3; weekday = 6),
#         statistic_function = StatsBase.mean,
#     )
#     ews_df = DataFrame(
#         "ews_metric" => String[],
#         "test_specification" => IndividualTestSpecification[],
#         "ews_metric_value" => Float64[],
#         "ews_metric_vector" => Vector{Float64}[],
#     )
#
#     for individual_test_specification in test_tuple
#         test_arr = create_testing_arrs(
#             cases_arr,
#             noise_arr,
#             1.0,
#             individual_test_specification,
#         )
#
#         obs_index = calculate_ews_enddate(
#             long_plotdata;
#             week_aggregation = week_aggregation,
#             obsdate = obsdate,
#         )
#
#         ewsmetric_specification = EWSMetricSpecification(
#             ews_method,
#             ews_aggregation,
#             Int(ews_bandwidth / week_aggregation),
#             ews_lag,
#         )
#
#         test_ewsmetrics =
#             map(
#             k -> EWSMetrics(
#                 ewsmetric_specification,
#                 test_arr[
#                     1:obs_index, 5, k,
#                 ],
#             ),
#             axes(test_arr, 3),
#         ) |>
#             x -> StructArray(x)
#
#         for ews_metric in ews_metrics
#             ews_metric_tau = ews_metric * "_tau"
#             ews_metric_tau_sym = Symbol(ews_metric_tau)
#
#             ews_tau = get_tau(
#                 test_ewsmetrics;
#                 tau_metric = ews_metric_tau_sym,
#                 statistic_function = statistic_function,
#             )
#
#             push!(
#                 ews_df,
#                 (
#                     ews_metric,
#                     individual_test_specification,
#                     ews_tau,
#                     getproperty(test_ewsmetrics, ews_metric_tau_sym),
#                 ),
#             )
#         end
#     end
#     return ews_df
# end

function get_tau(
        ews_metrics;
        tau_metric = :variance_tau,
        statistic_function = StatsBase.mean,
    )
    tau_vector = getproperty(ews_metrics, tau_metric)

    return statistic_function(tau_vector)
end

function calculate_ews_lead_time(
        ews_thresholds;
        week_aggregation = 1,
        consecutive_thresholds = 2,
        output_type = :days,
    )
    threshold_index = Try.@? calculate_ews_trigger_index(
        ews_thresholds, consecutive_thresholds
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
    output_multiplier = if output_type == :days
        7
    elseif output_type == :weeks
        1
    elseif output_type == :months
        7 / 30.5
    elseif output_type == :years
        7 / 365
    else
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
        ews_thresholds::T1,
        consecutive_thresholds = 2,
    ) where {T1 <: AbstractMatrix{<:Bool}}
    reshaped_ews_thresholds = reshape(ews_thresholds, :)

    @assert length(reshaped_ews_thresholds) == length(ews_thresholds[:, 1])

    return calculate_ews_trigger_index(
        reshaped_ews_thresholds,
        consecutive_thresholds,
    )
end

function calculate_ews_trigger_index(
        ews_thresholds::T1,
        consecutive_thresholds = 2,
    ) where {T1 <: AbstractVector{<:Bool}}
    cumulative_thresholds = cumsum(ews_thresholds)
    for (i, v) in pairs(cumulative_thresholds)
        if i > consecutive_thresholds &&
                v >=
                cumulative_thresholds[i - consecutive_thresholds] +
                consecutive_thresholds
            return Try.Ok(i)
        end
    end
    return Try.Err(nothing)
end

# function ews_lead_time_df!(
#         lead_time_df::DataFrame,
#         cases_arr,
#         noise_arr,
#         long_plotdata,
#         individual_test_specification;
#         noise_type = "poisson",
#         noise_magnitude = 1.0,
#         week_aggregation = 1,
#         ews_method = Main.Centered,
#         ews_aggregation = 1,
#         ews_bandwidth = 52,
#         ews_lag = 1,
#         ews_metric = "variance",
#         ews_threshold_window = Main.Expanding,
#         ews_threshold_percentile = 0.95,
#         consecutive_thresholds = 2,
#         obsdate = cdc_week_to_date(1990, 3; weekday = 6),
#         lead_time_units = :days,
#         lead_time_percentile = 0.95,
#         return_objects = false,
#     )
#     test_arr = create_testing_arrs(
#         cases_arr,
#         noise_arr,
#         1.0,
#         individual_test_specification,
#     )
#
#     obs_index = calculate_ews_enddate(
#         long_plotdata;
#         week_aggregation = week_aggregation,
#         obsdate = obsdate,
#     )
#
#     test_ewsmetrics =
#         map(
#         k -> EWSMetrics(
#             EWSMetricSpecification(
#                 ews_method,
#                 ews_aggregation,
#                 Int(ews_bandwidth / week_aggregation),
#                 ews_lag,
#             ),
#             test_arr[
#                 1:obs_index, 5, k,
#             ],
#         ),
#         axes(test_arr, 3),
#     ) |>
#         x -> StructArray(x)
#
#     ews_metric_sym = Symbol(ews_metric)
#
#     weekly_thresholds = Vector{Matrix{Bool}}(undef, length(test_ewsmetrics))
#     ews_lead_time = fill(NaN, length(test_ewsmetrics))
#     for sim in axes(test_ewsmetrics, 1)
#         weekly_thresholds[sim] = expanding_ews_thresholds(
#             test_ewsmetrics[sim],
#             ews_metric_sym,
#             ews_threshold_window;
#             percentiles = ews_threshold_percentile,
#         )[2]
#
#         lead_time = calculate_ews_lead_time(
#             weekly_thresholds[sim];
#             week_aggregation = week_aggregation,
#             consecutive_thresholds = consecutive_thresholds,
#             output_type = lead_time_units,
#         )
#         if !isnothing(lead_time)
#             ews_lead_time[sim] = lead_time
#         end
#     end
#
#     filter!(!isnan, ews_lead_time)
#
#     percentile_tail = (1 - lead_time_percentile) / 2
#
#     if !isempty(ews_lead_time)
#         median_ews_lead_time = StatsBase.median(ews_lead_time)
#         lower_ews_lead_time = StatsBase.quantile(
#             ews_lead_time, percentile_tail
#         )
#         upper_ews_lead_time = StatsBase.quantile(
#             ews_lead_time, 1.0 - percentile_tail
#         )
#     else
#         median_ews_lead_time = NaN
#         lower_ews_lead_time = NaN
#         upper_ews_lead_time = NaN
#     end
#
#     push!(
#         lead_time_df,
#         (
#             noise_type = noise_type,
#             noise_magnitude = noise_magnitude,
#             test_specification = individual_test_specification,
#             week_aggregation = week_aggregation,
#             ews_method = ews_method,
#             lead_time_dist = ews_lead_time,
#             lead_time_median = median_ews_lead_time,
#             lead_time_lower = lower_ews_lead_time,
#             lead_time_upper = upper_ews_lead_time,
#             lead_time_units = String(lead_time_units),
#         ),
#     )
#
#     if return_objects
#         return test_arr, test_ewsmetrics, weekly_thresholds
#     end
#
#     return nothing
# end

function calculate_auc(
        emergent_tau,
        null_tau,
    )
    combined_taus = vcat(null_tau, emergent_tau)

    ranks = StatsBase.tiedrank(-combined_taus)
    n_emergent = length(emergent_tau)
    n_null = length(null_tau)
    sum_null_ranks = sum(ranks[1:n_null])
    return (sum_null_ranks - n_null * (n_null + 1) / 2) /
        (n_emergent * n_null)
end
