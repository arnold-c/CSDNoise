#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using ColorSchemes
using CSV: CSV
using DataFrames: DataFrames
using RCall
using StructArrays
using SumTypes

using CSDNoise

include(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
ensemble_single_individual_test_spec_no_lag = IndividualTestSpecification(
    1.0, 1.0, 0
)

ensemble_single_individual_test_spec_lag = IndividualTestSpecification(
    1.0, 1.0, 14
)

ensemble_single_outbreak_detection_spec = OutbreakDetectionSpecification(
    10, 7, 0.6, 0.2, "inferred_movingavg"
)

ews_metric_specification = EWSMetricSpecification("backward", 35, 1)

#%%
ensemble_single_seir_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_seir_arr"]

ensemble_single_Reff_arr = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_arr"]

ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_periodsum_vecs = ensemble_single_scenario_inc_file["ensemble_thresholds_vec"]

ensemble_single_inc_ews = get_ensemble_file(
    ensemble_specification,
    ensemble_outbreak_specification,
    ews_metric_specification,
)["inc_ewsmetrics"]

#%%
mkpath(plotsdir("ensemble/single-scenario"))

ensemble_single_scenario_incidence_prevalence_plot = incidence_prevalence_plot(
    ensemble_single_incarr,
    ensemble_single_seir_arr,
    ensemble_single_periodsum_vecs,
    ensemble_time_specification;
    threshold = ensemble_outbreak_specification.outbreak_threshold,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_incidence_prevalence.png",
    ),
    ensemble_single_scenario_incidence_prevalence_plot,
)

#%%
for sim_num in [10, 20, 55]
    ensemble_single_scenario_Reffective_plot = Reff_plot(
        ensemble_single_incarr,
        ensemble_single_Reff_arr,
        ensemble_single_Reff_thresholds_vec,
        ensemble_single_periodsum_vecs,
        ensemble_time_specification;
        sim = sim_num,
        threshold = ensemble_outbreak_specification.outbreak_threshold,
    )

    save(
        plotsdir(
            "ensemble/single-scenario/ensemble_single_scenario_Reffective_sim_$(sim_num).png",
        ),
        ensemble_single_scenario_Reffective_plot,
    )
end

#%%
noise_specification = ensemble_noise_specification_vec[1]

scenario_specification_no_lag = ScenarioSpecification(
    ensemble_specification,
    ensemble_outbreak_specification,
    noise_specification,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_no_lag,
    ews_metric_specification,
)

# scenario_specification_lag = ScenarioSpecification(
#     ensemble_specification,
#     ensemble_outbreak_specification,
#     noise_specification,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_single_individual_test_spec_lag,
#     ews_metric_specification,
# )

ensemble_solution_dict_no_lag = get_ensemble_file(scenario_specification_no_lag)
# ensemble_solution_dict_lag = get_ensemble_file(scenario_specification_lag)

OT_chars_no_lag = ensemble_solution_dict_no_lag["OT_chars"]
# OT_chars_lag = ensemble_solution_dict_lag["OT_chars"]

ewsmetrics_no_lag = ensemble_solution_dict_no_lag["test_ewsmetrics"]
# ewsmetrics_lag = ensemble_solution_dict_lag["test_ewsmetrics"]

noisearr, poisson_noise_prop = create_noise_arr(
    noise_specification,
    ensemble_single_incarr;
    ensemble_specification = ensemble_specification,
    seed = 1234,
)
noisedir = getdirpath(noise_specification)

testarr_no_lag, test_movingvg_arr_no_lag, inferred_positives_arr_no_lag = create_testing_arrs(
    ensemble_single_incarr,
    noisearr,
    ensemble_single_outbreak_detection_spec,
    ensemble_single_individual_test_spec_no_lag,
    ensemble_time_specification,
    ews_metric_specification,
)[[1, 3, 4]]

# testarr_lag, test_movingvg_arr_lag, inferred_positives_arr_lag = create_testing_arrs(
#     ensemble_single_incarr,
#     noisearr,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_single_individual_test_spec_lag,
#     ensemble_time_specification,
#     ews_metric_specification,
# )[[1, 3, 4]]

#%%
inferred_positive_no_lag_plot = lines(
    ensemble_time_specification.trange,
    ensemble_single_incarr[:, 1, 1];
    color = :orange,
)
lines!(
    ensemble_time_specification.trange,
    inferred_positives_arr_no_lag[:, 1];
    color = :black,
)

save(
    plotsdir(
        "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_no_lag.png",
    ),
    inferred_positive_no_lag_plot,
)

#%%
# inferred_positive_lag_plot = lines(
#     ensemble_time_specification.trange,
#     ensemble_single_incarr[:, 1, 1];
#     color = :orange,
# )
# lines!(
#     ensemble_time_specification.trange,
#     inferred_positives_arr_lag[:, 1];
#     color = :black,
# )
#
# save(
#     plotsdir(
#         "ensemble/single-scenario/ensemble_single_scenario_inferred_positives_lag.png",
#     ),
#     inferred_positive_lag_plot,
# )

#%%
# compare_inferred_positives_plot = lines(
#     ensemble_time_specification.trange,
#     inferred_positives_arr_no_lag[:, 1];
#     color = :orange,
# )
# lines!(
#     ensemble_time_specification.trange,
#     inferred_positives_arr_lag[:, 1];
#     color = :black,
# )
#
# save(
#     plotsdir(
#         "ensemble/single-scenario/ensemble_single_scenario_compare_inferred_positives.png",
#     ),
#     compare_inferred_positives_plot,
# )

#%%
sim_num = 1

#%%
# Can use RCall to @rput, but write to csv so can put R code in own file
# that can be run independently (as well as provide linting etc)
CSV.write(
    outdir("incidence-array_sim-$(sim_num).csv"),
    DataFrames.DataFrame(
        ensemble_single_incarr[:, :, sim_num],
        [:incidence, :above_threshold, :outbreak],
    ),
)

#%%
R"""
source(here::here("scripts","spaero-ews.R"))
"""

@rget spaero_ews_backward spaero_ews_centered

#%%
function centered_mean(timeseries, bw)
    tlength = length(timeseries)
    mean_vec = zeros(Float64, length(timeseries))
    centered_mean!(mean_vec, tlength, timeseries, bw)
    return mean_vec
end
function centered_mean!(mean_vec, tlength, timeseries, bw)
    @inbounds for i in eachindex(timeseries)
        if i < bw
            mean_vec[i] = mean(@view(timeseries[begin:(i + bw - 1)]))
        elseif i + bw > tlength
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):end]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):(i + bw - 1)]))
        end
    end
    return nothing
end

function centered_moment(timeseries, moment, bw)
    return centered_moment(
        centered_mean(timeseries, bw), timeseries, moment, bw
    )
end
function centered_moment(mean_timeseries, timeseries, moment, bw)
    diff = (timeseries .- mean_timeseries) .^ moment
    return centered_mean(diff, bw)
end
@benchmark centered_moment(testinc, 2, 30)

function compare_against_spaero(spaero_ews, my_ews)
    df = DataFrames.DataFrame([spaero_ews my_ews], [:spaero, :mine])
    df.absdiff = abs.(df.spaero .- df.mine)
    return df
end

function filter_spaero_comparison(spaero_ews, my_ews; tolerance = 1e-13)
    df = compare_against_spaero(spaero_ews, my_ews)
    return filter_spaero_comparison(df; tolerance = tolerance)
end

function filter_spaero_comparison(df; tolerance = 1e-13)
    return DataFrames.subset(
        df, :absdiff => x -> x .> tolerance; skipmissing = true
    )
end

testinc = ensemble_single_incarr[:, 1, sim_num]

#%%
mean_df_centered = compare_against_spaero(
    spaero_ews_centered.mean, centered_mean(testinc, 30)
)
filter_spaero_comparison(mean_df_centered)

var_df_centered = compare_against_spaero(
    spaero_ews_centered.variance, centered_moment(testinc, 2, 30)
)
filter_spaero_comparison(var_df_centered)

#%%
function backward_mean(timeseries, bw)
    tlength = length(timeseries)
    mean_vec = zeros(Float64, tlength)
    @inbounds for i in eachindex(timeseries)
        if i < bw
            mean_vec[i] = mean(@view(timeseries[begin:i]))
        else
            mean_vec[i] = mean(@view(timeseries[(i - bw + 1):i]))
        end
    end
    return mean_vec
end

function backward_moment(timeseries, moment, bw)
    return backward_moment(
        timeseries, backward_mean(timeseries, bw), moment, bw
    )
end

function backward_moment(timeseries, mean_timeseries, moment, bw)
    diff = (timeseries .- mean_timeseries) .^ moment
    return backward_mean(diff, bw)
end

mean_df_backward = compare_against_spaero(
    spaero_ews_backward.mean,
    backward_mean(testinc, 30)
)
filter_spaero_comparison(mean_df_backward)

var_df_backward = compare_against_spaero(
    spaero_ews_backward.variance,
    centered_moment(testinc, 2, 30)
)
filter_spaero_comparison(var_df_backward)

#%%
function cov(mean_func, moment_func, timeseries, bandwidth)
    return sqrt.(moment_func(timeseries, 2, bandwidth)) ./
           mean_func(timeseries, bandwidth)
end
@benchmark cov(centered_mean, centered_moment, testinc, 30)

@sum_type EWSMethod begin
    Backward
    Centered
end

function window_functions(
    method::EWSMethod
)::Tuple{
    Union{typeof(centered_mean),typeof(backward_mean)},
    Union{typeof(centered_moment),typeof(backward_moment)},
}
    mean_func = @cases method begin
        Backward => backward_mean
        Centered => centered_mean
    end
    moment_func = @cases method begin
        Backward => backward_moment
        Centered => centered_moment
    end
    return mean_func, moment_func
end

@report_opt window_functions(Backward)
@code_warntype window_functions(Backward)

function cov(method::EWSMethod, timeseries, bandwidth)
    mean_func, moment_func = window_functions(method)
    return sqrt.(moment_func(timeseries, 2, bandwidth)) ./
           mean_func(timeseries, bandwidth)
end

@benchmark cov(Centered, testinc, 30)

cov_df_centered = compare_against_spaero(
    spaero_ews_centered.coefficient_of_variation,
    cov(Centered, testinc, 30)
)
filter_spaero_comparison(cov_df_centered)

#%%
function iod(method::EWSMethod, timeseries, bandwidth)
    mean_func, moment_func = window_functions(method)
    return (moment_func(timeseries, 2, bandwidth)) ./
           mean_func(timeseries, bandwidth)
end

@benchmark iod(Centered, testinc, 30)

iod_df_centered = compare_against_spaero(
    spaero_ews_centered.index_of_dispersion,
    iod(Centered, testinc, 30)
)
filter_spaero_comparison(iod_df_centered)

#%%
function skew(method::EWSMethod, timeseries, bandwidth)
    moment_func = window_functions(method)[2]
    return moment_func(timeseries, 3, bandwidth) ./
           moment_func(timeseries, 2, bandwidth) .^ 1.5
end

skew_df_centered = compare_against_spaero(
    spaero_ews_centered.skewness,
    skew(Centered, testinc, 30)
)
filter_spaero_comparison(skew_df_centered)

#%%
function kurtosis(method::EWSMethod, timeseries, bandwidth)
    moment_func = window_functions(method)[2]
    return moment_func(timeseries, 4, bandwidth) ./
           moment_func(timeseries, 2, bandwidth) .^ 2
end

kurtosis_df_centered = compare_against_spaero(
    spaero_ews_centered.kurtosis,
    kurtosis(Centered, testinc, 30)
)
filter_spaero_comparison(kurtosis_df_centered; tolerance = 1e-10)

#%%
function autocov(method::EWSMethod, timeseries, bandwidth)
    mean_func = window_functions(method)[1]
    meandiff = timeseries .- mean_func(timeseries, bandwidth)
    worker = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i == 1
            continue
        end
        worker[i] = meandiff[i] * meandiff[i - 1]
    end
    worker = mean_func(worker, bandwidth)
    worker[1] = NaN
    return worker
end

autocov_df_centered = compare_against_spaero(
    spaero_ews_centered.autocovariance,
    autocov(Centered, testinc, 30)
)
filter_spaero_comparison(autocov_df_centered; tolerance = 1e-10)

#%%
function autocor(method::EWSMethod, timeseries, bandwidth)
    mean_func, moment_func = window_functions(method)
    sd = sqrt.(moment_func(timeseries, 2, bandwidth))
    lagged_sd = zeros(Float64, length(timeseries))
    @inbounds for i in eachindex(timeseries)
        if i == 1
            lagged_sd[i] = NaN
            continue
        end
        lagged_sd[i] = sd[i - 1]
    end
    return autocov(method, timeseries, bandwidth) ./ (sd .* lagged_sd)
end

autocor_df_centered = compare_against_spaero(
    spaero_ews_centered.autocorrelation,
    autocor(Centered, testinc, 30)
)
filter_spaero_comparison(autocor_df_centered; tolerance = 1e-10)

#%%
spaero_ensemble_single_inc_ews = StructArray([
    EWSMetrics(
        ensemble_time_specification.tstep,
        ews_metric_specification,
        spaero_ensemble_single_inc_ews_arr.mean,
        spaero_ensemble_single_inc_ews_arr.variance,
        spaero_ensemble_single_inc_ews_arr.coefficient_of_variation,
        spaero_ensemble_single_inc_ews_arr.index_of_dispersion,
        spaero_ensemble_single_inc_ews_arr.skewness,
        spaero_ensemble_single_inc_ews_arr.kurtosis,
        convert(
            Vector{Float64},
            replace!(
                spaero_ensemble_single_inc_ews_arr.autocovariance,
                missing => NaN,
            ),
        ),
        convert(
            Vector{Float64},
            replace!(
                spaero_ensemble_single_inc_ews_arr.autocorrelation,
                missing => NaN,
            ),
        ),
    ),
])

#%%
mean_df = DataFrames.DataFrame(
    [
        ensemble_single_inc_ews[1].mean,
        spaero_ensemble_single_inc_ews[1].mean,
    ],
    [:mine, :spaero],
)

mean_df.diff = mean_df.mine .- mean_df.spaero;
mean_df

lines(mean_df.diff)

#%%
var_df = DataFrames.DataFrame(
    [
        ensemble_single_inc_ews[1].variance,
        spaero_ensemble_single_inc_ews[1].variance,
    ],
    [:mine, :spaero],
)

var_df.diff = var_df.mine .- var_df.spaero;
var_df

lines(var_df.diff)

#%%
ewsmetrics = [
    :autocorrelation,
    :autocovariance,
    :coefficient_of_variation,
    :index_of_dispersion,
    :kurtosis,
    :mean,
    :skewness,
    :variance,
]

basedir = plotsdir(
    "ensemble/single-scenario/ewsmetrics/$(ews_metric_specification.dirpath)"
)
mkpath(basedir)
testdir = joinpath(basedir, "test")
mkpath(testdir)
incdir = joinpath(basedir, "incidence")
spaero_comparison_dir = joinpath(incdir, "spaero-comparison")
mkpath(spaero_comparison_dir)

for ewsmetric in ewsmetrics
    inc_ews_metric_plot = Reff_ews_plot(
        ensemble_single_incarr,
        ensemble_single_Reff_arr,
        ensemble_single_Reff_thresholds_vec,
        ensemble_single_inc_ews,
        ewsmetric,
        ensemble_single_periodsum_vecs,
        ensemble_time_specification;
        sim = sim_num,
        threshold = ensemble_outbreak_specification.outbreak_threshold,
        plottitle = "Incidence $(ewsmetric): $(ews_metric_specification.method)",
    )

    save(
        joinpath(
            incdir,
            "ensemble_single_scenario_incidence_metric_$(String(ewsmetric)).png",
        ),
        inc_ews_metric_plot,
    )

    spaero_comparison_inc_ews_metric_plot = Reff_ews_plot(
        ensemble_single_incarr,
        ensemble_single_Reff_arr,
        ensemble_single_Reff_thresholds_vec,
        ensemble_single_inc_ews,
        spaero_ensemble_single_inc_ews,
        ewsmetric,
        ensemble_single_periodsum_vecs,
        ensemble_time_specification;
        sim = sim_num,
        threshold = ensemble_outbreak_specification.outbreak_threshold,
        plottitle = "Incidence $(ewsmetric): $(ews_metric_specification.method)",
    )

    save(
        joinpath(
            spaero_comparison_dir,
            "ensemble_single_scenario_incidence_metric_$(String(ewsmetric)).png",
        ),
        spaero_comparison_inc_ews_metric_plot,
    )

    for (ewsmetric_sa, lag_label) in
        zip((
        # ewsmetrics_lag,
        [ewsmetrics_no_lag]
    ), (
        # "lag",
        ["no_lag"]
    ))
        test_ews_metric_plot = Reff_ews_plot(
            ensemble_single_incarr,
            ensemble_single_Reff_arr,
            ensemble_single_Reff_thresholds_vec,
            ewsmetric_sa,
            ewsmetric,
            ensemble_single_periodsum_vecs,
            ensemble_time_specification;
            sim = sim_num,
            threshold = ensemble_outbreak_specification.outbreak_threshold,
            plottitle = "Test Positive $(ewsmetric)\t$(lag_label): $(ews_metric_specification.method)",
        )

        save(
            joinpath(
                testdir,
                "ensemble_single_scenario_testpositive_metric_$(String(ewsmetric))_$lag_label.png",
            ),
            test_ews_metric_plot,
        )
    end
end

# plot_all_single_scenarios(
#     noisearr,
#     poisson_noise_prop,
#     noisedir,
#     OT_chars,
#     ensemble_single_incarr,
#     testarr,
#     test_movingvg_arr,
#     ensemble_single_individual_test_spec_no_lag,
#     ensemble_single_outbreak_detection_spec,
#     ensemble_time_specification,
# )

# GC.gc(true)
# @info "Finished plotting the single scenario for $(noisedir)"
# println("=================================================================")
