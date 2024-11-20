#%%
using DrWatson
@quickactivate "CSDNoise"
using CSDNoise
using Try: Try
using DataFrames
using Dates
using StructArrays
using StatsBase: StatsBase
using CSV: CSV
using Match: @match

include(srcdir("cairomakie-plotting-setup.jl"))
CairoMakie.activate!(; type = "svg")

#%%
ensemble_model_type = ("seasonal-infectivity-import", "tau-leaping")

burnin_years = 5
nyears = 20
burnin_time = 365.0 * burnin_years
ensemble_time_specification = SimTimeParameters(;
    burnin = 365.0 * burnin_years, tmin = 0.0, tmax = 365.0 * nyears,
    tstep = 1.0,
)

ensemble_state_specification = StateParameters(
    500_000,
    Dict(:s_prop => 0.05, :e_prop => 0.0, :i_prop => 0.0, :r_prop => 0.95),
)

mu = calculate_mu(27)
beta_mean = calculate_beta(
    R0, GAMMA, mu, 1, ensemble_state_specification.init_states.N
)
epsilon = calculate_import_rate(
    mu, R0, ensemble_state_specification.init_states.N
)

min_burnin_vaccination_coverage = calculate_vaccination_rate_to_achieve_Reff(
    0.9,
    burnin_years * 2,
    ensemble_state_specification.init_states.S,
    ensemble_state_specification.init_states.N,
    R0,
    mu,
)

max_burnin_vaccination_coverage = 1.0
min_vaccination_coverage = 0.0
max_vaccination_coverage = 0.8

ensemble_dynamics_specification = DynamicsParameterSpecification(
    beta_mean,
    0.0,
    cos,
    SIGMA,
    GAMMA,
    mu,
    27,
    epsilon,
    R0,
    min_burnin_vaccination_coverage,
    max_burnin_vaccination_coverage,
    0.6,
    0.8,
)

null_dynamics_specification = DynamicsParameterSpecification(
    map(
        pn -> getproperty(ensemble_dynamics_specification, pn),
        filter(
            name ->
                name != :min_vaccination_coverage &&
                    name != :max_vaccination_coverage,
            propertynames(ensemble_dynamics_specification),
        ),
    )...,
    nothing,
    nothing,
)

#%%
ensemble_nsims = 100

ensemble_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    ensemble_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

null_specification = EnsembleSpecification(
    ensemble_model_type,
    ensemble_state_specification,
    null_dynamics_specification,
    ensemble_time_specification,
    ensemble_nsims,
)

ensemble_outbreak_specification = OutbreakSpecification(
    5, 30, 500
)

ensemble_single_scenario_inc_file = get_ensemble_file(
    ensemble_specification, ensemble_outbreak_specification
)

null_single_scenario_inc_file = get_ensemble_file(
    null_specification, ensemble_outbreak_specification
)

ensemble_single_incarr = ensemble_single_scenario_inc_file["ensemble_inc_arr"]
null_single_incarr = null_single_scenario_inc_file["ensemble_inc_arr"]
ensemble_single_Reff_thresholds_vec = get_ensemble_file(
    ensemble_specification
)["ensemble_Reff_thresholds_vec"]

#%%
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

ews_metric_specification = EWSMetricSpecification(
    Backward, Day(28), Day(364), 1
)
ews_threshold_burnin = Year(5)
ews_enddate_type = Reff_start
ews_threshold_window = ExpandingThresholdWindow
percent_tested = 1.0

optimization_filepath = get_most_recent_hyperparam_filepath(
    "ews-hyperparam-optimization.jld2",
    outdir("ensemble", "ews-hyperparam-optimization"),
)

previous_optimal_ews_df = load(Try.unwrap(optimization_filepath))["optimal_ews_df"]

test_specification_vec = [
    IndividualTestSpecification(0.8, 0.8, 0),
    IndividualTestSpecification(0.9, 0.9, 0),
    IndividualTestSpecification(0.95, 0.95, 0),
    IndividualTestSpecification(0.96, 0.96, 0),
    IndividualTestSpecification(0.97, 0.97, 0),
    IndividualTestSpecification(0.98, 0.98, 0),
    IndividualTestSpecification(0.99, 0.99, 0),
    IndividualTestSpecification(1.0, 1.0, 0),
]

subset_optimal_df = subset(
    previous_optimal_ews_df,
    :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
    :ews_metric_specification =>
        ByRow(==(ews_metric_specification)),
    :percent_tested => ByRow(==(percent_tested)),
    :ews_metric_specification => ByRow(==(ews_metric_specification)),
    :ews_enddate_type => ByRow(==(ews_enddate_type)),
    :test_specification => ByRow(in(test_specification_vec)),
)

#%%
function plot_noise_filename(
    noise_specification::T
) where {T<:PoissonNoiseSpecification}
    return "poisson_$(noise_specification.noise_mean_scaling)x"
end

function plot_noise_filename(
    noise_specification::T
) where {T<:DynamicalNoiseSpecification}
    mean_vaccination_coverage = mean([
        noise_specification.min_vaccination_coverage,
        noise_specification.max_vaccination_coverage,
    ])
    return "dynamical_$(mean_vaccination_coverage)"
end

function plot_test_filename(test_specification)
    return "sens-$(test_specification.sensitivity)_spec-$(test_specification.specificity)_lag-$(test_specification.test_result_lag)"
end

function extract_tau_auc_metric(
    df,
    test_specification,
    tau_auc_metric,
    header_label;
    rev = true,
)
    return combine(
        sort(
            subset(
                df,
                :test_specification => ByRow(==(test_specification)),
            ),
            order(tau_auc_metric; rev = rev),
        ),
        [:ews_metric, tau_auc_metric] =>
            ByRow(
                (x, y) ->
                    sentencecase(replace(x, "_" => " ")) *
                    " ($(round(y; digits = 2)))",
            ) => header_label,
    )
end

#%%
tables_path = projectdir("manuscript", "manuscript_files", "tables")
mkpath(tables_path)

gdfs = groupby(
    subset_optimal_df,
    [
        :noise_specification
    ],
)

tau_comparison_df = DataFrame(; Rank = collect(1:8))
auc_comparison_df = DataFrame(; Rank = collect(1:8))
auc_magnitude_comparison_df = DataFrame(; Rank = collect(1:8))
alert_auc_magnitude_comparison_df = DataFrame(; Rank = collect(1:8))
combined_survival_df = DataFrame()
for (noise_num, gdf) in enumerate(gdfs)
    @assert length(unique(gdf.noise_specification)) == 1
    noise_specification = gdf.noise_specification[1]
    noise_description = noise_table_description(noise_specification)

    optimal_heatmap_df = optimal_ews_heatmap_df(
        gdf;
        tiebreaker_preference = "specificity",
        optimal_grouping_parameters = [
            :noise_specification,
            :test_specification,
            :percent_tested,
            :ews_metric_specification,
            :ews_enddate_type,
            :ews_threshold_window,
            :ews_threshold_burnin,
            :ews_metric,
        ],
    )

    optimal_heatmap_plot = optimal_ews_heatmap_plot(
        optimal_heatmap_df;
        subtitle = noise_description,
        colormap = :Blues,
        colorrange = (0.5, 0.8),
        textcolorthreshold = 0.68,
    )

    optimal_heatmap_plot_name =
        "optimal_heatmap_" * plot_noise_filename(noise_specification) *
        ".svg"

    heatmap_plotdir = projectdir(
        "manuscript", "manuscript_files", "plots", "optimal-threshold-heatmaps"
    )

    mkpath(heatmap_plotdir)

    save(
        joinpath(heatmap_plotdir, optimal_heatmap_plot_name),
        optimal_heatmap_plot;
        size = (1700, 1000),
    )

    test_df = subset(
        optimal_heatmap_df,
        :ews_metric_specification =>
            ByRow(==(ews_metric_specification)),
        :ews_enddate_type => ByRow(==(ews_enddate_type)),
        :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
        :ews_threshold_window => ByRow(==(ews_threshold_window)),
        :noise_specification => ByRow(==(noise_specification)),
        :percent_tested => ByRow(==(percent_tested)),
    )

    for (i, survival_ews_metric) in pairs(ews_metrics)
        survival_df, (
        vec_of_testarr,
        vec_of_null_testarr,
        vec_of_ews_vals_vec,
        vec_of_null_ews_vals_vec,
        vec_of_exceed_thresholds,
        vec_of_null_exceed_thresholds,
        vec_of_threshold_percentiles,
        vec_of_null_threshold_percentiles,
        vec_of_detection_index_vec,
        vec_of_null_detection_index_vec
), noisearr = simulate_ews_survival_data(
            test_df,
            ensemble_specification,
            ensemble_single_incarr,
            null_single_incarr,
            ensemble_single_Reff_thresholds_vec;
            ews_metric = survival_ews_metric,
        )

        global combined_survival_df = vcat(
            combined_survival_df,
            survival_df;
            cols = :union,
        )

        if i == 1
            perfect_test_df = DataFrame(; Rank = collect(1:length(ews_metrics)))
            perfect_test_accuracy = extract_tau_auc_metric(
                optimal_heatmap_df,
                IndividualTestSpecification(1.0, 1.0, 0),
                :accuracy,
                "All Noise - Accuracy";
                rev = true,
            )

            for (start_ind, timelength_plotdir, timelength_label) in zip(
                [1, Int64(round(Dates.days(Year(5))))],
                ["full-length", "after-burnin"],
                ["Full Time Series", "After Burn-in Period"],
            )
                auc_df = DataFrame(
                    "ews_metric" => String[],
                    "test_specification" => IndividualTestSpecification[],
                    "auc" => Float64[],
                    "auc_magnitude" => Float64[],
                    "emergent_tau" => Vector{Float64}[],
                    "null_tau" => Vector{Float64}[],
                )

                unique_tests = unique(survival_df.test_specification)
                subset_start_ind = Int64(round(Dates.days(Year(5))))

                for test_ind in eachindex(vec_of_ews_vals_vec)
                    test_specification = unique_tests[test_ind]
                    subset_testarr = vec_of_testarr[test_ind][
                        :, 5, :,
                    ]
                    subset_null_testarr = vec_of_null_testarr[test_ind][
                        :, 5, :,
                    ]

                    enddate_vec = map(axes(subset_testarr, 2)) do sim
                        Try.unwrap(
                            calculate_ews_enddate(
                                ensemble_single_Reff_thresholds_vec[sim],
                                ews_enddate_type,
                            ),
                        )
                    end

                    emergent_sv =
                        map(axes(subset_testarr, 2)) do sim
                            enddate = enddate_vec[sim]
                            EWSMetrics(
                                ews_metric_specification,
                                subset_testarr[start_ind:enddate, sim],
                            )
                        end |>
                        v -> StructVector(v)

                    null_sv =
                        map(axes(subset_null_testarr, 2)) do sim
                            enddate = enddate_vec[sim]
                            EWSMetrics(
                                ews_metric_specification,
                                subset_null_testarr[start_ind:enddate, sim],
                            )
                        end |>
                        v -> StructVector(v)

                    for metric in ews_metrics
                        metric_tau = Symbol(metric * "_tau")
                        emergent_tau = getproperty(emergent_sv, metric_tau)
                        null_tau = getproperty(null_sv, metric_tau)
                        auc = calculate_auc(
                            emergent_tau,
                            null_tau,
                        )
                        auc_magnitude = abs(auc - 0.5)
                        push!(
                            auc_df,
                            (
                                ews_metric = metric,
                                test_specification = test_specification,
                                auc = auc,
                                auc_magnitude,
                                emergent_tau = emergent_tau,
                                null_tau = null_tau,
                            ),
                        )
                    end
                end

                auc_heatmap = tau_auc_heatmap(
                    auc_df,
                    :auc,
                    :auc;
                    colormap = :RdBu,
                    colorrange = [0, 1.0],
                    textcolorthreshold = (0.3, 0.7),
                )

                tau_auc_heatmap_plot_name =
                    "tau_auc-heatmap_" *
                    plot_noise_filename(noise_specification) *
                    ".svg"

                tau_auc_heatmap_plotdir = projectdir(
                    "manuscript",
                    "manuscript_files",
                    "plots",
                    "tau_auc-heatmaps",
                    timelength_plotdir,
                )

                mkpath(tau_auc_heatmap_plotdir)

                save(
                    joinpath(
                        tau_auc_heatmap_plotdir, tau_auc_heatmap_plot_name
                    ),
                    auc_heatmap,
                )

                auc_magnitude_heatmap = tau_auc_heatmap(
                    auc_df,
                    :auc_magnitude,
                    :auc_magnitude;
                    colormap = :Blues,
                    colorrange = [0, 0.5],
                    textcolorthreshold = 0.3,
                )

                tau_auc_magnitude_heatmap_plot_name =
                    "tau_auc-magnitude-heatmap_" *
                    plot_noise_filename(noise_specification) *
                    ".svg"

                tau_auc_magnitude_heatmap_plotdir = projectdir(
                    "manuscript",
                    "manuscript_files",
                    "plots",
                    "tau_auc-magnitude-heatmaps",
                    timelength_plotdir,
                )

                mkpath(tau_auc_magnitude_heatmap_plotdir)

                save(
                    joinpath(
                        tau_auc_magnitude_heatmap_plotdir,
                        tau_auc_magnitude_heatmap_plot_name,
                    ),
                    auc_magnitude_heatmap,
                )

                emergent_tau_df = hcat(
                    select(auc_df, Cols(1:2)),
                    combine(
                        auc_df,
                        :emergent_tau =>
                            ByRow(x -> mean(x)) => :ews_metric_value,
                    ),
                )

                null_tau_df = hcat(
                    select(auc_df, Cols(1:2)),
                    combine(
                        auc_df,
                        :null_tau => ByRow(x -> mean(x)) => :ews_metric_value,
                    ),
                )

                emergent_tau_heatmap = tycho_tau_heatmap_plot(
                    emergent_tau_df
                )

                emergent_tau_heatmap_plot_name =
                    "emergent-tau-heatmap_" *
                    plot_noise_filename(noise_specification) *
                    ".svg"

                emergent_tau_heatmap_plotdir = projectdir(
                    "manuscript",
                    "manuscript_files",
                    "plots",
                    "tau_heatmaps",
                    "emergent",
                    timelength_plotdir,
                )

                mkpath(emergent_tau_heatmap_plotdir)

                save(
                    joinpath(
                        emergent_tau_heatmap_plotdir,
                        emergent_tau_heatmap_plot_name,
                    ),
                    emergent_tau_heatmap,
                )

                null_tau_heatmap = tycho_tau_heatmap_plot(
                    null_tau_df
                )

                null_tau_heatmap_plot_name =
                    "null-tau-heatmap_" *
                    plot_noise_filename(noise_specification) *
                    ".svg"

                null_tau_heatmap_plotdir = projectdir(
                    "manuscript",
                    "manuscript_files",
                    "plots",
                    "tau_heatmaps",
                    "null",
                    timelength_plotdir,
                )

                mkpath(null_tau_heatmap_plotdir)

                save(
                    joinpath(
                        null_tau_heatmap_plotdir, null_tau_heatmap_plot_name
                    ),
                    null_tau_heatmap,
                )

                if noise_num == 1
                    perfect_test_tau = extract_tau_auc_metric(
                        emergent_tau_df,
                        IndividualTestSpecification(1.0, 1.0, 0),
                        :ews_metric_value,
                        "$timelength_label Tau";
                        rev = true,
                    )

                    perfect_test_auc = extract_tau_auc_metric(
                        auc_df,
                        IndividualTestSpecification(1.0, 1.0, 0),
                        :auc,
                        "$timelength_label AUC";
                        rev = true,
                    )

                    perfect_test_auc_magnitude = extract_tau_auc_metric(
                        auc_df,
                        IndividualTestSpecification(1.0, 1.0, 0),
                        :auc_magnitude,
                        "$timelength_label |AUC-0.5|";
                        rev = true,
                    )

                    perfect_test_df = hcat(
                        perfect_test_df,
                        perfect_test_tau,
                        perfect_test_auc_magnitude,
                    )

                    if timelength_plotdir == "after-burnin"
                        global tau_comparison_df = hcat(
                            tau_comparison_df,
                            rename(perfect_test_tau, 1 => "All Noise"),
                        )
                        global auc_comparison_df = hcat(
                            auc_comparison_df,
                            rename(perfect_test_auc, 1 => "All Noise"),
                        )
                        global auc_magnitude_comparison_df = hcat(
                            auc_magnitude_comparison_df,
                            rename(
                                perfect_test_auc_magnitude, 1 => "All Noise"
                            ),
                        )
                        global alert_auc_magnitude_comparison_df = hcat(
                            alert_auc_magnitude_comparison_df,
                            rename(
                                perfect_test_auc_magnitude,
                                1 => "All Noise - |AUC-0.5|",
                            ),
                            perfect_test_accuracy,
                        )
                    end
                end

                if timelength_plotdir == "after-burnin"
                    rdt_tau = extract_tau_auc_metric(
                        emergent_tau_df,
                        IndividualTestSpecification(0.9, 0.9, 0),
                        :ews_metric_value,
                        "$noise_description";
                        rev = true,
                    )

                    global tau_comparison_df = hcat(
                        tau_comparison_df,
                        rdt_tau,
                    )

                    global rdt_auc = extract_tau_auc_metric(
                        auc_df,
                        IndividualTestSpecification(0.9, 0.9, 0),
                        :auc,
                        "$noise_description";
                        rev = true,
                    )

                    global rdt_auc_magnitude = extract_tau_auc_metric(
                        auc_df,
                        IndividualTestSpecification(0.9, 0.9, 0),
                        :auc_magnitude,
                        "$noise_description";
                        rev = true,
                    )

                    global auc_comparison_df = hcat(
                        auc_comparison_df,
                        rdt_auc,
                    )

                    global auc_magnitude_comparison_df = hcat(
                        auc_magnitude_comparison_df,
                        rdt_auc_magnitude,
                    )
                end
            end

            if noise_num == 1
                select!(perfect_test_df, :Rank, Cols(r"Tau"), Cols(r"AUC"))

                CSV.write(
                    joinpath(tables_path, "perfect-test_tau-auc.csv"),
                    perfect_test_df,
                )
            end
        end
    end

    rdt_accuracy = extract_tau_auc_metric(
        optimal_heatmap_df,
        IndividualTestSpecification(0.9, 0.9, 0),
        :accuracy,
        "$noise_description";
        rev = true,
    )

    global alert_auc_magnitude_comparison_df = hcat(
        alert_auc_magnitude_comparison_df,
        rdt_accuracy,
    )
end

#%%
select!(tau_comparison_df, :Rank, Cols("All Noise"), All())
CSV.write(
    joinpath(tables_path, "tau-comparison.csv"),
    tau_comparison_df,
)

select!(auc_magnitude_comparison_df, :Rank, Cols("All Noise"), All())
CSV.write(
    joinpath(tables_path, "auc-magnitude-comparison.csv"),
    auc_magnitude_comparison_df,
)

select!(auc_comparison_df, :Rank, Cols("All Noise"), All())
CSV.write(
    joinpath(tables_path, "auc-comparison.csv"),
    auc_comparison_df,
)

select!(alert_auc_magnitude_comparison_df, :Rank, Cols(r"All Noise"), All())
CSV.write(
    joinpath(tables_path, "alert-accuracy-auc-comparison.csv"),
    alert_auc_magnitude_comparison_df,
)

#%%
lineplot_df = similar(gdfs[1], 0)
for gdf in gdfs,
    ewsmetric in ["autocovariance", "variance", "mean", "index_of_dispersion"]

    prepare_line_plot_df!(
        lineplot_df,
        gdf,
        ewsmetric,
        test_specification_vec,
    )
end

accuracy_line_plot = line_plot(lineplot_df)

line_plotdir = projectdir("manuscript", "manuscript_files", "plots")

save(
    joinpath(
        line_plotdir,
        "accuracy-line-plot.svg",
    ),
    accuracy_line_plot,
)

#%%
supplemental_lineplot_df = similar(gdfs[1], 0)
for gdf in gdfs,
    ewsmetric in
    ["autocorrelation", "coefficient_of_variation", "skewness", "kurtosis"]

    prepare_line_plot_df!(supplemental_lineplot_df, gdf, ewsmetric)
end

supplemental_accuracy_line_plot = line_plot(supplemental_lineplot_df)

supplemental_line_plotdir = projectdir(
    "manuscript", "supplemental_files", "plots"
)

save(
    joinpath(
        supplemental_line_plotdir,
        "accuracy-line-plot.svg",
    ),
    supplemental_accuracy_line_plot,
)

#%%
survival_plotdir = projectdir(
    "manuscript",
    "manuscript_files",
    "plots",
    "survival",
)
supplemental_survival_plotdir = projectdir(
    "manuscript",
    "supplemental_files",
    "plots",
    "survival",
)
mkpath(survival_plotdir)
mkpath(supplemental_survival_plotdir)

for metric in ews_metrics
    survival_plot = ews_survival_plot(
        subset(combined_survival_df, :ews_metric => ByRow(==(metric)));
        noise_specification_vec = [
            PoissonNoiseSpecification(1.0),
            PoissonNoiseSpecification(7.0),
            DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.8734),
            DynamicalNoiseSpecification(5.0, 7, 14, "in-phase", 0.15, 0.102),
        ],
        test_specification_vec = [
            IndividualTestSpecification(1.0, 1.0, 0),
            IndividualTestSpecification(0.9, 0.9, 0),
        ],
        linestyle_vec = [:solid, :dot],
    )

    survival_plot_name = "survival_" *
                         "ews-" * metric * ".svg"

    dir = if metric == "autocovariance"
        survival_plotdir
    else
        supplemental_survival_plotdir
    end
    save(
        joinpath(dir, survival_plot_name),
        survival_plot,
    )
end
