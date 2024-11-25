function simulate_and_plot_ews_survival(
    optimal_heatmap_df,
    ews_metric_specification,
    ews_enddate_type,
    ews_threshold_burnin,
    ews_threshold_window,
    noise_specification,
    ensemble_specification,
    individual_test_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec;
    ews_metric = "mean",
    plottitle = "Survival",
    subtitle = "",
)
    subset_df = subset(
        optimal_heatmap_df,
        :ews_metric_specification => ByRow(==(ews_metric_specification)),
        :ews_enddate_type => ByRow(==(ews_enddate_type)),
        :ews_threshold_burnin => ByRow(==(ews_threshold_burnin)),
        :ews_threshold_window => ByRow(==(ews_threshold_window)),
        :noise_specification => ByRow(==(noise_specification)),
    )

    return simulate_and_plot_ews_survival(
        subset_df,
        ews_metric_specification,
        ews_threshold_burnin,
        ensemble_specification,
        individual_test_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = ews_metric,
        plottitle = plottitle,
        subtitle = subtitle,
    )
end

function simulate_and_plot_ews_survival(
    subset_df,
    ews_metric_specification,
    ews_threshold_burnin,
    ensemble_specification,
    individual_test_specification,
    ensemble_single_incarr,
    null_single_incarr,
    ensemble_single_Reff_thresholds_vec;
    ews_metric = "mean",
    plottitle = "Survival",
    subtitle = "",
)
    test_subset_df = subset(
        subset_df,
        :test_specification => ByRow(==(individual_test_specification)),
    )

    survival_df = simulate_ews_survival_data(
        test_subset_df,
        ensemble_specification,
        ensemble_single_incarr,
        null_single_incarr,
        ensemble_single_Reff_thresholds_vec;
        ews_metric = ews_metric,
    )[1]

    detection_survival_vecs, null_survival_vecs = create_ews_survival_data(
        survival_df
    )

    @unpack aggregation = ews_metric_specification

    return ews_survival_plot(
        detection_survival_vecs,
        null_survival_vecs,
        survival_df.enddate;
        ews_aggregation = aggregation,
        burnin = ews_threshold_burnin,
        plottitle = plottitle,
        subtitle = subtitle,
    )
end

function create_ews_survival_data(
    ews_optimal_simulation_df
)
    nsims = nrow(ews_optimal_simulation_df)
    unique_detection_indices = sort(
        filter(
            !isnothing, unique(ews_optimal_simulation_df.detection_index)
        ),
    )

    unique_null_detection_indices = sort(
        filter(
            !isnothing,
            unique(
                ews_optimal_simulation_df.null_detection_index
            ),
        ),
    )

    detection_survival_vec, detection_indices_vec = calculate_detection_indices_and_survival(
        ews_optimal_simulation_df.detection_index,
        unique_detection_indices,
        nsims,
    )

    null_survival_vec, null_indices_vec = calculate_detection_indices_and_survival(
        ews_optimal_simulation_df.null_detection_index,
        unique_null_detection_indices,
        nsims,
    )

    return (
        (; detection_survival_vec, detection_indices_vec),
        (; null_survival_vec, null_indices_vec),
    )
end

function calculate_detection_indices_and_survival(
    detection_vec,
    detection_indices,
    nsims,
)
    if isempty(detection_indices)
        @assert nsims ==
            sum(detection_vec .== nothing)

        survival_vec = [nsims, nsims]
        indices_vec = Int64[]
    else
        detection_indices_counts = map(
            ((i, v),) -> sum(
                detection_vec .== v
            ),
            enumerate(detection_indices),
        )

        @assert nsims ==
            sum(detection_indices_counts) +
                sum(detection_vec .== nothing)

        survival_vec = nsims .- cumsum(detection_indices_counts)

        survival_vec = vcat(
            nsims, nsims, survival_vec..., survival_vec[end]
        )
        detection_indices = Int64.(detection_indices)

        indices_vec = vcat(
            detection_indices[1],
            detection_indices,
        )
    end

    return survival_vec, indices_vec
end

function simulate_ews_survival_data(
    optimal_ews_df,
    ensemble_specification,
    ensemble_single_incarr,
    null_single_incarr,
    thresholds;
    ews_metric = "mean",
    logfilepath = scriptsdir("ensemble-sim_ews-optimization.log.txt"),
)
    logfile = open(logfilepath, "a")

    survival_df = DataFrame(
    (
        ews_metric_specification = EWSMetricSpecification[],
        ews_threshold_burnin = Union{<:Dates.Day,<:Dates.Year}[],
        noise_specification = Union{
            <:PoissonNoiseSpecification,<:DynamicalNoiseSpecification
        }[],
        test_specification = IndividualTestSpecification[],
        ews_metric = String[],
        enddate = Vector{Union{Int64,Try.Err}}(),
        detection_index = Vector{Union{Nothing,Int64}}(),
        null_detection_index = Vector{Union{Nothing,Int64}}(),
    )
)

    subset_optimal_ews_df = subset(
        optimal_ews_df,
        :ews_metric => ByRow(==(ews_metric)),
    )

    for n in [
        :noise_specification,
        :percent_tested,
        :ews_metric_specification,
        :ews_enddate_type,
        :ews_threshold_window,
        :ews_threshold_burnin,
        :ews_metric,
    ]
        @assert length(unique(subset_optimal_ews_df[:, n])) == 1
    end

    noise_specification = subset_optimal_ews_df[1, :noise_specification]
    percent_tested = subset_optimal_ews_df[1, :percent_tested]
    ews_metric_specification = subset_optimal_ews_df[
        1, :ews_metric_specification
    ]
    ews_enddate_type = subset_optimal_ews_df[1, :ews_enddate_type]
    ews_threshold_window = subset_optimal_ews_df[1, :ews_threshold_window]
    ews_threshold_burnin = subset_optimal_ews_df[1, :ews_threshold_burnin]

    ews_metric = subset_optimal_ews_df[1, :ews_metric]

    noisearr = create_noise_arr(
        noise_specification,
        ensemble_single_incarr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )[1]

    enddate_vec = Vector{Union{Try.Ok,Try.Err}}(
        undef, size(ensemble_single_incarr, 3)
    )

    for sim in axes(ensemble_single_incarr, 3)
        enddate_vec[sim] = calculate_ews_enddate(
            thresholds[sim],
            ews_enddate_type,
        )
    end
    failed_sims = sum(Try.iserr.(enddate_vec))

    vec_of_testarr = Vector{Array{Float64,3}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_testarr = Vector{Array{Float64,3}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_ews_vals_vec = Vector{Vector{Union{Missing,EWSMetrics}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_ews_vals_vec = Vector{Vector{Union{Missing,EWSMetrics}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_exceed_thresholds = Vector{Vector{Matrix{Bool}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_exceed_thresholds = Vector{Vector{Matrix{Bool}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_threshold_percentiles = Vector{Vector{Matrix{Float64}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_threshold_percentiles = Vector{Vector{Matrix{Float64}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_detection_index_vec = Vector{Vector{Union{Nothing,Int64}}}(
        undef, nrow(subset_optimal_ews_df)
    )
    vec_of_null_detection_index_vec = Vector{Vector{Union{Nothing,Int64}}}(
        undef, nrow(subset_optimal_ews_df)
    )

    for (i, df_row) in pairs(eachrow(subset_optimal_ews_df))
        test_specification = df_row[:test_specification]
        ews_threshold_percentile = df_row[:ews_threshold_percentile]
        ews_consecutive_thresholds = df_row[:ews_consecutive_thresholds]

        testarr = create_testing_arrs(
            ensemble_single_incarr,
            noisearr,
            percent_tested,
            test_specification,
        )

        vec_of_testarr[i] = testarr

        null_testarr = create_testing_arrs(
            null_single_incarr,
            noisearr,
            percent_tested,
            test_specification,
        )
        vec_of_null_testarr[i] = null_testarr

        ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
            undef, size(testarr, 3)
        )
        null_ews_vals_vec = Vector{Union{Missing,EWSMetrics}}(
            undef, size(testarr, 3)
        )
        fill!(ews_vals_vec, missing)
        fill!(null_ews_vals_vec, missing)

        exceeds_threshold_vec = Vector{Matrix{Bool}}(
            undef, size(testarr, 3)
        )
        null_exceeds_threshold_vec = Vector{Matrix{Bool}}(
            undef, size(testarr, 3)
        )
        threshold_percentiles_vec = Vector{Matrix{Float64}}(
            undef, size(testarr, 3)
        )
        null_threshold_percentiles_vec = Vector{Matrix{Float64}}(
            undef, size(testarr, 3)
        )

        detection_index_vec = Vector{Union{Nothing,Int64}}(
            undef, size(testarr, 3)
        )
        null_detection_index_vec = Vector{Union{Nothing,Int64}}(
            undef, size(testarr, 3)
        )
        fill!(detection_index_vec, nothing)
        fill!(null_detection_index_vec, nothing)

        for sim in axes(testarr, 3)
            if Try.isok(enddate_vec[sim])
                enddate = Try.unwrap(enddate_vec[sim])
                ews_vals_vec[sim] = EWSMetrics(
                    ews_metric_specification,
                    @view(testarr[1:enddate, 5, sim])
                )

                null_ews_vals_vec[sim] = EWSMetrics(
                    ews_metric_specification,
                    @view(null_testarr[1:enddate, 5, sim])
                )

                exceeds_threshold_vec[sim] = expanding_ews_thresholds(
                    ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[2]

                threshold_percentiles_vec[sim] = expanding_ews_thresholds(
                    ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[1]

                detection_index_vec[sim] = calculate_ews_trigger_index(
                    exceeds_threshold_vec[sim];
                    consecutive_thresholds = ews_consecutive_thresholds,
                )

                null_exceeds_threshold_vec[sim] = expanding_ews_thresholds(
                    null_ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[2]

                null_threshold_percentiles_vec[sim] = expanding_ews_thresholds(
                    null_ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    percentiles = ews_threshold_percentile,
                    burn_in = ews_threshold_burnin,
                )[1]

                null_detection_index_vec[sim] = calculate_ews_trigger_index(
                    null_exceeds_threshold_vec[sim];
                    consecutive_thresholds = ews_consecutive_thresholds,
                )

            else
                write(
                    logfile,
                    "$(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))\n",
                )
                write(logfile, "Error:\t$(Try.unwrap_err(enddate_vec[sim]))\n")
                write(logfile, "Simulation:\t$(sim)\n\n")
            end
        end

        vec_of_ews_vals_vec[i] = ews_vals_vec
        vec_of_null_ews_vals_vec[i] = null_ews_vals_vec
        vec_of_exceed_thresholds[i] = exceeds_threshold_vec
        vec_of_null_exceed_thresholds[i] = null_exceeds_threshold_vec
        vec_of_threshold_percentiles[i] = threshold_percentiles_vec
        vec_of_null_threshold_percentiles[i] = null_threshold_percentiles_vec
        vec_of_detection_index_vec[i] = detection_index_vec
        vec_of_null_detection_index_vec[i] = null_detection_index_vec

        append!(
            survival_df,
            DataFrame(;
                ews_metric_specification = repeat(
                    [ews_metric_specification], length(detection_index_vec)
                ),
                ews_threshold_burnin = repeat(
                    [ews_threshold_burnin], length(detection_index_vec)
                ),
                noise_specification = repeat(
                    [noise_specification], length(detection_index_vec)
                ),
                test_specification = repeat(
                    [test_specification], length(detection_index_vec)
                ),
                enddate = Try.unwrap_or_else.(Try.Err, enddate_vec),
                ews_metric = repeat([ews_metric], length(detection_index_vec)),
                detection_index = detection_index_vec,
                null_detection_index = null_detection_index_vec,
            );
            promote = true,
        )
    end

    return survival_df,
    (
        vec_of_testarr,
        vec_of_null_testarr,
        vec_of_ews_vals_vec,
        vec_of_null_ews_vals_vec,
        vec_of_exceed_thresholds,
        vec_of_null_exceed_thresholds,
        vec_of_threshold_percentiles,
        vec_of_null_threshold_percentiles,
        vec_of_detection_index_vec,
        vec_of_null_detection_index_vec,
    ),
    noisearr
end
