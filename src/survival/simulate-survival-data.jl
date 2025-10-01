export simulate_ews_survival_data

function simulate_ews_survival_data(
        optimal_ews_df,
        ensemble_specification,
        emergent_incidence_arr,
        null_incidence_arr,
        thresholds;
        ews_metric = "mean",
        logfilename = "ensemble-sim_ews-optimization.log.txt",
        logfiledir = DrWatson.scriptsdir()
    )
    !isdir(logfiledir) && mkpath(logfiledir)

    logfilepath = joinpath(logfiledir, logfilename)
    logfile = open(logfilepath, "a")

    survival_df = DF.DataFrame(
        (
            ews_metric_specification = EWSMetricSpecification[],
            ews_threshold_burnin = Union{<:Dates.Day, <:Dates.Year}[],
            noise_specification = Union{
                <:PoissonNoiseSpecification, <:DynamicalNoiseSpecification,
            }[],
            test_specification = IndividualTestSpecification[],
            ews_metric = String[],
            enddate = Vector{Union{Int64, Try.Err}}(),
            detection_index = Vector{Union{Nothing, Int64}}(),
            null_detection_index = Vector{Union{Nothing, Int64}}(),
        )
    )

    subset_optimal_ews_df = DF.subset(
        optimal_ews_df,
        :ews_metric => DF.ByRow(==(ews_metric)),
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
        1, :ews_metric_specification,
    ]
    ews_enddate_type = subset_optimal_ews_df[1, :ews_enddate_type]
    ews_threshold_window = subset_optimal_ews_df[1, :ews_threshold_window]
    ews_threshold_burnin = subset_optimal_ews_df[1, :ews_threshold_burnin]

    ews_metric = subset_optimal_ews_df[1, :ews_metric]

    noisearr = create_noise_arr(
        noise_specification,
        emergent_incidence_arr;
        ensemble_specification = ensemble_specification,
        seed = 1234,
    )[1]

    enddate_vec = Vector{Union{Try.Ok, Try.Err}}(
        undef, size(emergent_incidence_arr, 3)
    )

    for sim in axes(emergent_incidence_arr, 3)
        enddate_vec[sim] = calculate_ews_enddate(
            thresholds[sim],
            ews_enddate_type,
        )
    end
    failed_sims = sum(Try.iserr.(enddate_vec))

    vec_of_testarr = Vector{Array{Float64, 3}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_null_testarr = Vector{Array{Float64, 3}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_ews_vals_vec = Vector{Vector{Union{Missing, EWSMetrics}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_null_ews_vals_vec = Vector{Vector{Union{Missing, EWSMetrics}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_exceed_thresholds = Vector{Vector{Matrix{Bool}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_null_exceed_thresholds = Vector{Vector{Matrix{Bool}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_threshold_quantiles = Vector{Vector{Matrix{Float64}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_null_threshold_quantiles = Vector{Vector{Matrix{Float64}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_detection_index_vec = Vector{Vector{Union{Nothing, Int64}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )
    vec_of_null_detection_index_vec = Vector{Vector{Union{Nothing, Int64}}}(
        undef, DF.nrow(subset_optimal_ews_df)
    )

    for (i, df_row) in pairs(eachrow(subset_optimal_ews_df))
        test_specification = df_row[:test_specification]
        ews_threshold_quantile = df_row[:ews_threshold_quantile]
        ews_consecutive_thresholds = df_row[:ews_consecutive_thresholds]

        testarr = create_testing_arrs(
            emergent_incidence_arr,
            noisearr,
            percent_tested,
            test_specification,
        )

        vec_of_testarr[i] = testarr

        null_testarr = create_testing_arrs(
            null_incidence_arr,
            noisearr,
            percent_tested,
            test_specification,
        )
        vec_of_null_testarr[i] = null_testarr

        ews_vals_vec = Vector{Union{Missing, EWSMetrics}}(
            undef, size(testarr, 3)
        )
        null_ews_vals_vec = Vector{Union{Missing, EWSMetrics}}(
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
        threshold_quantiles_vec = Vector{Matrix{Float64}}(
            undef, size(testarr, 3)
        )
        null_threshold_quantiles_vec = Vector{Matrix{Float64}}(
            undef, size(testarr, 3)
        )

        detection_index_vec = Vector{Union{Nothing, Int64}}(
            undef, size(testarr, 3)
        )
        null_detection_index_vec = Vector{Union{Nothing, Int64}}(
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
                    quantiles = ews_threshold_quantile,
                    burn_in = ews_threshold_burnin,
                )[2]

                threshold_quantiles_vec[sim] = expanding_ews_thresholds(
                    ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    quantiles = ews_threshold_quantile,
                    burn_in = ews_threshold_burnin,
                )[1]

                detection_index_vec[sim] = Try.@? calculate_ews_trigger_index(
                    exceeds_threshold_vec[sim],
                    ews_consecutive_thresholds,
                )

                null_exceeds_threshold_vec[sim] = expanding_ews_thresholds(
                    null_ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    quantiles = ews_threshold_quantile,
                    burn_in = ews_threshold_burnin,
                )[2]

                null_threshold_quantiles_vec[sim] = expanding_ews_thresholds(
                    null_ews_vals_vec[sim],
                    Symbol(ews_metric),
                    ews_threshold_window;
                    quantiles = ews_threshold_quantile,
                    burn_in = ews_threshold_burnin,
                )[1]

                null_detection_index_vec[sim] = Try.@? calculate_ews_trigger_index(
                    null_exceeds_threshold_vec[sim],
                    ews_consecutive_thresholds,
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
        vec_of_threshold_quantiles[i] = threshold_quantiles_vec
        vec_of_null_threshold_quantiles[i] = null_threshold_quantiles_vec
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
            vec_of_threshold_quantiles,
            vec_of_null_threshold_quantiles,
            vec_of_detection_index_vec,
            vec_of_null_detection_index_vec,
        ),
        noisearr
end
