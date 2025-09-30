using GLMakie
using DataFrames
using Printf

export optimal_ews_heatmap_df,
    optimal_ews_heatmap_plot

function optimal_ews_heatmap_df(
        optimal_ews_df;
        tiebreaker_preference = "speed",
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
    tiebreaker_args = tiebreaker_args = if tiebreaker_preference == "speed"
        (:ews_consecutive_thresholds, false)
    elseif tiebreaker_preference == "specificity"
        (:specificity, true)
    else
        error(
            "Invalid preference: $tiebreaker_preference. Please choose either \"speed\" or \"specificity\"."
        )
    end

    return map(
        collect(
            groupby(
                optimal_ews_df,
                optimal_grouping_parameters,
            ),
        ),
    ) do df
        sort(df, order(tiebreaker_args[1]; rev = tiebreaker_args[2]))[1, :]
    end |>
        x -> vcat(DataFrame.(x)...; cols = :union)
end

function optimal_ews_heatmap_plot(
        df;
        outcome = :accuracy,
        baseline_test = IndividualTestSpecification(1.0, 1.0, 0),
        colormap = :RdBu,
        colorrange = [0.2, 0.8],
        textcolorthreshold = (0.4, 0.68),
        xlabel = "Test Sensitivity & Specificity",
        ylabel = "EWS Metric",
        colorlabel = "Accuracy",
        accuracy_fontsize = 12,
        rest_fontsize = 12,
        legendsize = 22,
        xlabelsize = 22,
        ylabelsize = 22,
        xticklabelsize = 22,
        yticklabelsize = 22,
        legendticklabelsize = 22,
        legendwidth = 20,
        kwargs...,
    )
    kwargs_dict = Dict{Symbol, Any}(kwargs)

    outcome_str = titlecase(string(outcome))

    if !haskey(kwargs_dict, :plottitle)
        plottitle =
            "Heatmap: " * titlecase(string(outcome))
    else
        plottitle = kwargs_dict[:plottitle]
    end

    if !haskey(kwargs_dict, :subtitle)
        ews_enddate_type_str = split(string(df[1, :ews_enddate_type]), "::")[1]

        subtitle =
            "Noise: $(get_noise_magnitude_description(df[1, :noise_specification])), $(ews_enddate_type_str)" *
            "\nP = Quantile Threshold, C = Consecutive Thresholds, S = Specificity"
    else
        subtitle = kwargs_dict[:subtitle]
    end

    df[!, :test_sens] =
        getproperty.(df.test_specification, :sensitivity)
    df[!, :test_spec] = getproperty.(df.test_specification, :sensitivity)
    df[!, :test_result_lag] =
        getproperty.(df.test_specification, :test_result_lag)
    sort!(
        df,
        [:test_sens, :test_spec, :test_result_lag];
        rev = [true, true, false],
    )

    unique_tests = unique(df.test_specification)
    if mapreduce(x -> x == baseline_test, +, unique_tests) == 0
        error("default_test $baseline_test must be in unique_tests")
    end

    default_test_metric_order, mat = _create_ews_heatmap_matrix(
        df,
        outcome,
    )
    default_test_metric_order_labels = clean_ews_metric_names(
        default_test_metric_order
    )

    if minimum(mat) < colorrange[1]
        colorrange[1] = floor(minimum(mat); digits = 1)
    end
    if maximum(mat) > colorrange[2]
        colorrange[2] = ceil(maximum(mat); digits = 1)
    end

    threshold_quantile_matrix = _create_ews_heatmap_matrix(
        df,
        :ews_threshold_quantile,
        default_test_metric_order,
    )

    consecutive_thresholds_df = _create_ews_heatmap_matrix(
        df,
        :ews_consecutive_thresholds,
        default_test_metric_order,
    )

    specificity_df = _create_ews_heatmap_matrix(
        df,
        :specificity,
        default_test_metric_order,
    )

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        subtitle = subtitle,
        xlabel = xlabel,
        ylabel = ylabel,
        xlabelsize = xlabelsize,
        ylabelsize = ylabelsize,
        xticks = (1:length(unique_tests), test_axis_label.(unique_tests)),
        yticks = (
            1:length(default_test_metric_order),
            default_test_metric_order_labels,
        ),
        xticklabelsize = xticklabelsize,
        yticklabelsize = yticklabelsize,
    )

    hmap = heatmap!(
        ax,
        mat;
        colormap = colormap,
        colorrange = colorrange,
    )

    for j in axes(mat, 2), i in axes(mat, 1)
        val = mat[i, j]
        if length(textcolorthreshold) == 1
            textcolor = val <= textcolorthreshold ? :black : :white
        elseif length(textcolorthreshold) == 2
            textcolor =
            if val >= textcolorthreshold[1] && val <= textcolorthreshold[2]
                :black
            else
                :white
            end
        else
            error(
                "variable `textcolorthreshold` should be length 1 or 2. Instead received length $(length(textcolorthreshold))"
            )
        end
        acc = @sprintf("%.2f", mat[i, j])
        perc = @sprintf("%.2f", threshold_quantile_matrix[i, j])
        spec = @sprintf("%.2f", specificity_df[i, j])
        consec = consecutive_thresholds_df[i, j]
        label = rich("accuracy")
        acc_label = "Acc: $(acc)\n"
        acc_label_length = length(acc_label)
        rest_label = "Q: $(perc) C: $(consec)\nS: $(spec)"
        rest_label_length = length(rest_label)
        label = acc_label * rest_label
        text!(
            ax,
            label;
            position = (i, j),
            color = textcolor,
            align = (:center, :center),
            fontsize = [
                fill(accuracy_fontsize, acc_label_length);
                fill(rest_fontsize, rest_label_length)
            ],
            font = [
                fill("TeX Gyre Heros Makie Bold", acc_label_length);
                fill("TeX Gyre Heros Makie", rest_label_length)
            ],
        )
    end

    limits!(
        ax,
        (0, length(unique_tests) + 1),
        (0, length(default_test_metric_order) + 1),
    )
    #
    Colorbar(
        fig[1, 2],
        hmap;
        label = colorlabel,
        width = legendwidth,
        labelsize = legendsize,
        ticklabelsize = legendticklabelsize,
    )
    return fig
end

function _create_ews_heatmap_matrix(
        df, outcome::Symbol, ews_metric_order
    )
    ordered_df =
        unstack(
        select(df, [:ews_metric, :test_specification, outcome]),
        :test_specification,
        outcome,
    ) |>
        df -> df[indexin(ews_metric_order, df.ews_metric), :]

    return Matrix(ordered_df[:, 2:end])'
end

function _create_ews_heatmap_matrix(
        df, outcome::Symbol
    )
    ordered_df =
        unstack(
        select(df, [:ews_metric, :test_specification, outcome]),
        :test_specification,
        outcome,
    ) |>
        df -> sort(df, order(2; rev = false))

    default_test_metric_order = ordered_df.ews_metric

    return default_test_metric_order, Matrix(ordered_df[:, 2:end])'
end
