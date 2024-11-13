function tau_auc_heatmap(
    df,
    outcome = :auc_magnitude;
    baseline_test = IndividualTestSpecification(1.0, 1.0, 0),
    plottitle = "Kendall's Tau $(replace(titlecase(replace(string(outcome), "_" => " ")), "Auc" => "AUC")) Heatmap",
    colormap = :RdBu,
    colorrange = [0.2, 0.8],
    textcolorthreshold = (0.4, 0.68),
    kwargs...,
)
    df[!, :test_sens] = getproperty.(df.test_specification, :sensitivity)
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
        error(
            "default_test $baseline_test must be in unique_tests. Unique tests are $unique_tests"
        )
    end

    ordered_df =
        unstack(
            select(df, [:ews_metric, :test_specification, outcome]),
            :test_specification,
            outcome,
        ) |>
        df -> sort(df, order(2; rev = false))

    default_test_metric_order = ordered_df.ews_metric

    mat = Matrix(ordered_df[:, 2:end])'

    if minimum(mat) < colorrange[1]
        colorrange[1] = floor(minimum(mat); digits = 1)
    end
    if maximum(mat) > colorrange[2]
        colorrange[2] = ceil(maximum(mat); digits = 1)
    end

    function test_axis_label(test)
        return "($(test.sensitivity), $(test.specificity), $(test.test_result_lag))"
    end

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = plottitle,
        xlabel = "Test Specification (Sensitivity, Specificity, Lag)",
        ylabel = "EWS Metric",
        xticks = (1:length(unique_tests), test_axis_label.(unique_tests)),
        yticks = (
            1:length(default_test_metric_order), default_test_metric_order
        ),
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
        text!(
            ax,
            "$(round(mat[i,j], digits = 5))";
            position = (i, j),
            color = textcolor,
            align = (:center, :center),
        )
    end

    limits!(
        ax,
        (0, length(unique_tests) + 1),
        (0, length(default_test_metric_order) + 1),
    )
    Colorbar(fig[1, 2], hmap; label = "values", width = 15, ticksize = 15)
    return fig
end
