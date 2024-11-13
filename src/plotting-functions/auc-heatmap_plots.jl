function tau_auc_heatmap(
    df;
    baseline_test = IndividualTestSpecification(1.0, 1.0, 0),
    colormap = :RdBu,
    textcolorthreshold = 0.6,
    plottitle = "Kendall's Tau AUC Heatmap",
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
            select(df, [:ews_metric, :test_specification, :auc]),
            :test_specification,
            :auc,
        ) |>
        df -> sort(df, order(2; rev = false))

    default_test_metric_order = ordered_df.ews_metric

    mat = Matrix(ordered_df[:, 2:end])'

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
        colorrange = (-1, 1),
    )

    for j in axes(mat, 2), i in axes(mat, 1)
        val = mat[i, j]
        textcolor = abs(val) < textcolorthreshold ? :black : :white
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
