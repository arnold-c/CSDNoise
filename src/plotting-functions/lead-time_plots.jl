export ews_lead_time_plot

function ews_lead_time_plot(
        lead_time_df;
        week_aggregation = 1,
        ews_method = EWSMethod(Backward()),
        lead_time_units = :days,
    )
    unique_tests = unique(lead_time_df.test_specification)
    lead_time_units_string = string(lead_time_units)

    unique_noise = unique(lead_time_df.noise_magnitude)

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Lead Time By Test and Noise Magnitude",
        subtitle = "$week_aggregation Week Aggregation and $(titlecase(method_string(ews_method))) EWS",
        xlabel = "Noise Magnitude",
        ylabel = "Lead Time ($(lead_time_units))",
        xticks = minimum(unique_noise):maximum(unique_noise),
    )

    for (i, test) in pairs(unique_tests)
        filtered_df = subset(
            lead_time_df,
            :test_specification => DF.ByRow(==(test)),
            :week_aggregation => DF.ByRow(==(week_aggregation)),
            :ews_method => DF.ByRow(==(ews_method)),
            :lead_time_units => DF.ByRow(==(lead_time_units_string)),
        )

        band!(
            ax,
            filtered_df.noise_magnitude,
            filtered_df.lead_time_lower,
            filtered_df.lead_time_upper;
            color = (GLMakie.wong_colors()[i], 0.1),
            # label = get_test_description(test),
        )

        scatterlines!(
            ax,
            filtered_df.noise_magnitude,
            filtered_df.lead_time_median;
            color = GLMakie.wong_colors()[i],
            label = get_test_description(test),
        )
    end

    Legend(fig[1, 2], ax, "Test")

    return fig
end
