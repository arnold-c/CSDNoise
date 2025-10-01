using Dates: Dates
using UnPack: @unpack

function prepare_survival_facet_params(
        detection_survival_vecs,
        null_survival_vecs,
        enddate_vec;
        ews_aggregation = Dates.Day(7),
        endpoint_aggregation = Dates.Day(30),
    )
    @unpack detection_survival_vec, detection_indices_vec =
        detection_survival_vecs
    @unpack null_survival_vec, null_indices_vec = null_survival_vecs

    filtered_enddate_vec = filter(isinteger, enddate_vec)

    times =
        collect(1:Dates.days(ews_aggregation):maximum(filtered_enddate_vec)) ./
        365

    detection_survival_times = vcat(
        0,
        times[detection_indices_vec]...,
        times[end],
    )
    null_survival_times = vcat(
        0,
        times[null_indices_vec]...,
        times[end],
    )

    nsims = detection_survival_vec[1]

    enddate_vec = div.(enddate_vec, Dates.days(endpoint_aggregation))
    unique_enddate_vec = sort(unique(enddate_vec))
    enddate_times =
        (unique_enddate_vec .* Dates.days(endpoint_aggregation)) ./ 365

    enddate_counts = map(unique_enddate_vec) do enddate
        sum(enddate_vec .== enddate)
    end

    return times, enddate_times,
        enddate_counts,
        detection_survival_times,
        detection_survival_vec,
        null_survival_times,
        null_survival_vec,
        nsims
end
