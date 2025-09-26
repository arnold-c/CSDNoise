# module DetectionThresholds
#
# export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds

using ProgressMeter
using FLoops
using StatsBase: rle
using UnPack

function create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_specification::OutbreakSpecification
    )
    emergent_incidence_arr = zeros(
        Int64, size(ensemble_inc_vecs, 1), 3, size(ensemble_inc_vecs, 2)
    )

    emergent_outbreak_threshold_vecs = Vector{Array{Int64, 2}}(
        undef, size(ensemble_inc_vecs, 2)
    )

    create_inc_infec_arr!(
        emergent_incidence_arr,
        emergent_outbreak_threshold_vecs,
        ensemble_inc_vecs,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    return emergent_incidence_arr, emergent_outbreak_threshold_vecs
end

function create_inc_infec_arr!(
        emergent_incidence_arr, emergent_outbreak_threshold_vecs, ensemble_inc_vecs,
        outbreakthreshold, minoutbreakdur,
        minoutbreaksize,
    )
    @inbounds for sim in axes(ensemble_inc_vecs, 2)
        convert_svec_to_matrix!(
            @view(emergent_incidence_arr[:, 1, sim]),
            @view(ensemble_inc_vecs[:, sim])
        )

        emergent_incidence_arr[:, 2, sim] .=
            @view(emergent_incidence_arr[:, 1, sim]) .>= outbreakthreshold

        abovethresholdrle = rle(@view(emergent_incidence_arr[:, 2, sim]))

        outbreak_thresholds = calculate_outbreak_thresholds(
            abovethresholdrle; ncols = 4
        )

        emergent_outbreak_threshold_vecs[sim] = classify_all_outbreaks!(
            @view(emergent_incidence_arr[:, 1, sim]),
            @view(emergent_incidence_arr[:, 3, sim]),
            outbreak_thresholds,
            minoutbreakdur,
            minoutbreaksize,
        )
    end
    return nothing
end

function calculate_outbreak_thresholds(outbreakrle; ncols = 4)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])

    outbreak_thresholds = zeros(Int64, length(upperbound_indices), ncols)

    @inbounds outbreak_thresholds[:, 2] .= @view(
        outbreakaccum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1,
        @view(outbreak_thresholds[:, 1]),
        upperbound_indices,
    )

    return outbreak_thresholds
end

function classify_all_outbreaks!(
        incidence_vec,
        alertstatus_vec,
        all_thresholds_arr,
        minoutbreakdur,
        minoutbreaksize,
    )
    for (row, (lower, upper)) in pairs(eachrow(all_thresholds_arr[:, 1:2]))
        all_thresholds_arr[row, 3] = sum(
            @view(incidence_vec[lower:upper])
        )

        all_thresholds_arr[row, 4] = classify_outbreak(
            all_thresholds_arr[row, 3],
            upper,
            lower,
            minoutbreakdur,
            minoutbreaksize,
        )

        @view(alertstatus_vec[lower:upper]) .= @view(
            all_thresholds_arr[row, 4]
        )
    end

    return @view(
        all_thresholds_arr[(all_thresholds_arr[:, 1] .!= 0), :]
    )
end

function calculate_period_sum(incvec)
    return sum(incvec)
end

function classify_outbreak(
        periodsumvec,
        upper_time,
        lower_time,
        minoutbreakdur,
        minoutbreaksize
    )
    if upper_time - lower_time >= minoutbreakdur &&
            periodsumvec >= minoutbreaksize
        return 1
    end
    return 0
end

function Reff_ge_than_one(Reff_vec; ncols = 3)
    @assert ncols >= 3 "ncols must be at least 3"

    Reff_rle = rle(Reff_vec .>= 1)
    Reff_rle_bounds = calculate_outbreak_thresholds(
        Reff_rle; ncols = ncols
    )

    # Amount of time with Reff greater than 1
    @views Reff_rle_bounds[:, 3] .=
        Reff_rle_bounds[:, 2] .- Reff_rle_bounds[:, 1] .+ 1

    return Reff_rle_bounds
end
# end
