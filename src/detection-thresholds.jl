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

        outbreak_thresholds = calculate_outbreak_thresholds(abovethresholdrle)

        emergent_outbreak_threshold_vecs[sim] = classify_all_outbreaks(
            @view(emergent_incidence_arr[:, 1, sim]),
            @view(emergent_incidence_arr[:, 3, sim]),
            outbreak_thresholds,
            minoutbreakdur,
            minoutbreaksize,
        )
    end
    return nothing
end

function calculate_outbreak_thresholds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])

    upper_bounds = SizedVector{length(upperbound_indices), Int64}(undef)
    lower_bounds = similar(upper_bounds)
    duration = similar(upper_bounds)

    @inbounds upper_bounds .= @view(
        outbreakaccum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1,
        upper_bounds,
        upperbound_indices,
    )
    duration .= upper_bounds .- lower_bounds .+ 1

    return Thresholds(lower_bounds, upper_bounds, duration)
end

function classify_all_outbreaks(
        incidence_vec,
        alertstatus_vec,
        incidence_thresholds::Thresholds,
        minoutbreakdur,
        minoutbreaksize,
    )
    outbreak_status_vec = similar(incidence_thresholds.lower_bounds)
    ninf_during_bounds_vec = similar(incidence_thresholds.lower_bounds)

    for i in eachindex(incidence_thresholds.lower_bounds)
        @inbounds begin
            local lower_bound = incidence_thresholds.lower_bounds[i]
            local upper_bound = incidence_thresholds.upper_bounds[i]
            local duration = incidence_thresholds.duration[i]
            local ninf_during_bounds = sum(
                @view(incidence_vec[lower:upper])
            )
        end

        local oubreak_status = classify_outbreak(
            ninf_during_bounds,
            duration,
            minoutbreakdur,
            minoutbreaksize,
        )

        outbreak_status_vec[i] = outbreak_status
        ninf_during_bounds_vec[i] = ninf_during_bounds

        @view(alertstatus_vec[lower:upper]) .= outbreak_status
    end

    outbreak_idx = findall(==(1), outbreak_status_vec)

    return OutbreakThresholds(
        incidence_thresholds.lower_bounds[outbreak_idx],
        incidence_thresholds.upper_bounds[outbreak_idx],
        incidence_thresholds.duration[outbreak_idx],
        ninf_during_bounds_vec[outbreak_idx]
    )
end

function classify_outbreak(
        ninf_during_bounds,
        duration,
        minoutbreakdur,
        minoutbreaksize
    )
    if duration >= minoutbreakdur && ninf_during_bounds >= minoutbreaksize
        return 1
    end
    return 0
end

function Reff_ge_than_one(Reff_vec)

    Reff_rle = rle(Reff_vec .>= 1)
    return calculate_outbreak_thresholds(Reff_rle)
end
# end
