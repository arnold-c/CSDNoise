# module DetectionThresholds
#
# export create_inc_infec_arr, create_inc_infec_arr!, calculate_outbreak_thresholds

using ProgressMeter
using StatsBase: rle
using UnPack
using FixedSizeArrays: FixedSizeVector
using StructArrays: StructVector

function calculate_outbreak_thresholds(
        seir_results::StructVector{SEIRRun},
        outbreak_specification::OutbreakSpecification
    )
    return calculate_outbreak_thresholds(
        seir_results.incidence,
        outbreak_specification
    )
end

function calculate_outbreak_thresholds(
        incidence_vecs::VSV,
        outbreak_specification::OutbreakSpecification
    ) where {VSV <: AbstractVector{<:AbstractVector}}

    emergent_outbreak_threshold_vecs = Vector{OutbreakThresholds}(
        undef, length(incidence_vecs)
    )

    calculate_outbreak_thresholds!(
        emergent_outbreak_threshold_vecs,
        incidence_vecs,
        outbreak_specification.outbreak_threshold,
        outbreak_specification.minimum_outbreak_duration,
        outbreak_specification.minimum_outbreak_size,
    )

    return StructVector(emergent_outbreak_threshold_vecs)
end

function calculate_outbreak_thresholds!(
        emergent_outbreak_threshold_vecs,
        incidence_vecs,
        outbreakthreshold,
        minoutbreakdur,
        minoutbreaksize,
    )
    above_threshold_worker_vec = FixedSizeVector{Bool}(undef, length(incidence_vecs[1]))
    @inbounds for (sim, incidence_vec) in pairs(incidence_vecs)
        above_threshold_worker_vec .= incidence_vec .>= outbreakthreshold

        abovethresholdrle = rle(above_threshold_worker_vec)

        threshold_bounds = calculate_above_threshold_bounds(abovethresholdrle)

        local outbreak_thresholds = classify_all_outbreaks(
            incidence_vec,
            threshold_bounds,
            minoutbreakdur,
            minoutbreaksize,
        )
        emergent_outbreak_threshold_vecs[sim] = outbreak_thresholds
    end
    return nothing
end

function calculate_above_threshold_bounds(outbreakrle)
    # Calculate upper and lower indices of consecutive days of infection
    outbreakaccum = accumulate(+, outbreakrle[2])
    upperbound_indices = findall(isequal(1), outbreakrle[1])

    upper_bounds = FixedSizeVector{Int64}(undef, length(upperbound_indices))
    lower_bounds = similar(upper_bounds)
    duration = similar(upper_bounds)

    @inbounds upper_bounds .= @view(
        outbreakaccum[upperbound_indices]
    )
    map!(
        x -> x - 1 == 0 ? 1 : outbreakaccum[x - 1] + 1,
        lower_bounds,
        upperbound_indices,
    )
    duration .= upper_bounds .- lower_bounds .+ 1

    return Thresholds(lower_bounds, upper_bounds, duration)
end

function classify_all_outbreaks(
        incidence_vec,
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
                @view(incidence_vec[lower_bound:upper_bound])
            )
        end

        local outbreak_status = classify_outbreak(
            ninf_during_bounds,
            duration,
            minoutbreakdur,
            minoutbreaksize,
        )

        outbreak_status_vec[i] = outbreak_status
        ninf_during_bounds_vec[i] = ninf_during_bounds

    end

    outbreak_idx = findall(isequal(1), outbreak_status_vec)
    num_outbreaks = length(outbreak_idx)

    return OutbreakThresholds(
        FixedSizeVector{Int64}(incidence_thresholds.lower_bounds[outbreak_idx]),
        FixedSizeVector{Int64}(incidence_thresholds.upper_bounds[outbreak_idx]),
        FixedSizeVector{Int64}(incidence_thresholds.duration[outbreak_idx]),
        FixedSizeVector{Int64}(ninf_during_bounds_vec[outbreak_idx])
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
    return calculate_above_threshold_bounds(Reff_rle)
end
# end
