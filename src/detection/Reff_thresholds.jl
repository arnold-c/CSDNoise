export calculate_all_Reff_thresholds,
    Reff_ge_than_one

function calculate_all_Reff_thresholds(ensemble_run::EnsembleSEIRRun)
    nsims = length(ensemble_run.emergent_seir_run)

    Reff_vecs = Vector{Thresholds}(undef, nsims)
    for i in eachindex(Reff_vecs)
        Reff_vecs[i] = Reff_ge_than_one(ensemble_run.emergent_seir_run[i].Reff)
    end

    return StructVector(Reff_vecs)
end

function Reff_ge_than_one(Reff_vec)
    Reff_rle = StatsBase.rle(Reff_vec .>= 1)
    return calculate_above_threshold_bounds(Reff_rle)
end
