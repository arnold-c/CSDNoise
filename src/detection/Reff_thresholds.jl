export Reff_ge_than_one

function Reff_ge_than_one(Reff_vec)
    Reff_rle = StatsBase.rle(Reff_vec .>= 1)
    return calculate_above_threshold_bounds(Reff_rle)
end
