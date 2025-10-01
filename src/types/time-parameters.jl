export SimTimeParameters

Base.@kwdef struct SimTimeParameters
    burnin::Float64
    tmin::Float64
    tmax::Float64
    tstep::Float64
    trange::StepRangeLen{Float64, Float64, Float64, Int64}
    tspan::Tuple{Float64, Float64}
    tlength::Int64
end

function SimTimeParameters(;
        burnin = 0.0, tmin = 0.0, tmax = 365.0 * 100.0, tstep = 1.0
    )
    @assert burnin <= tmax
    return SimTimeParameters(
        burnin, tmin, tmax, tstep, tmin:tstep:tmax, (tmin, tmax),
        length(tmin:tstep:tmax),
    )
end
