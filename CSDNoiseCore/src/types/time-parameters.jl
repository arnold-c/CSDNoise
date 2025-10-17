export SimTimeParameters

"""
    SimTimeParameters

Parameters defining the temporal structure of a simulation.

# Fields
- `burnin::Float64`: Duration of the burn-in period in days (must be ≤ tmax and > 0)
- `tmin::Float64`: Start time of the simulation
- `tmax::Float64`: End time of the simulation in days (must be > tmin + tstep)
- `tstep::Float64`: Time step size in days (must be > 0)
- `trange::StepRangeLen{Float64, Float64, Float64, Int64}`: Range from tmin to tmax by tstep
- `tspan::Tuple{Float64, Float64}`: Tuple of (tmin, tmax)
- `tlength::Int64`: Number of time steps in trange

# Constructor
    SimTimeParameters(; burnin=365.0*5, tmin=0.0, tmax=365.0*20.0, tstep=1.0)

Creates a `SimTimeParameters` instance with default values suitable for a 20-year simulation
with a 5-year burn-in period and daily time steps.
"""
struct SimTimeParameters
    burnin::Dates.Day
    tmin::Float64
    tmax::Float64
    tstep::Float64
    trange::StepRangeLen{Float64, Float64, Float64, Int64}
    tspan::Tuple{Float64, Float64}
    tlength::Int64
    function SimTimeParameters(
            burnin,
            tmin,
            tmax,
            tstep,
            trange,
            tspan,
            tlength
        )
        @assert Dates.days(burnin) <= tmax
        @assert tmax > tmin + tstep
        for var in (Dates.days(burnin), tstep)
            @assert var > 0.0
        end
        return new(
            burnin,
            tmin,
            tmax,
            tstep,
            trange,
            tspan,
            tlength
        )
    end
end

function SimTimeParameters(;
        burnin = 365.0 * 5,
        tmin = 0.0,
        tmax = 365.0 * 20.0,
        tstep = 1.0
    )
    return SimTimeParameters(
        Dates.Day(burnin),
        tmin,
        tmax,
        tstep,
        tmin:tstep:tmax,
        (tmin, tmax),
        length(tmin:tstep:tmax),
    )
end
