# module SEIRModel
#
# export calculate_beta_amp, seir_mod, seir_mod!, seir_mod_loop!
#
"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""

using DifferentialEquations
using Statistics
using Distributions
using Random
using UnPack
using LoopVectorization
using StaticArrays

"""
    seir_mod(states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model with a vaccinations going directly to the R compartment and produce the transmission rate array.
"""
function seir_mod(
    states, dynamics_params, time_params; seed = 1234
)
    state_vec = Vector{typeof(states)}(undef, time_params.tlength)
    beta_vec = Vector{Float64}(undef, time_params.tlength)
    inc_vec = Vector{typeof(SVector(0))}(undef, time_params.tlength)

    seir_mod!(
        state_vec,
        inc_vec,
        beta_vec,
        states,
        dynamics_params,
        time_params;
        seed = seed,
    )

    return state_vec, inc_vec, beta_vec
end

"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_mod!(
    state_vec,
    inc_vec,
    beta_vec,
    states,
    dynamics_params,
    time_params;
    seed = 1234,
)
    Random.seed!(seed)

    @inbounds begin
        mu = dynamics_params.mu
        epsilon = dynamics_params.epsilon
        sigma = dynamics_params.sigma
        gamma = dynamics_params.gamma
        R_0 = dynamics_params.R_0
        vaccination_coverage = dynamics_params.vaccination_coverage
        timestep = time_params.tstep
        beta_mean = dynamics_params.beta_mean
        beta_force = dynamics_params.beta_force
        trange = time_params.trange

        state_vec[1] = states
        inc_vec[1] = SVector(0)
    end

    @. beta_vec = calculate_beta_amp(
        beta_mean, beta_force, trange; seasonality = dynamics_params.seasonality
    )

    @inbounds for i in 2:(time_params.tlength)
        state_vec[i], inc_vec[i] = seir_mod_loop!(
            state_vec[i - 1],
            beta_vec[i],
            mu,
            epsilon,
            sigma,
            gamma,
            R_0,
            vaccination_coverage,
            timestep,
        )
    end

    return nothing
end

"""
    seir_mod_loop!(state_arr, change_arr, jump_arr, j, params, t, tstep; type = type)

The inner loop that is called by `seir_mod!()` function.
"""
function seir_mod_loop!(
    state_vec,
    beta_t,
    mu,
    epsilon,
    sigma,
    gamma,
    R_0,
    vaccination_coverage,
    timestep,
)

    # TODO: Benchmak StaticArrays implementation as potentially much faster.
    # Would need to use permutedims(reshape(reinterperate(Float64, SVector), (...), (...))
    # to get it into an array that could be used later on.
    # Create views of the state variables for easier use
    @inbounds begin
        S = state_vec[1]
        E = state_vec[2]
        I = state_vec[3]
        R = state_vec[4]
        N = state_vec[5]

        contact_inf = rand(Poisson(beta_t * S * I * timestep)) # Contact: S -> E
        S_births = rand(Poisson(mu * (1 - vaccination_coverage) * N * timestep)) # Birth -> S
        S_death = rand(Poisson(mu * S * timestep)) # S -> death
        R_death = rand(Poisson(mu * R * timestep)) # R -> death
        import_inf = rand(Poisson((epsilon * N / R_0) * timestep)) # Import: S -> E
        R_births = rand(Poisson(mu * vaccination_coverage * N * timestep)) # Birth -> R
        latent = rand(Binomial(E, sigma * timestep)) # E -> I
        E_death = rand(Binomial(E - latent, mu * timestep)) # E -> death
        recovery = rand(Binomial(I, gamma * timestep)) # I -> R
        I_death = rand(Binomial(I - recovery, mu * timestep)) # I -> death

        dS = S_births - (contact_inf + import_inf + S_death)
        dE = (contact_inf + import_inf) - (latent + E_death)
        dI = latent - (recovery + I_death)
        dR = (recovery + R_births) - R_death
        dN = dS + dE + dI + dR
    end

    return (
        SVector(S + dS, E + dE, I + dI, R + dR, N + dN), SVector(contact_inf)
    )
end

function convert_svec_to_matrix(svec)
    arr = Matrix{Int64}(undef, size(svec, 1), length(svec[1]))
    convert_svec_to_matrix!(arr, svec)
    return arr
end
function convert_svec_to_matrix!(arr, svec)
    @inbounds for state in eachindex(svec[1]), time in axes(svec, 1)
        arr[time, state] = svec[time][state]
    end
    return nothing
end

function convert_svec_to_array(svec)
    arr = Array{Int64}(undef, size(svec, 1), length(svec[1]), size(svec, 2))
    @inbounds for sim in axes(svec, 2)
        convert_svec_to_matrix!(@view(arr[:, :, sim]), @view(svec[:, sim]))
    end
    return arr
end

# end
