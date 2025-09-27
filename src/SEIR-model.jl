# module SEIRModel
#
# export calculate_beta_amp, seir_mod, seir_mod!, seir_mod_loop!
#
"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""

using StatsBase
using Distributions: Poisson, Binomial
using Random
using UnPack
using StaticArrays
using LabelledArrays: SLVector

"""
    seir_mod(states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model with a vaccinations going directly to the R compartment and produce the transmission rate array.
"""
function seir_mod(
        states::SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters;
        seed::Int64 = 1234
    )
    tlength = time_params.tlength
    state_vec = FixedSizeVector{typeof(states)}(undef, tlength)
    inc_vec = FixedSizeVector{Int64}(undef, tlength)
    Reff_vec = FixedSizeVector{Float64}(undef, tlength)
    beta_vec = FixedSizeVector{Float64}(undef, tlength)

    for i in eachindex(beta_vec)
        beta_vec[i] = calculate_beta_amp(
            dynamics_params.beta_mean,
            dynamics_params.beta_force,
            time_params.trange[i];
            seasonality = dynamics_params.seasonality
        )
    end

    seir_mod!(
        state_vec,
        inc_vec,
        Reff_vec,
        beta_vec,
        states,
        dynamics_params,
        time_params,
        seed,
    )

    return SEIRRun(
        state_vec,
        inc_vec,
        Reff_vec,
    )
end

function seir_mod(
        states::SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        beta_vec::Vector{Float64},
        time_params::SimTimeParameters;
        seed::Int64 = 1234
    )
    tlength = time_params.tlength
    state_vec = FixedSizeVector{typeof(states)}(undef, tlength)
    inc_vec = FixedSizeVector{Int64}(undef, tlength)
    Reff_vec = FixedSizeVector{Float64}(undef, tlength)

    seir_mod!(
        state_vec,
        inc_vec,
        Reff_vec,
        beta_vec,
        states,
        dynamics_params,
        time_params,
        seed,
    )

    return SEIRRun(
        state_vec,
        inc_vec,
        Reff_vec,
    )
end


"""
    seir_mod!(state_arr, change_arr, jump_arr, beta_arr, states, dynamics_params, trange; tstep, type = "stoch")

The in-place function to run the SEIR model and produce the transmission rate array.
"""
function seir_mod!(
        state_vec::ASV,
        inc_vec::AI,
        Reff_vec::AF,
        beta_vec::Vector{Float64},
        states::SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters,
        seed::Int64,
    ) where {
        ASV <: AbstractVector{SVector{5, Int64}},
        AI <: AbstractVector{Int64},
        AF <: AbstractVector{Float64},
    }
    Random.seed!(seed)

    @inbounds begin
        mu = dynamics_params.mu
        epsilon = dynamics_params.epsilon
        sigma = dynamics_params.sigma
        gamma = dynamics_params.gamma
        R_0 = dynamics_params.R_0
        timestep = time_params.tstep
        tlength = time_params.tlength
        burnin_days = time_params.burnin

        state_vec[1] = states
        inc_vec[1] = 0
        Reff_vec[1] = calculateReffective(
            beta_vec[1],
            dynamics_params,
            states[1],
            states[5],
        )
    end

    # Pre-compute constant values to avoid repeated calculations
    mu_timestep = mu * timestep
    sigma_timestep = sigma * timestep
    gamma_timestep = gamma * timestep
    epsilon_over_R0_timestep = (epsilon / R_0) * timestep

    @inbounds for i in 2:(tlength)
        vaccination_coverage = if i < burnin_days
            dynamics_params.burnin_vaccination_coverage
        else
            dynamics_params.vaccination_coverage
        end

        state_vec[i], inc_vec[i] = seir_mod_loop(
            state_vec[i - 1],
            beta_vec[i - 1],
            mu_timestep,
            epsilon_over_R0_timestep,
            sigma_timestep,
            gamma_timestep,
            vaccination_coverage,
            timestep,
        )
        Reff_vec[i] = calculateReffective(
            beta_vec[i],
            dynamics_params,
            state_vec[i][1],
            state_vec[i][5],
        )
    end

    return nothing
end

"""
    seir_mod_loop!(state_arr, change_arr, jump_arr, j, params, t, tstep; type = type)

The inner loop that is called by `seir_mod!()` function.
"""
@inline function seir_mod_loop(
        state_vec::SVector{5, Int64},
        beta_t::Float64,
        mu_timestep::Float64,
        epsilon_over_R0_timestep::Float64,
        sigma_timestep::Float64,
        gamma_timestep::Float64,
        vaccination_coverage::Float64,
        timestep::Float64,
    )::Tuple{SVector{5, Int64}, Int64}

    @inbounds begin
        S = state_vec[1]
        E = state_vec[2]
        I = state_vec[3]
        R = state_vec[4]
        N = state_vec[5]

        # Pre-compute common terms
        beta_S_I_timestep = beta_t * S * I * timestep
        mu_N_timestep = mu_timestep * N

        contact_inf = rand(Poisson(beta_S_I_timestep)) # Contact: S -> E
        S_births = rand(Poisson(mu_N_timestep * (1 - vaccination_coverage))) # Birth -> S
        S_death = rand(Poisson(mu_timestep * S)) # S -> death
        R_death = rand(Poisson(mu_timestep * R)) # R -> death
        import_inf = rand(Poisson(epsilon_over_R0_timestep * N)) # Import: S -> E
        R_births = rand(Poisson(mu_N_timestep * vaccination_coverage)) # Birth -> R
        latent = rand(Binomial(E, sigma_timestep)) # E -> I
        E_death = rand(Binomial(E - latent, mu_timestep)) # E -> death
        recovery = rand(Binomial(I, gamma_timestep)) # I -> R
        I_death = rand(Binomial(I - recovery, mu_timestep)) # I -> death

        dS = S_births - (contact_inf + import_inf + S_death)
        dE = (contact_inf + import_inf) - (latent + E_death)
        dI = latent - (recovery + I_death)
        dR = (recovery + R_births) - R_death
        dN = dS + dE + dI + dR
    end

    return (
        SVector(S + dS, E + dE, I + dI, R + dR, N + dN),
        contact_inf,
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
