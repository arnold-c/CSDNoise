# module SEIRModel
#
# export calculate_beta_amp, seir_mod, seir_mod!, seir_mod_loop!
#
"""
This is a simulation of an SIR model that uses Tau-leaping, with commuter
imports. All jumps are manually defined.
"""

using StatsBase
using Random: Random
using Distributions: Poisson, Binomial
using StaticArrays

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
    Random.seed!(seed)

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
    Random.seed!(seed)

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
        Reff_vec::AF1,
        beta_vec::AF2,
        states::SVector{5, Int64},
        dynamics_params::DynamicsParameters,
        time_params::SimTimeParameters,
    ) where {
        ASV <: AbstractVector{SVector{5, Int64}},
        AI <: AbstractVector{Int64},
        AF1 <: AbstractVector{Float64},
        AF2 <: AbstractVector{Float64},
    }

    @inbounds begin
        mu = dynamics_params.mu
        epsilon = dynamics_params.epsilon
        sigma = dynamics_params.sigma
        gamma = dynamics_params.gamma
        R_0 = dynamics_params.R_0
        timestep = time_params.tstep
        tlength = length(inc_vec)
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
    epsilon_timestep = epsilon * timestep

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
            epsilon_timestep,
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
@noinline function seir_mod_loop(
        state_vec::SVector{5, Int64},
        beta_t::Float64,
        mu_timestep::Float64,
        epsilon_timestep::Float64,
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

        # Pre-compute common terms and probabilities (reuse these!)
        mu_N_timestep = mu_timestep * N

        # Pre-compute all probabilities once
        contact_prob = convert_rate_to_prob(beta_t * I * timestep)
        import_prob = convert_rate_to_prob(epsilon_timestep)
        death_prob = convert_rate_to_prob(mu_timestep)  # Reuse for all death processes
        latent_prob = convert_rate_to_prob(sigma_timestep)
        recovery_prob = convert_rate_to_prob(gamma_timestep)

        # Births (Poisson since external arrivals)
        S_births = rand(Poisson(mu_N_timestep * (1 - vaccination_coverage)))
        R_births = rand(Poisson(mu_N_timestep * vaccination_coverage))

        # For S compartment: handle competing outflows sequentially
        remaining_S = S

        # Contact infections: S -> E (use smart transition)
        contact_rate = remaining_S * contact_prob
        contact_inf = remaining_S > 0 ?
            smart_transition(remaining_S, contact_prob, contact_rate) :
            0
        remaining_S = max(0, remaining_S - contact_inf)

        # Import infections: S -> E (use smart transition)
        import_rate = remaining_S * import_prob
        import_inf = remaining_S > 0 ?
            smart_transition(remaining_S, import_prob, import_rate) :
            0
        remaining_S = max(0, remaining_S - import_inf)

        # S deaths (use smart transition with pre-computed death_prob)
        S_death_rate = remaining_S * death_prob
        S_death = remaining_S > 0 ?
            smart_transition(remaining_S, death_prob, S_death_rate) :
            0

        # E compartment transitions
        latent_rate = E * latent_prob
        latent = smart_transition(E, latent_prob, latent_rate)

        remaining_E = E - latent
        E_death_rate = remaining_E * death_prob  # Reuse death_prob
        E_death = smart_transition(remaining_E, death_prob, E_death_rate)

        # I compartment transitions
        recovery_rate = I * recovery_prob
        recovery = smart_transition(I, recovery_prob, recovery_rate)

        remaining_I = I - recovery
        I_death_rate = remaining_I * death_prob  # Reuse death_prob
        I_death = smart_transition(remaining_I, death_prob, I_death_rate)

        # R deaths (use smart transition with pre-computed death_prob)
        R_death_rate = R * death_prob  # Reuse death_prob
        R_death = smart_transition(R, death_prob, R_death_rate)

        # Calculate changes
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

@inline function convert_rate_to_prob(rate)
    return min(rate, 1.0)
end

@inline function smart_transition(n::Int64, prob::Float64)
    rate = n * prob
    return smart_transition(n, prob, rate)
end

# Smart transition function that chooses between Poisson and Binomial
@inline function smart_transition(n::Int64, prob::Float64, rate::Float64)
    # Use Poisson when:
    # - Population is large (n > 30)
    # - Expected value is small (rate < 10)
    # This gives good approximation with better performance

    if n == 0
        return 0
    end

    if n > 30 && rate < 10.0
        # Use Poisson approximation for large n, small rate
        return min(rand(Poisson(rate)), n)  # Cap at population size
    else
        # Use exact Binomial for small n or large rate
        return rand(Binomial(n, prob))
    end
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
