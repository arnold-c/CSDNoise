# module TransmissionFunctions

using LinearAlgebra
using LightSumTypes: @sumtype, variant

# Define seasonality function variants
struct CosineSeasonality end
struct SineSeasonality end

# Create sum type for seasonality functions
@sumtype SeasonalityFunction(CosineSeasonality, SineSeasonality)

# export calculate_beta, calculateR0, calculate_import_rate

"""
    calculate_beta(
    	R_0::AbstractFloat,
    	gamma::AbstractFloat,
    	mu::AbstractFloat,
		contact_mat::Matrix{AbstractFloat},
		pop_matrix::Vector{AbstractFloat}
    )

Calculate the value beta for a given set of parameters, a contact matrix, and the population matrix, using the Next-Generation Matrix (K).
The contact matrix represents the contact rates between an individual in group 𝒾 and an individual in group 𝒿.
The population matrix is a vector that represents the number of individual in each group.

Uses the property of eigenvalues that for any scalar c and matrix A, the spectral radius ρ(c⋅A)=∣c∣⋅ρ(A).
Applying this to our equation for R_0:

R_0 = ρ(K) where K = F⋅V⁻¹ = β⋅K′, K′ = K / β = Q⋅V⁻¹, Q = contact_mat * pop_matrix

R_0 = ρ(β K′)

R_0 = β ρ(K′)

β = R_0 / ρ(Q⋅V⁻¹)
"""
function calculate_beta(
        R_0::T, gamma::T, mu::T, contact_mat::Matrix{T}, pop_matrix::Vector{T}
    ) where {T <: AbstractFloat}
    # TODO: Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.
    if size(contact_mat, 1) == size(contact_mat, 2)
        nothing
    else
        error("contact_mat must be square")
    end
    if size(contact_mat, 1) == size(pop_matrix, 1)
        nothing
    else
        error("contact_mat and pop_matrix must have the same number of rows")
    end

    Q = contact_mat .* pop_matrix
    V = Diagonal(repeat([gamma + mu], size(contact_mat, 1)))

    K_prime = Q * inv(V)
    eigenvals = eigen(K_prime).values

    beta = R_0 / maximum(real(eigenvals))

    return beta
end

function calculate_beta(
        R_0,
        gamma,
        mu,
        contact_mat::T1,
        pop_matrix::T2
    ) where {T1 <: Union{<:AbstractFloat, <:Integer}, T2 <: Union{<:AbstractFloat, <:Integer}}
    return calculate_beta(
        convert(Float64, R_0),
        convert(Float64, gamma),
        convert(Float64, mu),
        fill(convert(Float64, contact_mat), 1, 1),
        fill(convert(Float64, pop_matrix), 1)
    )
end

"""
    calculate_gamma(
    	R_0::AbstractFloat,
    	beta::AbstractFloat,
    	mu::AbstractFloat,
		contact_mat::Matrix{AbstractFloat},
		pop_matrix::Vector{AbstractFloat}
    )

Calculate the value gamma for a given set of parameters, a contact matrix, and the population matrix, using the Next-Generation Matrix (K).
The contact matrix represents the contact rates between an individual in group 𝒾 and an individual in group 𝒿.
The population matrix is a vector that represents the number of individual in each group.

Uses the inverse relationship from calculate_beta. Given that:
R_0 = β ρ(Q⋅V⁻¹) where V = Diagonal([gamma + mu, ...])

Since V is diagonal with identical entries, V⁻¹ = (1/(gamma + mu)) * I
Therefore: R_0 = β * ρ(Q) / (gamma + mu)
Solving for gamma: γ = (β * ρ(Q) / R_0) - μ
"""
function calculate_gamma(
        R_0::T, beta::T, mu::T, contact_mat::Matrix{T}, pop_matrix::Vector{T}
    ) where {T <: AbstractFloat}
    # Validate input dimensions (same as calculate_beta)
    if size(contact_mat, 1) != size(contact_mat, 2)
        error("contact_mat must be square")
    end
    if size(contact_mat, 1) != size(pop_matrix, 1)
        error("contact_mat and pop_matrix must have the same number of rows")
    end

    Q = contact_mat .* pop_matrix

    # We need to solve for gamma such that:
    # R_0 = beta * max_eigenvalue(Q * inv(V))
    # where V = Diagonal([gamma + mu, gamma + mu, ...])

    # This becomes: R_0 = beta * max_eigenvalue(Q) / (gamma + mu)
    # Rearranging: gamma = (beta * max_eigenvalue(Q) / R_0) - mu

    eigenvals_Q = eigen(Q).values
    max_eigenval_Q = maximum(real(eigenvals_Q))

    gamma = (beta * max_eigenval_Q / R_0) - mu

    # Validate the result
    if gamma <= 0
        error("Calculated gamma is non-positive. Check parameter consistency.")
    end

    return gamma
end

function calculate_gamma(
        R_0,
        beta,
        mu,
        contact_mat::T1,
        pop_matrix::T2
    ) where {T1 <: Union{<:AbstractFloat, <:Integer}, T2 <: Union{<:AbstractFloat, <:Integer}}
    return calculate_gamma(
        convert(Float64, R_0),
        convert(Float64, beta),
        convert(Float64, mu),
        fill(convert(Float64, contact_mat), 1, 1),
        fill(convert(Float64, pop_matrix), 1)
    )
end


"""
    calculate_beta_amp(beta_mean, beta_force, t; seasonality)

Calculate the amplitude of the transmission rate beta as a function of time.
`beta_mean` is the mean transmission rate, `beta_force` is the amplitude of the `seasonality` function.
`seasonality` should be a SeasonalityFunction sum type.
"""
function calculate_beta_amp(beta_mean, beta_force, t; seasonality::SeasonalityFunction)
    return _calculate_beta_amp(beta_mean, beta_force, t, variant(seasonality))
end

# Dispatch on the extracted variant
_calculate_beta_amp(beta_mean, beta_force, t, ::CosineSeasonality) =
    beta_mean * (1 + beta_force * cos(2π * t / 365))

_calculate_beta_amp(beta_mean, beta_force, t, ::SineSeasonality) =
    beta_mean * (1 + beta_force * sin(2π * t / 365))

"""
    calculateReffective_t!(Reff_vec, beta_vec, dynamics_params, contact_mat, pop_matrix, seir_arr)

Calculate the effective reproduction number, R_eff, at each time step for a given set of parameters and contact matrix.
"""
function calculateReffective_t!(
        Reff_vec::AbstractVector{Float64},
        beta_vec::Vector{Float64},
        dynamics_params,
        contact_mat::Int64,
        seir_arr
    )::Nothing
    for i in eachindex(Reff_vec)
        Reff_vec[i] = calculateReffective(
            beta_vec[i],
            dynamics_params,
            contact_mat,
            seir_arr[i][1],
            seir_arr[i][5],
        )
    end

    return nothing
end

"""
    calculateReffective!(beta_t, dynamics_params, contact_mat, pop_matrix, S, N)

Calculate the effective reproduction number, R_eff, for a given set of parameters and contact matrix.
"""
function calculateReffective(
        beta_t::Float64,
        dynamics_params,
        contact_mat::Int64,
        S::Int64,
        N::Int64
    )::Float64
    s_prop = S / N

    Reff =
        calculateR0(
        beta_t, dynamics_params.gamma, dynamics_params.mu, contact_mat,
        N
    ) * s_prop

    return Reff
end

"""
    calculateR0(beta, gamma, mu, contact_mat, pop_matrix)

Calculate the basic reproduction number R_0 for a given set of parameters and contact matrix.

```jldoctest
julia> calculateR0(0.00025, 1 / 8, 0.0, ones(1, 1), [1_000])
2.0
```

* * *

**TODO** Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.

* * *
"""
function calculateR0(
        beta::T, gamma::T, mu::T, contact_mat::Matrix{T}, pop_matrix::Matrix{T}
    )::Float64 where {T <: AbstractFloat}
    if size(contact_mat, 1) == size(contact_mat, 2)
        nothing
    else
        error("contact_mat must be square")
    end
    if size(contact_mat, 1) == size(pop_matrix, 1)
        nothing
    else
        error("contact_mat and pop_matrix must have the same number of rows")
    end

    B = beta * contact_mat

    F = B .* pop_matrix
    V = Diagonal(repeat([gamma + mu], size(contact_mat, 1)))

    FV⁻¹ = F * inv(V)
    eigenvals = eigen(FV⁻¹).values

    R_0 = maximum(real(eigenvals))

    return R_0
end

function calculateR0(beta::Float64, gamma::Float64, mu::Float64, contact_mat::Int64, pop_matrix::Int64)::Float64
    return calculateR0(
        convert(Float64, beta),
        convert(Float64, gamma),
        convert(Float64, mu),
        reshape([Float64(contact_mat)], 1, 1),
        reshape([Float64(pop_matrix)], 1, 1),
    )
end

function calculate_mu(annual_births_per_k)
    life_expectancy_years = 1000 / annual_births_per_k
    return 1 / (life_expectancy_years * 365)
end

"""
    calculate_import_rate(mu, R_0, N)

Calulate the rate of new infectious individuals imported into the simulation using the commuter import formula defined in p210 of Keeling & Rohani
"""
function calculate_import_rate(mu, R_0, N)
    return (1.06 * mu * (R_0 - 1)) / sqrt(N)
end

# end
