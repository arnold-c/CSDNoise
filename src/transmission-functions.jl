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
    calculate_beta(R_0, gamma, mu, contact_mat, pop_matrix)

Calculate the value beta for a given set of parameters and contact matrix.

```jldoctest
julia> calculate_beta(2.0, 1 / 8, 0.0, ones(1, 1), [1_000])
0.00025
```

"""
# TODO: Currently only works when the populations are the same size as each other, and doesn't account for an exposed state.
function calculate_beta(
        R_0::T, gamma::T, mu::T, contact_mat::Array{T}, pop_matrix::Array{T}
    ) where {T <: AbstractFloat}
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

    F = contact_mat .* pop_matrix
    V = Diagonal(repeat([gamma + mu], size(contact_mat, 1)))

    FV⁻¹ = F * inv(V)
    eigenvals = eigen(FV⁻¹).values
    beta = R_0 / maximum(real(eigenvals))

    return beta
end

function calculate_beta(R_0, gamma, mu, contact_mat, pop_matrix)
    return calculate_beta(
        convert(Float64, R_0),
        convert(Float64, gamma),
        convert(Float64, mu),
        convert(Array{Float64}, [contact_mat]),
        convert(Array{Float64}, [pop_matrix]),
    )
end


"""
    calculate_beta_amp(beta_mean, beta_force, t; seasonality)

Calculate the amplitude of the transmission rate beta as a function of time.
`beta_mean` is the mean transmission rate, `beta_force` is the amplitude of the `seasonality` function.
`seasonality` should be a SeasonalityFunction sum type.
"""
function calculate_beta_amp(beta_mean, beta_force, t; seasonality)
    return calculate_beta_amp_impl(beta_mean, beta_force, t, variant(seasonality))
end

# Dispatch on the extracted variant
calculate_beta_amp_impl(beta_mean, beta_force, t, ::CosineSeasonality) =
    beta_mean * (1 + beta_force * cos(2π * t / 365))

calculate_beta_amp_impl(beta_mean, beta_force, t, ::SineSeasonality) =
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
        seir_arr::AbstractMatrix{Int64}
    )::Nothing
    for i in eachindex(Reff_vec)
        Reff_vec[i] = calculateReffective(
            beta_vec[i],
            dynamics_params,
            contact_mat,
            seir_arr[i, 1],
            seir_arr[i, 5],
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
