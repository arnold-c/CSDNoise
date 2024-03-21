# module TransmissionFunctions

using LinearAlgebra

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
) where {T<:AbstractFloat}
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
    calculate_beta_amp(beta_mean, beta_force, t)

Calculate the amplitude of the transmission rate beta as a function of time.
`beta_mean` is the mean transmission rate, `beta_force` is the amplitude of the `seasonality` function.
`seasonality` defaults to using the `cosine` function
"""
function calculate_beta_amp(beta_mean, beta_force, t; seasonality = cos)
    return beta_mean * (1 + beta_force * seasonality(2pi * t / 365))
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
    beta::T, gamma::T, mu::T, contact_mat::Array{T}, pop_matrix::Array{T}
) where {T<:AbstractFloat}
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

function calculateR0(beta, gamma, mu, contact_mat, pop_matrix)
    return calculateR0(
        convert(Float64, beta),
        convert(Float64, gamma),
        convert(Float64, mu),
        convert(Array{Float64}, [contact_mat]),
        convert(Array{Float64}, [pop_matrix]),
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
