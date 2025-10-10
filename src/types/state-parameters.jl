export StateParameters

"""
    StateParameters

Initial population state configuration for SEIR disease modeling.

This struct defines the initial population distribution across the SEIR compartments
(Susceptible, Exposed, Infected, Recovered) and stores both absolute counts and
proportions. The structure uses LabelledArrays for efficient access to compartment
values by name rather than index.

# Fields
- `init_states::LabelledArrays.SLArray{Tuple{5}, Int64, 1, 5, (:S, :E, :I, :R, :N)}`: Initial population counts for each compartment (S, E, I, R) plus total population N
- `init_state_props::LabelledArrays.SLArray{Tuple{4}, Float64, 1, 4, (:s_prop, :e_prop, :i_prop, :r_prop)}`: Proportions of population in each compartment (must sum to 1.0)

# Constructors
    StateParameters(N::Int64, init_state_props::Dict)
    StateParameters(; N=500_000, s_prop=0.1, e_prop=0.01, i_prop=0.01)

# Example
```julia
# Using keyword constructor with default population
state_params = StateParameters(
    s_prop = 0.8,   # 80% susceptible
    e_prop = 0.05,  # 5% exposed
    i_prop = 0.05   # 5% infected
    # r_prop calculated automatically as 1 - (s_prop + e_prop + i_prop) = 0.1
)

# Using dictionary constructor
props_dict = Dict(:s_prop => 0.9, :e_prop => 0.05, :i_prop => 0.05)
state_params = StateParameters(1_000_000, props_dict)

# Accessing values
state_params.init_states.S  # Number of susceptible individuals
state_params.init_state_props.s_prop  # Proportion susceptible
```

# See Also
- [`DynamicsParameters`](@ref): Disease dynamics parameters that work with initial states
- [`simulate_ensemble_seir_results`](@ref): Uses StateParameters for simulation initialization
"""
struct StateParameters
    init_states::LabelledArrays.SLArray{Tuple{5}, Int64, 1, 5, (:S, :E, :I, :R, :N)}
    init_state_props::LabelledArrays.SLArray{Tuple{4}, Float64, 1, 4, (:s_prop, :e_prop, :i_prop, :r_prop)}
end

"""
    StateParameters(N::Int64, init_state_props::Dict)

Construct StateParameters from total population and a dictionary of initial state proportions.

This constructor provides an alternative interface using a dictionary to specify the
initial proportions of the population in each compartment. The recovered proportion
is calculated automatically to ensure all proportions sum to 1.0.

# Arguments
- `N::Int64`: Total population size
- `init_state_props::Dict`: Dictionary containing initial state proportions with keys `:s_prop`, `:e_prop`, `:i_prop`

# Returns
- `StateParameters`: Configured state parameters with calculated absolute counts and proportions

# Example
```julia
props = Dict(
    :s_prop => 0.85,
    :e_prop => 0.05,
    :i_prop => 0.10
)
state_params = StateParameters(500_000, props)
```
"""
function StateParameters(N::Int64, init_state_props::Dict)
    return StateParameters(;
        N = N,
        s_prop = init_state_props[:s_prop],
        e_prop = init_state_props[:e_prop],
        i_prop = init_state_props[:i_prop],
    )
end

"""
    StateParameters(; N=500_000, s_prop=0.1, e_prop=0.01, i_prop=0.01)

Construct StateParameters using keyword arguments for population size and initial state proportions.

This is the primary constructor that creates a StateParameters instance with sensible defaults
for a typical disease modeling scenario. The recovered proportion is automatically calculated
to ensure all proportions sum to 1.0, and absolute population counts are derived by rounding
the proportional values.

# Keyword Arguments
- `N::Int64=500_000`: Total population size
- `s_prop::Float64=0.1`: Proportion of population initially susceptible (0.0 to 1.0)
- `e_prop::Float64=0.01`: Proportion of population initially exposed (0.0 to 1.0)
- `i_prop::Float64=0.01`: Proportion of population initially infected (0.0 to 1.0)

# Returns
- `StateParameters`: Configured state parameters with calculated recovered proportion and absolute counts

# Notes
- The recovered proportion `r_prop` is calculated as `1 - (s_prop + e_prop + i_prop)`
- All proportions must be non-negative and their sum must not exceed 1.0
- Absolute counts are calculated by rounding `proportion * N` to the nearest integer

# Example
```julia
# Default parameters (10% susceptible, 1% exposed, 1% infected, 88% recovered)
default_state = StateParameters()

# Custom endemic equilibrium state
endemic_state = StateParameters(
    N = 1_000_000,
    s_prop = 0.2,    # 20% susceptible
    e_prop = 0.001,  # 0.1% exposed
    i_prop = 0.001   # 0.1% infected
    # r_prop = 0.798 (calculated automatically)
)
```
"""
function StateParameters(;
        N = 500_000,
        s_prop = 0.1,
        e_prop = 0.01,
        i_prop = 0.01
    )
    r_prop = 1 - (s_prop + e_prop + i_prop)

    states = LabelledArrays.SLVector(;
        S = Int64(round(s_prop * N)),
        E = Int64(round(e_prop * N)),
        I = Int64(round(i_prop * N)),
        R = Int64(round(r_prop * N)),
        N = N,
    )
    state_props = LabelledArrays.SLVector(;
        s_prop = s_prop,
        e_prop = e_prop,
        i_prop = i_prop,
        r_prop = r_prop,
    )

    return StateParameters(
        states, state_props
    )
end
