export create_combinations_vec

"""
    create_combinations_vec(custom_function, combinations; init = custom_function[])

Apply a custom function to all combinations of input parameters and concatenate results.

This utility function generates the Cartesian product of all provided parameter combinations
and applies a custom function to each combination, then concatenates all results into a
single vector. It's particularly useful for parameter sweeps and grid searches where you
need to evaluate a function across all possible combinations of input parameters.

# Arguments
- `custom_function`: Function to apply to each parameter combination. Should return a vector or collection that can be concatenated
- `combinations`: Tuple or collection of parameter vectors to combine. Each element represents the possible values for one parameter
- `init = custom_function[]`: Initial value for the reduction operation. Defaults to an empty result from `custom_function`

# Returns
- Vector containing concatenated results from applying `custom_function` to all parameter combinations

# Example
```julia
# Define a function that creates parameter sets
function create_param_set(R0, period)
    return [DiseaseParams(R0=R0, latent_period=period)]
end

# Define parameter ranges
R0_values = [1.5, 2.0, 2.5]
period_values = [3.0, 5.0, 7.0]

# Generate all combinations
all_params = create_combinations_vec(
    create_param_set,
    (R0_values, period_values)
)

# Results in 9 parameter sets (3 × 3 combinations)
length(all_params) # 9
```

# Notes
- The function uses `Iterators.product` to generate all combinations efficiently
- Results are concatenated using `mapreduce` with `vcat`
- Currently marked as TODO for type stability improvements
- The `custom_function` should return a collection that can be concatenated with `vcat`

# See Also
- [`Iterators.product`](@ref): For generating Cartesian products
- [`mapreduce`](@ref): For applying function and reducing results
"""
function create_combinations_vec(custom_function, combinations; init = custom_function[])
    combs = Iterators.product(combinations...)

    # TODO: make type stable
    return mapreduce(combination -> custom_function(combination...), vcat, combs; init = init)
end
