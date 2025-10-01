export create_combinations_vec

function create_combinations_vec(custom_function, combinations; init = custom_function[])
    combs = Iterators.product(combinations...)

    # TODO: make type stable
    return mapreduce(combination -> custom_function(combination...), vcat, combs; init = init)
end
