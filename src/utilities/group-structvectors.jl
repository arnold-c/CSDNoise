export group_structvector

"""
    group_structvector(sv::StructVector{T}, props::Symbol...) where {T}

Group a StructVector by specified properties, returning a dictionary of grouped views.

Takes a StructVector and groups its rows by the values of the specified properties.
Returns a dictionary where keys are NamedTuples containing the property values,
and values are views of the original StructVector containing all rows that share
those property values.

# Arguments
- `sv::StructVector{T}`: The StructVector to group
- `props::Symbol...`: Property names to group by

# Returns
- `Dict{NamedTuple, SubArray}`: Dictionary mapping property value combinations to grouped data views

# Example
```julia
using StructArrays
sv = StructArray(a=[1,1,2,2], b=[:x,:y,:x,:y], c=[10,20,30,40])
groups = group_structvector(sv, :a, :b)
# Returns groups like (a=1, b=:x) => view of matching rows

"""
function group_structvector(sv::StructVector{T}, props::Symbol...) where {T}
    prop_arrays = tuple((getproperty(sv, p) for p in props)...)

    groups = Dict{NamedTuple{props}, Vector{Int}}()

    for i in eachindex(sv)
        key = NamedTuple{props}(tuple((prop_arrays[j][i] for j in 1:length(props))...))
        if haskey(groups, key)
            push!(groups[key], i)
        else
            groups[key] = [i]
        end
    end

    return Dict(k => view(sv, idxs) for (k, idxs) in groups)
end
