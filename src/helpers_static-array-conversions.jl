export convert_svec_to_matrix,
    convert_svec_to_matrix!,
    convert_svec_to_array

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
