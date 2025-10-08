export StateParameters

struct StateParameters
    init_states::LabelledArrays.SLArray{Tuple{5}, Int64, 1, 5, (:S, :E, :I, :R, :N)}
    init_state_props::LabelledArrays.SLArray{Tuple{4}, Float64, 1, 4, (:s_prop, :e_prop, :i_prop, :r_prop)}
end

function StateParameters(N::Int64, init_state_props::Dict)
    return StateParameters(;
        N = N,
        s_prop = init_state_props[:s_prop],
        e_prop = init_state_props[:e_prop],
        i_prop = init_state_props[:i_prop],
    )
end

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
