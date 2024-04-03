# module CleaningFunctions
#
# export create_sir_df, create_sir_beta_dfs, create_sir_sim_array!,
#     create_sir_all_sims_array, create_sir_all_sims_array!,
#     create_sir_all_sim_quantiles, create_sir_all_sim_quantiles!

using DataFrames, DataFramesMeta

function create_sir_df(sir_array::Matrix, trange, states = [:S, :I, :R, :N])
    if size(sir_array, 1) == length(states)
        sir_array = sir_array'
    end
    return create_sir_df_inner(
        hcat(trange, DataFrame(Tables.table(sir_array))), states
    )
end

function create_sir_df_inner(sir_df::DataFrame, states)
    @chain sir_df begin
        rename!([:time, states...])
        if :N in states
            stack(_, [states...]; variable_name = :State, value_name = :Number)
        else
            transform!(_, states => (+) => :N)
            stack(
                _, [states..., :N]; variable_name = :State, value_name = :Number
            )
        end
    end
end

function create_sir_df(sol, states)
    return create_sir_df_inner(DataFrame(sol), states)
end

# end
