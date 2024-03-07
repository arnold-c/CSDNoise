# module CleaningFunctions
#
# export create_sir_df, create_sir_beta_dfs, create_sir_sim_array!,
#     create_sir_all_sims_array, create_sir_all_sims_array!,
#     create_sir_all_sim_quantiles, create_sir_all_sim_quantiles!

using DataFrames, DataFramesMeta
using ModelingToolkit, DifferentialEquations
using FLoops

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

function create_sir_df(sol::RODESolution, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        if "value1" in names(_)
            rename!(_, [:time, states...])
        else
            rename!(s -> replace(s, "(t)" => ""), _)
            rename!(:timestamp => :time)
        end
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

function create_sir_df(sol::ODESolution, states = [:S, :I, :R])
    @chain DataFrame(sol) begin
        if "value1" in names(DataFrame(_))
            rename!(_, [:time, states...])
        else
            rename!(s -> replace(s, "(t)" => ""), _)
            rename!(:timestamp => :time)
        end
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

function create_sir_beta_dfs(sol, states = [:S, :I, :R])
    state_df = create_sir_df(sol, states)

    beta_df = select(state_df, [:time, :beta])
    unique!(beta_df)

    select!(state_df, Not(:beta))

    return state_df, beta_df
end

function create_sir_sim_array!(; jump_sol)
    sir_array[1:3, :] = Array(jump_sol)
    sir_array[4, :] = sum(sir_array[1:3, :]; dims = 1)

    return nothing
end

function create_sir_all_sim_quantiles!(all_sims_array, sim_quantiles; quantiles)
    for state in axes(all_sims_array, 2), time in axes(all_sims_array, 1)
            sim_quantiles[:, time, state] = quantile(
                skipmissing(
                    replace(@view(all_sims_array[time, state, :]), NaN => missing)
                ),
                quantiles,
            )
    end
end

function create_sir_all_sim_quantiles(all_sims_array; quantiles)
    quantile_array = zeros(
        length(quantiles), size(all_sims_array, 1), size(all_sims_array, 2)
    )

    create_sir_all_sim_quantiles!(
        all_sims_array, quantile_array; quantiles = quantiles
    )

    return quantile_array
end

# end
