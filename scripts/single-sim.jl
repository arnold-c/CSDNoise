#%%
using DrWatson
@quickactivate "CSDNoise"

using UnPack
using JLD2

using CSDNoise

#%%
@unpack singlesim_states_p, singlesim_time_p, singlesim_dynamics_p = load(
    joinpath(outdir("singlesim"), "single-sim_setup.jld2")
)

#%%
seir_vec, inc_vec, beta_vec = seir_mod(
    singlesim_states_p.init_states,
    singlesim_dynamics_p,
    singlesim_time_p; seed = 1234,
)

seir_array = convert_svec_to_matrix(seir_vec)

seir_df = create_sir_df(
    seir_array, singlesim_time_p.trange, [:S, :E, :I, :R, :N]
)

#%%
jldsave(
    "out/singlesim/single-sim_arrays.jld2";
    seir_array,
    inc_vec,
    beta_vec,
    seir_df,
)
