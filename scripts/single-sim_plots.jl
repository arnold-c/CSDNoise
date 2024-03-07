#%%
using DrWatson
@quickactivate "CSDNoise"

using JLD2
using DataFrames
using DataFramesMeta

using OutbreakDetection

include(srcdir("makie-plotting-setup.jl"))

@unpack singlesim_time_p = load(
    joinpath(outdir("singlesim"), "single-sim_setup.jld2")
)
@unpack trange = singlesim_time_p;

@unpack seir_array, inc_vec, beta_vec, seir_df = load(
    joinpath(outdir("singlesim"), "single-sim_arrays.jld2")
)

#%%
singlesim_timeseries_plot = draw_sir_plot(
    seir_df;
    annual = true,
    colors = seircolors,
    labels = ["S", "E", "I", "R", "N"],
)

mkpath(plotsdir("singlesim"))

save(plotsdir("singlesim/single-sim_timeseries.png"), singlesim_timeseries_plot)

#%%
singlesim_si_state_space_plot = @chain DataFrame(Tables.table(seir_array)) begin
    hcat(trange, _)
    rename!(["time", ["S", "E", "I", "R", "N"]...])
    data(_) *
    mapping(:I, :S; color = :time) *
    visual(Lines)
    draw
end

save(
    plotsdir("singlesim/single-sim_SI-state-space.png"),
    singlesim_si_state_space_plot,
)

#%%
singlesim_beta_plot = Figure()
beta_ax = Axis(
    singlesim_beta_plot[1, 1]; xlabel = "Time (years)", ylabel = "Beta"
)

lines!(beta_ax, trange / 365, beta_vec; linewidth = 1)

xlims!(beta_ax, (0, 3))

singlesim_beta_plot

save(plotsdir("singlesim/single-sim_beta.png"), singlesim_beta_plot)
