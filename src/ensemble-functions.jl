# module EnsembleFunctions
#
# export run_ensemble_jump_prob, run_jump_prob, summarize_ensemble_jump_prob,
#     jump_prob_summary, get_ensemble_file

using DrWatson
using UnPack
using FLoops
using ProgressMeter
using StructArrays

# include("transmission-functions.jl")
# # using .TransmissionFunctions
#
# include("cleaning-functions.jl")
# # using .CleaningFunctions
#
# include("SEIR-model.jl")
# # using .SEIRModel
#
# include("structs.jl")
# using .ODStructs

function create_combinations_vec(custom_function, combinations; init = custom_function[])
    combs = Iterators.product(combinations...)

    # TODO: make type stable
    return mapreduce(combination -> custom_function(combination...), vcat, combs; init = init)
end

# function create_ensemble_spec_combinations(
#         beta_force_vec,
#         seasonality_vec,
#         sigma_vec,
#         gamma_vec,
#         annual_births_per_k_vec,
#         R_0_vec,
#         burnin_vaccination_coverage_params_vec,
#         vaccination_coverage_params_vec,
#         N_vec,
#         init_states_prop_dict,
#         model_types_vec,
#         time_p_vec,
#         nsims_vec,
#     )
#     ensemble_spec_combinations = Iterators.product(
#         beta_force_vec,
#         seasonality_vec,
#         sigma_vec,
#         gamma_vec,
#         annual_births_per_k_vec,
#         R_0_vec,
#         burnin_vaccination_coverage_params_vec,
#         vaccination_coverage_params_vec,
#         N_vec,
#         init_states_prop_dict,
#         model_types_vec,
#         time_p_vec,
#         nsims_vec,
#     )
#
#     ensemble_spec_vec = Vector(undef, length(ensemble_spec_combinations))
#
#     for (
#             i,
#             (
#                 beta_force,
#                 seasonality,
#                 sigma,
#                 gamma,
#                 annual_births_per_k,
#                 R_0,
#                 burnin_vaccination_coverage_params,
#                 vaccination_coverage_pairs,
#                 N,
#                 init_states_prop,
#                 model_type,
#                 time_p,
#                 nsims,
#             ),
#         ) in enumerate(ensemble_spec_combinations)
#         mu = calculate_mu(annual_births_per_k)
#         beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
#         epsilon = calculate_import_rate(mu, R_0, N)
#
#         states_params = StateParameters(
#             N, init_states_prop
#         )
#
#         vaccination_coverage_params = if burnin_vaccination_coverage_params isa Tuple{Nothing, Float64, Float64}
#             (
#                 (
#                     burnin_vaccination_coverage_params...,
#                     Int64((time_p.burnin / (365.0 / time_p.tstep)) * 2),
#                 ),
#                 vaccination_coverage_pairs,
#                 states_params.init_states,
#             )
#         elseif burnin_vaccination_coverage_params isa Tuple{Nothing, Float64, Float64, Int64}
#             (
#                 burnin_vaccination_coverage_params,
#                 vaccination_coverage_pairs,
#                 states_params.init_states,
#             )
#         elseif burnin_vaccination_coverage_params isa Tuple{Nothing, Float64}
#             (
#                 burnin_vaccination_coverage_params,
#                 vaccination_coverage_pairs,
#             )
#         else
#             (
#                 burnin_vaccination_coverage_params...,
#                 vaccination_coverage_pairs...,
#             )
#         end
#
#         dynamics_params_specification = DynamicsParameterSpecification(
#             beta_mean,
#             beta_force,
#             seasonality,
#             sigma,
#             gamma,
#             mu,
#             annual_births_per_k,
#             epsilon,
#             R_0,
#             vaccination_coverage_params...,
#         )
#
#         ensemble_spec_vec[i] = EnsembleSpecification(
#             model_type,
#             states_params,
#             dynamics_params_specification,
#             time_p,
#             nsims,
#         )
#     end
#
#     return ensemble_spec_vec
# end
#
# function run_ensemble_jump_prob(dict_of_ensemble_params; force = false)
#     prog = Progress(length(dict_of_ensemble_params))
#     for ensemble_params in dict_of_ensemble_params
#         @produce_or_load(
#             run_jump_prob,
#             ensemble_params,
#             "$(ensemble_params[:ensemble_spec].dirpath)";
#             filename = "ensemble-solution",
#             loadfile = false,
#             force = force
#         )
#         next!(prog)
#     end
#     return
# end
#
# """
#     run_jump_prob(ensemble_param_dict)
# """
# function run_jump_prob(ensemble_param_dict)
#     @unpack ensemble_spec,
#         seed,
#         executor,
#         outbreak_spec_dict = ensemble_param_dict
#
#     @unpack state_parameters,
#         dynamics_parameter_specification, time_parameters,
#         nsims =
#         ensemble_spec
#
#     @unpack tstep, tlength, trange = time_parameters
#
#     ensemble_seir_vecs = Array{typeof(state_parameters.init_states), 2}(
#         undef,
#         tlength,
#         nsims,
#     )
#
#     ensemble_inc_vecs = Array{typeof(SVector(0)), 2}(
#         undef,
#         tlength,
#         nsims,
#     )
#
#     ensemble_beta_arr = zeros(Float64, tlength)
#
#     ensemble_Reff_arr = zeros(Float64, tlength, nsims)
#     ensemble_Reff_thresholds_vec = Vector{Array{Int64, 2}}(
#         undef, size(ensemble_inc_vecs, 2)
#     )
#
#     dynamics_parameters = Vector{DynamicsParameters}(undef, nsims)
#
#     for sim in axes(ensemble_inc_vecs, 2)
#         run_seed = seed + (sim - 1)
#
#         dynamics_parameters[sim] = DynamicsParameters(
#             dynamics_parameter_specification; seed = run_seed
#         )
#
#         seir_mod!(
#             @view(ensemble_seir_vecs[:, sim]),
#             @view(ensemble_inc_vecs[:, sim]),
#             ensemble_beta_arr,
#             state_parameters.init_states,
#             dynamics_parameters[sim],
#             time_parameters,
#             run_seed,
#         )
#     end
#
#     ensemble_seir_arr = convert_svec_to_array(ensemble_seir_vecs)
#
#     for sim in axes(ensemble_inc_vecs, 2)
#         calculateReffective_t!(
#             @view(ensemble_Reff_arr[:, sim]),
#             ensemble_beta_arr,
#             dynamics_parameters[sim],
#             1,
#             @view(ensemble_seir_arr[:, :, sim]),
#         )
#
#         ensemble_Reff_thresholds_vec[sim] = Reff_ge_than_one(
#             @view(ensemble_Reff_arr[:, sim])
#         )
#     end
#
#     for dict in outbreak_spec_dict
#         dict[:dirpath] = joinpath(
#             ensemble_spec.dirpath, dict[:outbreak_spec].dirpath
#         )
#         dict[:ensemble_spec] = ensemble_spec
#         dict[:ensemble_inc_vecs] = ensemble_inc_vecs
#     end
#
#     run_define_outbreaks(outbreak_spec_dict; executor = executor)
#
#     return @strdict ensemble_seir_arr ensemble_spec dynamics_parameters ensemble_Reff_arr ensemble_Reff_thresholds_vec
# end
#
# function run_define_outbreaks(
#         dict_of_outbreak_spec_params; executor = ThreadedEx()
#     )
#     return @floop executor for outbreak_spec_params in dict_of_outbreak_spec_params
#         @produce_or_load(
#             define_outbreaks,
#             outbreak_spec_params,
#             "$(outbreak_spec_params[:dirpath])";
#             filename = "ensemble-incidence-array",
#             loadfile = false,
#             force = true
#         )
#     end
# end
#
# function define_outbreaks(incidence_param_dict)
#     @unpack ensemble_spec,
#         ensemble_inc_vecs,
#         outbreak_spec = incidence_param_dict
#
#     ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
#         ensemble_inc_vecs, outbreak_spec
#     )
#
#     return @strdict ensemble_inc_arr ensemble_thresholds_vec
# end
#
# function get_ensemble_file() end
#
# function get_ensemble_file(
#         ensemble_spec::EnsembleSpecification, outbreak_spec::OutbreakSpecification
#     )
#     dirpath = joinpath(ensemble_spec.dirpath, outbreak_spec.dirpath)
#
#     return load(collect_ensemble_file("incidence-array", dirpath)...)
# end
#
# function get_ensemble_file(spec::EnsembleSpecification)
#     return load(collect_ensemble_file("solution", spec.dirpath)...)
# end
#
# function get_ensemble_file(
#         ensemble_spec::EnsembleSpecification,
#         outbreak_spec::OutbreakSpecification,
#         ews_spec::EWSMetricSpecification,
#     )
#     return load(
#         collect_ensemble_file(
#             "incidence-ews-metrics",
#             joinpath(
#                 ensemble_spec.dirpath, outbreak_spec.dirpath, ews_spec.dirpath
#             ),
#         )...,
#     )
# end
#
# function get_ensemble_file(spec::ScenarioSpecification)
#     return load(collect_ensemble_file("scenario", spec.dirpath)...)
# end
#
# function collect_ensemble_file(type, dirpath)
#     filecontainer = []
#     for f in readdir(dirpath)
#         match_ensemble_file!(type, dirpath, filecontainer, f)
#     end
#     if length(filecontainer) != 1
#         println("Matched $(length(filecontainer)) files, when should be 1")
#     end
#     return filecontainer
# end
#
# function match_ensemble_file!(criteria, dirpath, container, file)
#     return if occursin(criteria, file)
#         push!(container, joinpath(dirpath, file))
#     end
# end
#
# # end
