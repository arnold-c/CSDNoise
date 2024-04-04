# module EnsembleFunctions
#
# export run_ensemble_jump_prob, run_jump_prob, summarize_ensemble_jump_prob,
#     jump_prob_summary, get_ensemble_file

using DrWatson
using UnPack
using FLoops
using ProgressMeter

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

function create_combinations_vec(custom_function, combinations)
    combs = Iterators.product(combinations...)

    return vec(map(combination -> custom_function(combination...), combs))
end

function create_ensemble_spec_combinations(
    beta_force_vec,
    seasonality_vec,
    sigma_vec,
    gamma_vec,
    annual_births_per_k_vec,
    R_0_vec,
    vaccination_coverage_vec,
    N_vec,
    init_states_prop_dict,
    model_types_vec,
    time_p_vec,
    nsims_vec,
)
    ensemble_spec_combinations = Iterators.product(
        beta_force_vec,
        seasonality_vec,
        sigma_vec,
        gamma_vec,
        annual_births_per_k_vec,
        R_0_vec,
        vaccination_coverage_vec,
        N_vec,
        init_states_prop_dict,
        model_types_vec,
        time_p_vec,
        nsims_vec,
    )

    ensemble_spec_vec = Vector(undef, length(ensemble_spec_combinations))

    for (
        i,
        (
            beta_force,
            seasonality,
            sigma,
            gamma,
            annual_births_per_k,
            R_0,
            vaccination_coverage,
            N,
            init_states_prop,
            model_type,
            time_p,
            nsims,
        ),
    ) in enumerate(ensemble_spec_combinations)
        mu = calculate_mu(annual_births_per_k)
        beta_mean = calculate_beta(R_0, gamma, mu, 1, N)
        epsilon = calculate_import_rate(mu, R_0, N)

        ensemble_spec_vec[i] = EnsembleSpecification(
            model_type,
            StateParameters(
                N, init_states_prop
            ),
            DynamicsParameters(
                beta_mean,
                beta_force,
                seasonality,
                sigma,
                gamma,
                mu,
                annual_births_per_k,
                epsilon,
                R_0,
                vaccination_coverage,
            ),
            time_p,
            nsims,
        )
    end

    return ensemble_spec_vec
end

function run_ensemble_jump_prob(dict_of_ensemble_params; force = false)
    prog = Progress(length(dict_of_ensemble_params))
    for ensemble_params in dict_of_ensemble_params
        @produce_or_load(
            run_jump_prob,
            ensemble_params,
            "$(ensemble_params[:ensemble_spec].dirpath)";
            filename = "ensemble-solution",
            loadfile = false,
            force = force
        )
        next!(prog)
    end
end

"""
    run_jump_prob(ensemble_param_dict)
"""
function run_jump_prob(ensemble_param_dict)
    @unpack ensemble_spec,
    seed,
    executor,
    quantile_vec,
    outbreak_spec_dict,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec,
    ews_spec_vec = ensemble_param_dict

    @unpack state_parameters, dynamics_parameters, time_parameters, nsims =
        ensemble_spec

    @unpack tstep, tlength, trange = time_parameters

    ensemble_seir_vecs = Array{typeof(state_parameters.init_states),2}(
        undef,
        tlength,
        nsims
    )

    ensemble_inc_vecs = Array{typeof(SVector(0)),2}(
        undef,
        tlength,
        nsims,
    )

    ensemble_beta_arr = zeros(Float64, tlength)

    ensemble_Reff_arr = zeros(Float64, tlength, nsims)
    ensemble_Reff_thresholds_vec = Vector{Array{Int64,2}}(
        undef, size(ensemble_inc_vecs, 2)
    )

    for sim in axes(ensemble_inc_vecs, 2)
        run_seed = seed + (sim - 1)

        seir_mod!(
            @view(ensemble_seir_vecs[:, sim]),
            @view(ensemble_inc_vecs[:, sim]),
            ensemble_beta_arr,
            state_parameters.init_states,
            dynamics_parameters,
            time_parameters;
            seed = run_seed,
        )
    end

    ensemble_seir_arr = convert_svec_to_array(ensemble_seir_vecs)

    for sim in axes(ensemble_inc_vecs, 2)
        calculateReffective_t!(
            @view(ensemble_Reff_arr[:, sim]),
            ensemble_beta_arr,
            dynamics_parameters,
            1,
            @view(ensemble_seir_arr[:, :, sim]),
        )

        ensemble_Reff_thresholds_vec[sim] = Reff_ge_than_one(
            @view(ensemble_Reff_arr[:, sim])
        )
    end

    for dict in outbreak_spec_dict
        dict[:dirpath] = joinpath(
            ensemble_spec.dirpath, dict[:outbreak_spec].dirpath
        )
        dict[:ensemble_spec] = ensemble_spec
        dict[:ensemble_inc_vecs] = ensemble_inc_vecs
        dict[:noise_spec_vec] = noise_spec_vec
        dict[:outbreak_detection_spec_vec] = outbreak_detection_spec_vec
        dict[:test_spec_vec] = test_spec_vec
        dict[:ews_spec_vec] = ews_spec_vec
        dict[:seed] = seed
        dict[:executor] = executor
    end

    run_define_outbreaks(outbreak_spec_dict)

    return @strdict ensemble_seir_arr ensemble_spec ensemble_Reff_arr ensemble_Reff_thresholds_vec
end

function run_define_outbreaks(dict_of_outbreak_spec_params)
    @floop for outbreak_spec_params in dict_of_outbreak_spec_params
        @produce_or_load(
            define_outbreaks,
            outbreak_spec_params,
            "$(outbreak_spec_params[:dirpath])";
            filename = "ensemble-incidence-array",
            loadfile = false,
            force = true
        )
    end
end

function define_outbreaks(incidence_param_dict)
    @unpack ensemble_spec,
    ensemble_inc_vecs,
    outbreak_spec,
    noise_spec_vec,
    outbreak_detection_spec_vec,
    test_spec_vec, ews_spec_vec, seed, executor =
        incidence_param_dict

    ensemble_inc_arr, ensemble_thresholds_vec = create_inc_infec_arr(
        ensemble_inc_vecs, outbreak_spec
    )

    non_clinical_case_test_spec_vec = filter(
        spec -> !(spec in CLINICAL_TEST_SPECS),
        test_spec_vec
    )

    non_clinical_case_outbreak_detection_spec_vec = filter(
        spec -> spec.percent_clinic_tested !== 1.0,
        outbreak_detection_spec_vec
    )

    non_clinical_case_ensemble_scenarios = create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_spec],
            [outbreak_spec],
            noise_spec_vec,
            non_clinical_case_outbreak_detection_spec_vec,
            non_clinical_case_test_spec_vec,
            ews_spec_vec,
        ),
    )

    clinical_case_outbreak_detection_spec_vec = filter(
        spec -> spec.percent_clinic_tested == 1.0,
        outbreak_detection_spec_vec
    )

    clinical_case_ensemble_scenarios = create_combinations_vec(
        ScenarioSpecification,
        (
            [ensemble_spec],
            [outbreak_spec],
            noise_spec_vec,
            clinical_case_outbreak_detection_spec_vec,
            CLINICAL_TEST_SPECS,
            ews_spec_vec,
        ),
    )

    ensemble_scenarios = vcat(
        non_clinical_case_ensemble_scenarios, clinical_case_ensemble_scenarios
    )

    scenario_param_dict = dict_list(
        @dict(
            scenario_spec = ensemble_scenarios,
            ensemble_inc_arr,
            thresholds_vec = [ensemble_thresholds_vec],
            seed = seed,
        )
    )

    run_OutbreakThresholdChars_creation(
        scenario_param_dict; executor = executor
    )

    return @strdict ensemble_inc_arr ensemble_thresholds_vec
end

function run_OutbreakThresholdChars_creation(
    dict_of_OTchars_params; executor = ThreadedEx()
)
    @floop executor for OTChars_params in dict_of_OTchars_params
        @produce_or_load(
            OutbreakThresholdChars_creation,
            OTChars_params,
            "$(OTChars_params[:scenario_spec].dirpath)";
            filename = "ensemble-scenario",
            loadfile = false
        )
    end
end

function OutbreakThresholdChars_creation(OT_chars_param_dict)
    @unpack scenario_spec,
    ensemble_inc_arr,
    thresholds_vec,
    seed = OT_chars_param_dict

    @unpack noise_specification,
    ensemble_specification,
    outbreak_specification,
    outbreak_detection_specification,
    individual_test_specification,
    ews_specification = scenario_spec

    noise_array, noise_rubella_prop = create_noise_arr(
        noise_specification,
        ensemble_inc_arr;
        ensemble_specification = scenario_spec.ensemble_specification,
        seed = seed,
    )

    testarr, ewsarr = create_testing_arrs(
        ensemble_inc_arr,
        noise_array,
        outbreak_detection_specification,
        individual_test_specification,
        ensemble_specification.time_parameters,
        ews_specification,
    )[1:2]

    OT_chars = calculate_OutbreakThresholdChars(
        testarr, ensemble_inc_arr, thresholds_vec, noise_rubella_prop
    )

    return @strdict OT_chars ewsarr
end

function get_ensemble_file() end

function get_ensemble_file(
    ensemble_spec::EnsembleSpecification, outbreak_spec::OutbreakSpecification
)
    dirpath = joinpath(ensemble_spec.dirpath, outbreak_spec.dirpath)

    return load(collect_ensemble_file("incidence-array", dirpath)...)
end

function get_ensemble_file(spec::EnsembleSpecification)
    return load(collect_ensemble_file("solution", spec.dirpath)...)
end

function get_ensemble_file(spec::EnsembleSpecification, quantile::Int64)
    return load(collect_ensemble_file("quantiles_$(quantile)", spec.dirpath)...)
end

function get_ensemble_file(spec::ScenarioSpecification)
    return load(collect_ensemble_file("scenario", spec.dirpath)...)
end

function collect_ensemble_file(type, dirpath)
    filecontainer = []
    for f in readdir(dirpath)
        match_ensemble_file!(type, dirpath, filecontainer, f)
    end
    if length(filecontainer) != 1
        println("Matched $(length(filecontainer)) files, when should be 1")
    end
    return filecontainer
end

function match_ensemble_file!(criteria, dirpath, container, file)
    if occursin(criteria, file)
        push!(container, joinpath(dirpath, file))
    end
end

# end
