using StructArrays: StructVector
using UnPack: @unpack
using DataFrames: DataFrame

export create_scenarios_structvector,
    create_scenarios_dataframe,
    create_optimization_scenarios

# Keep old function for backward compatibility
"""
    create_scenarios_dataframe(specification_vecs)

Create a DataFrame containing all scenario combinations for optimization.
This replaces the nested loop approach with a DataFrame-based approach for better caching.
"""
function create_scenarios_dataframe(specification_vecs)
    scenarios = create_scenarios_structvector(specification_vecs)
    return DataFrame(scenarios)
end

function create_scenarios_structvector(specification_vecs)
    scenarios = create_optimization_scenarios(specification_vecs)
    return StructVector(scenarios)
end

"""
    create_optimization_scenarios(specification_vecs)

Create all scenario combinations for optimization.
Validates specification vectors before creating combinations.
"""
@unstable function create_optimization_scenarios(specification_vecs)
    # Validate specification vectors before processing
    _validate_specification_vectors(specification_vecs)

    @unpack ensemble_specification_vec,
        null_specification_vec,
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec = specification_vecs

    # Use Iterators.product instead of nested loops
    combinations = Iterators.product(
        ensemble_specification_vec,
        null_specification_vec,
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec
    ) |> collect

    scenarios_vec = mapreduce(vcat, combinations; init = OptimizationScenario[]) do (
                ensemble_spec,
                null_spec,
                noise_spec,
                test_spec,
                percent_tested,
                ews_metric_spec,
                ews_enddate_type,
                ews_window,
                ews_burnin,
                ews_metric,
            )
        OptimizationScenario(
            ensemble_spec,
            null_spec,
            noise_spec,
            test_spec,
            percent_tested,
            ews_metric_spec,
            ews_enddate_type,
            ews_window,
            ews_burnin,
            ews_metric,
        )
    end

    return scenarios_vec
end

"""
    validate_specification_vectors(specification_vecs)

Validate that specification vectors contain valid values and types.
Throws ArgumentError for invalid specifications and warns for empty vectors.
"""
function _validate_specification_vectors(specification_vecs)
    # Define required fields
    required_fields = Set(
        [
            :ensemble_specification_vec,
            :null_specification_vec,
            :noise_specification_vec,
            :test_specification_vec,
            :percent_tested_vec,
            :ews_metric_specification_vec,
            :ews_enddate_type_vec,
            :ews_threshold_window_vec,
            :ews_threshold_burnin_vec,
            :ews_metric_vec,
        ]
    )

    # Get provided fields
    provided_fields = Set(keys(specification_vecs))

    # Check for missing required fields
    missing_fields = setdiff(required_fields, provided_fields)
    if !isempty(missing_fields)
        throw(ArgumentError("Missing required fields: $(join(sort(collect(missing_fields)), ", "))"))
    end

    # Check for extra fields
    extra_fields = setdiff(provided_fields, required_fields)
    if !isempty(extra_fields)
        throw(ArgumentError("Unexpected extra fields provided: $(join(sort(collect(extra_fields)), ", ")). Only these fields are allowed: $(join(sort(collect(required_fields)), ", "))"))
    end

    @unpack ensemble_specification_vec,
        null_specification_vec,
        noise_specification_vec,
        test_specification_vec,
        percent_tested_vec,
        ews_metric_specification_vec,
        ews_enddate_type_vec,
        ews_threshold_window_vec,
        ews_threshold_burnin_vec,
        ews_metric_vec = specification_vecs

    # Validate ensemble specifications
    if !all(spec -> spec isa EnsembleSpecification, vcat(ensemble_specification_vec, null_specification_vec))
        throw(ArgumentError("All emergent and null specifications must be of type EnsembleSpecification"))
    end

    # Validate noise specifications
    if !all(spec -> spec isa NoiseSpecification, noise_specification_vec)
        throw(ArgumentError("All noise specifications must be of type NoiseSpecification"))
    end

    # Validate test specifications
    if !all(spec -> spec isa IndividualTestSpecification, test_specification_vec)
        throw(ArgumentError("All test specifications must be of type IndividualTestSpecification"))
    end

    # Validate percent tested values
    if !all(p -> 0.0 <= p <= 1.0, percent_tested_vec)
        throw(ArgumentError("All percent_tested values must be between 0.0 and 1.0"))
    end

    # Validate EWS metric specifications
    # if !all(spec -> spec isa EWSMetricSpecification, ews_metric_specification_vec)
    #     throw(ArgumentError("All EWS metric specifications must be of type EWSMetricSpecification"))
    # end

    # Validate EWS end date types
    if !all(enddate -> enddate isa EWSEndDateType, ews_enddate_type_vec)
        throw(ArgumentError("All EWS end date types must be of type EWSEndDateType"))
    end

    # Validate threshold window types
    if !all(window -> window isa EWSThresholdWindowType, ews_threshold_window_vec)
        throw(ArgumentError("All threshold windows must be either ExpandingThresholdWindow or RollingThresholdWindow"))
    end

    # Validate burnin periods
    if !all(burnin -> burnin isa Union{Dates.Day, Dates.Year}, ews_threshold_burnin_vec)
        throw(ArgumentError("All burnin periods must be of type Day or Year"))
    end

    # Validate EWS metric names
    valid_metrics = [
        "autocorrelation", "autocovariance", "coefficient_of_variation",
        "index_of_dispersion", "kurtosis", "mean", "skewness", "variance",
    ]
    if !all(metric -> metric in valid_metrics, ews_metric_vec)
        invalid_metrics = setdiff(ews_metric_vec, valid_metrics)
        throw(ArgumentError("Invalid EWS metrics: $(invalid_metrics). Valid metrics are: $(valid_metrics)"))
    end

    # Check for empty vectors (warn but don't error)
    empty_fields = String[]
    for field in required_fields
        if isempty(getfield(specification_vecs, field))
            push!(empty_fields, string(field))
        end
    end

    if !isempty(empty_fields)
        @warn "Empty specification vectors detected: $(join(empty_fields, ", ")). This will result in zero scenarios."
    end

    return nothing
end
