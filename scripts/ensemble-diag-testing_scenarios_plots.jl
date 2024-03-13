#%%
using DrWatson
@quickactivate "CSDNoise"

using ProgressMeter
using FLoops
using NaNMath: NaNMath

using CSDNoise

includet(srcdir("makie-plotting-setup.jl"))
includet(srcdir("ensemble-parameters.jl"))

#%%
threshold_comparison_test_spec_vec = [
    IndividualTestSpecification(0.85, 0.85, 0),
    CLINICAL_CASE_TEST_SPEC,
    IndividualTestSpecification(1.0, 1.0, 0),
    IndividualTestSpecification(1.0, 1.0, 3),
    IndividualTestSpecification(1.0, 1.0, 7),
    IndividualTestSpecification(1.0, 1.0, 14),
]

threshold_comparison_alertthreshold_vec = [
    1, collect(2:2:14)...
]

alert_method_vec = ["movingavg", "dailythreshold_movingavg"]

#%%
for (ensemble_noise_specification, alertmethod) in
    Iterators.product(ensemble_noise_specification_vec, alert_method_vec)
    threshold_comparison_params = (
        test_spec_vec = threshold_comparison_test_spec_vec,
        alertthreshold_vec = threshold_comparison_alertthreshold_vec,
        ensemble_specification = ensemble_specification,
        noise_specification = ensemble_noise_specification,
        outbreak_specification = ensemble_outbreak_specification,
        moving_avg_detection_lag = ensemble_moving_avg_detection_lag,
        percent_visit_clinic = ensemble_percent_visit_clinic,
        alertmethod = alertmethod,
    )

    @showprogress for percent_clinic_tested in
                      ensemble_percent_clinic_tested_vec
        plot_all_threshold_comparisons(
            percent_clinic_tested, threshold_comparison_params
        )
    end

    GC.gc(true)

    @info "All plots saved for $(ensemble_noise_specification.noise_type) and $(alertmethod)"
    println("==============================================")
end
