#!/usr/bin/env julia
"""
Test script to identify JET errors in ews_objective_function_with_tracking
"""

using DrWatson
@quickactivate("CSDNoise")

using CSDNoise
using JET
using Dates
using InteractiveUtils: @code_warntype

println("🔍 Testing JET errors in ews_objective_function_with_tracking")
println("="^60)

#%%
println("📊 Creating test data for profiling...")

# Create minimal but realistic test scenarios (similar to benchmark_ensemble_size.jl)
@report_opt target_modules = (CSDNoise,) create_single_test_scenario()
specification_vecs = create_single_test_scenario()

# Generate small ensemble data for profiling (use small number of simulations)
nsims = 50  # Small for profiling
println("  Generating ensemble specifications for $nsims simulations...")

ensemble_spec, null_spec, outbreak_spec = create_ensemble_specs(nsims)
@report_opt target_modules = (CSDNoise,) create_ensemble_specs(nsims)

println("  Generating ensemble data...")
@report_opt target_modules = (CSDNoise,) generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)
@report_call target_modules = (CSDNoise,) generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)
@code_warntype generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)
data_arrs = generate_ensemble_data(ensemble_spec, null_spec, outbreak_spec)


#%%
# Create test data
specification_vecs, data_arrs = create_test_data_for_profiling()
@report_opt target_modules = (CSDNoise,) create_test_data_for_profiling()

# Create scenario for testing
scenarios_vec = create_optimization_scenarios(specification_vecs)
scenario = scenarios_vec[1]

# Create cached data and tracker
cached_data = create_cached_simulation_data(scenario, data_arrs);
tracker = OptimizationTracker()
test_params = [0.9, 5.0]

#%%
println("Running JET analysis on ews_objective_function_with_tracking...")
println("-"^60)

#%%
# This is the call from line 63 that's causing errors
println("JET @report_opt analysis:")
result = @report_opt target_modules = (CSDNoise,) ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker);
println("JET result: ", result)

if !isempty(JET.get_reports(result))
    println(@code_warntype(ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker)))
end

#%%
println("\nRunning @report_call analysis:")
result2 = @report_call target_modules = (CSDNoise,) ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker);
println("JET call result: ", result2)

if !isempty(JET.get_reports(result2))
    println(@code_warntype(ews_objective_function_with_tracking(test_params, scenario, cached_data, tracker)))
end

#%%
println("\nDone!")
