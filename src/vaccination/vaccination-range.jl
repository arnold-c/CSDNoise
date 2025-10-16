export calculate_min_max_vaccination_range,
    calculate_post_burnin_vaccination_range

"""
    calculate_min_max_vaccination_range(mean_vaccination_coverage, max_vaccination_range = 0.2)

Calculate the minimum and maximum vaccination coverage range around a mean value.

This function determines a symmetric range around the mean vaccination coverage,
constrained by the maximum allowed range and the bounds [0, 1]. The actual range
used is the minimum of:
- The specified max_vaccination_range
- The distance to the upper bound (1.0 - mean_vaccination_coverage)
- The distance to the lower bound (mean_vaccination_coverage)

# Arguments
- `mean_vaccination_coverage`: Mean vaccination coverage (must be ≤ 1.0)
- `max_vaccination_range`: Maximum allowed range around the mean (default: 0.2)

# Returns
- `Tuple{Float64, Float64}`: (min_vaccination_coverage, max_vaccination_coverage)
"""
function calculate_min_max_vaccination_range(
        mean_vaccination_coverage,
        max_vaccination_range = 0.2,
    )
    @assert mean_vaccination_coverage <= 1.0

    min_vaccination_range = minimum(
        [
            max_vaccination_range,
            1.0 - mean_vaccination_coverage,
            mean_vaccination_coverage,
        ]
    )

    min_vaccination_coverage = round(
        mean_vaccination_coverage - min_vaccination_range; digits = 4
    )
    max_vaccination_coverage = round(
        mean_vaccination_coverage + min_vaccination_range; digits = 4
    )
    return min_vaccination_coverage, max_vaccination_coverage
end

"""
    calculate_post_burnin_vaccination_range(
        time_parameters, state_parameters, R_0;
        births_per_k_pop = 27.0, target_Reff = 1.0,
        estimated_post_burnin_Reff = 0.8, adjust_vaccination_coverage = -0.1,
        vaccination_bounds_spread = 0.2, reach_target_prop_of_remaining_simulation = 0.75,
        digits = 1
    )

Calculate the minimum and maximum vaccination coverage range for the post-burnin period.

This function determines vaccination coverage bounds needed to achieve a target effective
reproduction number (Reff) during the post-burnin simulation period. It calculates the
maximum vaccination rate required to reach the target Reff within a specified proportion
of the remaining simulation time, then applies adjustments and creates a range around
this value.

The calculation involves:
1. Computing time available to reach target Reff (proportion of post-burnin period)
2. Estimating required vaccination rate using current susceptible population
3. Applying floor adjustment and coverage adjustment
4. Creating min/max bounds with specified spread

# Arguments
- `time_parameters`: Simulation time parameters containing burnin and total length
- `state_parameters`: State parameters containing initial population states
- `R_0`: Basic reproduction number

# Keyword Arguments
- `births_per_k_pop`: Annual births per thousand population (default: 27.0)
- `target_Reff`: Target effective reproduction number (default: 1.0)
- `estimated_post_burnin_Reff`: Estimated Reff after burnin period (default: 0.8)
- `adjust_vaccination_coverage`: Adjustment to apply to calculated coverage (default: -0.1)
- `vaccination_bounds_spread`: Range spread around calculated coverage (default: 0.2)
- `reach_target_prop_of_remaining_simulation`: Proportion of post-burnin time to reach target (default: 0.75)
- `digits`: Number of decimal places for rounding (default: 1)

# Returns
- `Tuple{Float64, Float64}`: (min_post_burnin_vaccination_coverage, max_post_burnin_vaccination_coverage)
"""
function calculate_post_burnin_vaccination_range(
        time_parameters,
        state_parameters,
        R_0;
        births_per_k_pop = 27.0,
        target_Reff = 1.0,
        estimated_post_burnin_Reff = 0.8,
        adjust_vaccination_coverage = -0.1,
        vaccination_bounds_spread = 0.2,
        reach_target_prop_of_remaining_simulation = 0.75,
        digits = 1
    )
    time_to_reach_Reff = reach_target_prop_of_remaining_simulation *
        (time_parameters.tlength - Dates.days(time_parameters.burnin))

    N = state_parameters.init_states.N

    max_post_burnin_vaccination = calculate_vaccination_rate_to_achieve_Reff(
        time_to_reach_Reff,
        N * estimated_post_burnin_Reff / R_0,
        N,
        target_Reff,
        R_0,
        calculate_mu(births_per_k_pop)
    )
    max_post_burnin_vaccination_coverage = round(
        clamp(
            floor(max_post_burnin_vaccination; digits = 1) + adjust_vaccination_coverage,
            0.0,
            1.0
        );
        digits = digits
    )
    min_post_burnin_vaccination_coverage = round(
        clamp(
            max_post_burnin_vaccination_coverage - vaccination_bounds_spread,
            0.0,
            1.0
        );
        digits = digits
    )
    return (min_post_burnin_vaccination_coverage, max_post_burnin_vaccination_coverage)
end
