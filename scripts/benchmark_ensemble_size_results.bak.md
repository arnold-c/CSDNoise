  Activating project at `~/Documents/Repos/CSDNoise`
Ensemble Size Benchmark (100 vs 1000 vs 10000 simulations)
======================================================================
Configuration:
  Scenario: Autocovariance, Poisson(1.0x), Test(0.9,0.9)
  Sobol points: 100
  Ensemble sizes: 100, 1000, 10000


Testing with 100 simulations
--------------------------------------------------
Generating ensemble data...
  Generating ensemble dynamics...
  Generating null dynamics...
  Processing outbreak data...
  Data generation: 1.75s
Running multistart optimization...
  Optimization: 8.62s
  Best accuracy: 0.675
  Total time: 10.37s

Testing with 1000 simulations
--------------------------------------------------
Generating ensemble data...
  Generating ensemble dynamics...
  Generating null dynamics...
  Processing outbreak data...
  Data generation: 13.57s
Running multistart optimization...
  Optimization: 20.74s
  Best accuracy: 0.6685
  Total time: 34.31s

Testing with 10000 simulations
--------------------------------------------------
Generating ensemble data...
  Generating ensemble dynamics...
  Generating null dynamics...
  Processing outbreak data...
  Data generation: 134.55s
Running multistart optimization...
  Optimization: 208.84s
  Best accuracy: 0.676
  Total time: 343.39s

COMPARISON ANALYSIS
======================================================================

Accuracy Comparison
----------------------------------------
100 simulations:   0.675
1000 simulations:  0.6685
10000 simulations: 0.676

Improvements:
  1000 vs 100:   -0.0065 (-0.96%)
  10000 vs 1000: 0.0075 (1.13%)
  10000 vs 100:  0.001 (0.16%)

→ 10000 simulations provides the best accuracy

Performance Analysis
----------------------------------------
Total Times:
  100 simulations:   10.37s
  1000 simulations:  34.31s
  10000 simulations: 343.39s

Time per simulation:
  100 simulations:   103.68ms
  1000 simulations:  34.31ms
  10000 simulations: 34.34ms

Scaling factors:
  1000 vs 100:   3.31x
  10000 vs 1000: 10.01x
  10000 vs 100:  33.12x

Efficiency (accuracy/second):
  100 simulations:   0.0651
  1000 simulations:  0.0195
  10000 simulations: 0.002

→ Most efficient: 100 simulations
→ Excellent scaling: 0.33x of linear

Parameter Comparison
----------------------------------------
Threshold Percentile:
  100 sims:   0.987
  1000 sims:  0.988
  10000 sims: 0.99
  Differences:
    1000 vs 100:   0.002
    10000 vs 1000: 0.002
    10000 vs 100:  0.003

Consecutive Thresholds:
  100 sims:   5
  1000 sims:  2
  10000 sims: 2
  Differences:
    1000 vs 100:   3
    10000 vs 1000: 0
    10000 vs 100:  3

→ Parameter estimates are reasonably consistent
→ Parameters are converging with larger ensemble sizes

Results saved to: ensemble_size_comparison_100sobol_2025-08-29_12-02-29.csv

