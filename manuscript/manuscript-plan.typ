#set table(
  fill: (x, y) => {
    if y == 0 {gray}
  }
)

#show table.cell: it => {
  if it.y == 0 {
    strong(it)
  } else {
    it
  }
}

#align(center)[#text(size: 16pt)[#smallcaps("Outline")]]
= Introduction

- Infectious disease surveillance has primarily focussed on converting cases into alerts that trigger an outbreak response: reactive in nature
- Growing interest and research on the converting case data into early warning signals: proactive & can signal that a state is at risk of an outbreak in the future
  - Wide range of fields, but primarily originated in climatic research, thinking about ecosystem collapse and ecological systems @schefferEarlywarningSignalsCritical2009 @schefferForeseeingTippingPoints2010 @dakosSlowingEarlyWarning2008 @drakeEarlyWarningSignals2010 @boettigerQuantifyingLimitsDetection2012
- Prior work has demonstrated that EWS metrics are theoretically correlated with a critical transition for infectious disease systems ($R#sub[effective] = 1$) @oreganTheoryEarlyWarning2013 @drakeStatisticsEpidemicTransitions2019 @brettAnticipatingEmergenceInfectious2017 @brettAnticipatingEpidemicTransitions2018
- Much of this work has focussed on identifying metrics that are correlated with transitions (a necessary first step to prove that it is an avenue worth exploring further). Does this using Kendall's Tau @kendallNEWMEASURERANK1938 @brettAnticipatingEpidemicTransitions2018 @brettDetectingCriticalSlowing2020
  - Useful, but doesn't tell us about how far away we are from the transition, or when to act
  - Less predictive when only including data up to critical transition (i.e., not including cases from early stages of outbreak) @southallHowEarlyCan2022
- To address this, people have used various threshold based approaches:
  - Build distribution of metric values, and when crossed 2 #sym.sigma, trigger flag @drakeEarlyWarningSignals2010 @brettDynamicalFootprintsEnable2020
    - Sometimes require multiple consecutive flags for an alert @delecroixPotentialResilienceIndicators2023 @southallHowEarlyCan2022
    - Can be calculated on a weighted composite of multiple metrics @brettDetectingCriticalSlowing2020 @clementsIncludingTraitbasedEarly2016 @drakeEarlyWarningSignals2010 @clementsIndicatorsTransitionsBiological2018
      - Can use statistical model to define threshold that maximizes accuracy
  - Calculate p-value of Kendall's Tau @southallHowEarlyCan2022 @harrisEarlyWarningSignals2020
    - Bootstrap time series to produce null distribution to compare observed against
- While some has addressed under-reporting, haven't addressed uncertainty from imperfect tests
  - both under and over-reporting
- 


= Results

- Kendall's Tau:
  - Ranked order for perfect test - descending Tau (correlation)
    - #table(
      columns: 2,
      [Full Length],[After 5yr Burnin],
      [Autocovariance (+)], [Variance (+++)],
      [Variance (+)], [Autocovariance (+++)],
      [Mean (+)],[ Kurtosis (+++)],
      [Iod (+)], [Iod (+++)],
      [Autocorrelation (+)], [Mean (+++)],
      [Kurtosis (-)], [Skewness (++)],
      [CoV (-)], [Autocorrelation (+)],
      [Skewness (-)], [CoV (-)]
    )
  - Full length:
    - As test accuracy decreases, tau generally decreases with high Poisson noise, but not much change when noise is low
    - For dynamical noise, much more drastic decrease in tau, particularly with high noise
  - After burnin:
    - No consistency to change in tau as test accuracy decreases, particularly for high values of (either) noise
- Thresholds:
  - For a perfect test, detection isn't affected by noise structure or magnitude
    - Variance and mean both produce the highest accuracy (73%), with autocovariance close behind (71%)
    - Mean is more specific (71% vs 65%), with a longer delay (more consecutive thresholds required)
    - 7 day aggregation
      - #table(
        columns: 3,
        [Full Length Tau],[After 5yr Burnin Tau], [Threshold Accuracy (5yr burnin)],
        [Autocovariance (0.24)], [Variance (0.56)], [Variance (0.73)],
        [Variance (0.23)], [Autocovariance (0.53)], [Mean (0.73)],
        [Mean (0.22)], [Kurtosis (0.52)], [Autocovariance (0.71)],
        [Iod (0.21)], [Iod (0.52)], [Iod (0.66)],
        [Autocorrelation (0.20)], [Mean (0.44)], [Autocorrelation (0.62)],
        [Kurtosis (-0.12)], [Skewness (0.36)], [Skewness (0.52)],
        [CoV (-0.13)], [Autocorrelation (0.24)], [Kurtosis (0.5)],
        [Skewness (-0.14)], [CoV (-0.15)], [CoV (0.5)]
      )
    - 28 day aggregation
      - #table(
        columns: 3,
        [Full Length Tau],[After 5yr Burnin Tau], [Threshold Accuracy (5yr burnin)],
        [Iod (0.26)], [Autocovariance (0.78)], [Mean (0.72)],
        [Autocovariance (0.25)], [Variance (0.69)], [Variance (0.71)],
        [Variance (0.25)], [Autocorrelation (0.68)], [Autocovariance (0.70)],
        [Mean (0.22)], [Iod (0.67)], [Iod (0.63)],
        [Autocorrelation (0.17)], [Mean (0.47)], [Autocorrelation (0.62)],
        [CoV (-0.01)], [Kurtosis (0.40)], [Skewness (0.59)],
        [Kurtosis (-0.04)], [Skewness (0.20)], [Kurtosis (0.58)],
        [Skewness (-0.08)], [CoV (-0.01)], [CoV (0.5)]
      )
    - Variance and autocovariance
  - Poisson noise:


= Discussion

- As in other studies, autocovariance and variance seem to be robust under the assumption of a perfect test (i.e., no false positive test results) @brettAnticipatingEpidemicTransitions2018
- Kendall's tau most affected by length of time, not test accuracy
  -  Consistent with prior work that illustrates that window size is important @southallEarlyWarningSignals2021 @lentonEarlyWarningClimate2012 @kaurAnticipatingNovelCoronavirus2020
  - Doing it too early degrades correlation, so if it is going to be used to evaluate CSD, choice of start time is critical
  - Imperfect tests exacerbate effects and demonstrate no consistency in change of correlation
  - No consistency in the pattern of the metrics, particularly across noise structures and length of time series e.g.,
    - Absolute magnitude tends to increase with a shorter time series, but not in a consistent direction
  - In an outbreak situation, would not be able to use correlation as a strong sign that may be approaching critical transition, as no confidence whether the correlation is a result of an inaccurate test, and as more data is collected, even under the assumption of a perfect test and there is an eventual approach to the critical transition, the correlation will decline
    - May require a rolling window, rather than an expanding window
- Weekly vs monthly aggregation makes little difference to the qualitative results, but does allow for marginally higher accuracies with imperfect tests

== Limitations and Strengths
= Materials & Methods

- Simulated measles with parameters ...
- Noise simulated as Poisson or dynamical
  - Dynamical with rubella-like parameters
- Diagnostic tests applied to produce time series of test positives
  - Perfect test
  - 90/90 RDT
  - 80/80 RDT
- 100 simulations
  - Create paired simulations where Reff crosses 1, and null
    - Paired null simulations use same end point
  - Simulations have vaccination burn-in period of 5 years
    - Between 92.69% and 100% at birth
    - Ensures Reff won't cross threshold until 10 years
  - Vaccination rates identical between null and example time series for burn-in periods (5 years)
    - Also simulated for a burn-in of 50 days (see supplement)
    - In null simulations, vaccination rate set to same as in burnin period
- EWS use backward looking method
  - Calculated on test positive time series
  - Bandwidth of 52 weeks worth of data
  - Aggregate either weekly or monthly data

$$$
hat(mu)_t &= sum_(s = t-(2b-1) delta)^(t) X_s / (2b - 1)\
hat(sigma)^2_t &= sum_(s = t-(2b-1) delta)^(t) (X_s - hat(mu)_s)^2 / (2b - 1)
$$$

- Thresholds:
  - Build distribution of EWS metric during burn in period (5 years)
  - If EWS at time $t$ exceeds percentile (P) of distribution until $t-1$, considered a flag
    - where $t gt.eq 5$ years
  - Calculate sensitivity, specificity, and accuracy
    - Sensitivity: % of outbreak series that flag
    - Specificity: 100 - % of null series that flag
    - Accuracy: mean(sens + spec)
- Optimal threshold parameters:
  - Select combinations that produce the highest accuracy:
    - Distribution threshold percentile (P) $in [0.5, 1.0)$
    - Number of consecutive flags (C) to trigger an alert $in [2, 30]$
    - Multiple combinations of hyperparameters may produce the same accuracy
      - Show the results of the hyperparameters that are the most specific
        - Combination that is the most sensitive (fastest) shown in supplement

= Results Plots
== Tau Heatmaps
=== Full Length

#figure(
  image("manuscript_files/plots/tau-heatmaps/full-length/tau-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau-heatmaps/full-length/tau-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x noise]
)

#figure(
  image("manuscript_files/plots/tau-heatmaps/full-length/tau-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau-heatmaps/full-length/tau-heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x noise]
)

=== After 5yr Burn in

#figure(
  image("manuscript_files/plots/tau-heatmaps/after-burnin/tau-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau-heatmaps/after-burnin/tau-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x noise]
)

#figure(
  image("manuscript_files/plots/tau-heatmaps/after-burnin/tau-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau-heatmaps/after-burnin/tau-heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x noise]
)


== Optimal Threshold Accuracies

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x noise]
)

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x noise]
)

#pagebreak()

#set bibliography(style: "elsevier-vancouver")

#bibliography("CSD.bib")

