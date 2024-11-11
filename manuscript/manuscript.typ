#import "template.typ": article

#show: article.with(
  title: "Diagnostic Uncertainty Limits the Potential of Early Warning Signals for Epidemic Transitions",
  header-title: "true",
  authors: (
    "Callum R.K. Arnold": (
      affiliation: ("PSU-Bio", "CIDD"),
      corresponding: "true",
      email: "contact\@callumarnold.com",
    ),
    "Matthew J. Ferrari": (
      affiliation: ("PSU-Bio", "CIDD"),
    ),
  ),
  affiliations: (
    "PSU-Bio": "Department of Biology, Pennsylvania State University, University Park, PA, USA 16802",
    "CIDD": "Center for Infectious Disease Dynamics, Pennsylvania State University, University Park, PA, USA 16802",
  ),
  keywords: ("Rapid-Diagnostic Tests","ELISA","Infectious Disease Surveillance","Outbreak Detection"),
  abstract: [
300 words\
  ],
  author-summary: [
    200 words\
  ],
  word-count: true,
  line-numbers: true
)

#[
#align(center)[#text(size: 16pt)[#smallcaps("Outline")]]
= Introduction

- Infectious disease surveillance has primarily focussed on converting cases into alerts that trigger an outbreak response: reactive in nature
- Growing interest and research on the converting case data into early warning signals: proactive & can signal that a state is at risk of an outbreak in the future
  - Wide range of fields, but primarily originated in climatic research, thinking about ecosystem collapse and ecological systems @schefferEarlywarningSignalsCritical2009 @schefferForeseeingTippingPoints2010 @dakosSlowingEarlyWarning2008 @drakeEarlyWarningSignals2010 @boettigerQuantifyingLimitsDetection2012
- Prior work has demonstrated that EWS metrics are theoretically correlated with a critical transition for infectious disease systems ($R#sub[effective] = 1$) @oreganTheoryEarlyWarning2013 @drakeStatisticsEpidemicTransitions2019 @brettAnticipatingEmergenceInfectious2017 @brettAnticipatingEpidemicTransitions2018
- Much of this work has focussed on identifying metrics that are correlated with transitions (a necessary first step to prove that it is an avenue worth exploring further). Does this using Kendall's Tau @kendallNEWMEASURERANK1938 @brettAnticipatingEpidemicTransitions2018 @brettDetectingCriticalSlowing2020
  - Useful, but doesn't tell us about how far away we are from the transition, or when to act
- To address this, people have used various threshold based approaches:
  - Build distribution of metric values, and when crossed 2 #sym.sigma, trigger flag @drakeEarlyWarningSignals2010
    - Sometimes require multiple consecutive flags for an alert @delecroixPotentialResilienceIndicators2023 @southallHowEarlyCan2022
    - Can be calculated on a weighted composite of multiple metrics @brettDetectingCriticalSlowing2020 @clementsIncludingTraitbasedEarly2016 @drakeEarlyWarningSignals2010 @clementsIndicatorsTransitionsBiological2018
      - Can use statistical model to define threshold that maximizes accuracy
  - Calculate p-value of Kendall's Tau @southallHowEarlyCan2022 @harrisEarlyWarningSignals2020
    - Bootstrap time series to produce null distribution to compare observed against
- While some has addressed under-reporting, haven't addressed uncertainty from imperfect tests
  - both under and over-reporting


= Results

- Kendall's Tau:
  - 
- Thresholds:
  -

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


= Discussion

- Kendall's tau most affected by length of time, not test accuracy
  - Doing it too early degrades performance


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
    - Distribution threshold percentile (P) $in [0.9, 1.0)$
    - Number of consecutive flags (C) to trigger an alert $in [2, 30]$

]<additional-info>


#pagebreak()

= Introduction

= Results

= Discussion
== Limitations and Strengths

= Materials & Methods


#pagebreak()

#set par.line(
  numbering: none
)

#[
= Funding
- #emph[What to put here?]


= Acknowledgements
== Author Contributions
#emph[Conceptualization:] CA, MJF

#emph[Data curation:] MJF, CA

#emph[Formal analysis:] CA, MJF

#emph[Funding acquisition:] ???

#emph[Investigation:] CA, MJF

#emph[Methodology:] CA, MJF

#emph[Project administration:] MJF

#emph[Software:] CA

#emph[Supervision:] MJF

#emph[Validation:] CA, MJF

#emph[Visualization:] CA

#emph[Writing - original draft:] CA

#emph[Writing - review and editing:] all authors.

== Conflicts of Interest and Financial Disclosures
The authors declare no conflicts of interest.

== Data Access, Responsibility, and Analysis
Callum Arnold and Dr. Matthew J. Ferrari had full access to all the data in the study and take responsibility for the integrity of the data and the accuracy of the data analysis. Callum Arnold and Dr. Matthew J. Ferrari (Department of Biology, Pennsylvania State University) conducted the data analysis.

== Data Availability
All code and data for the simulations can be found at #link("https://github.com/arnold-c/CSDNoise")

#pagebreak()

#set bibliography(style: "elsevier-vancouver")

#bibliography("CSD.bib")
]<additional-info>
