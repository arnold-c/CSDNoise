#import "template.typ": article, two_header_table

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
  line-numbers: false
)

= Introduction

Despite sustained advances over decades, infectious diseases still pose a substantial threat to human life, causing an estimated X number of infections, and Y deaths, per annum _*[REF]*_.
For many diseases, effective and affordable vaccines have played a substantial role in reducing this burden, averting 154 million deaths since the introduction of the Expanded Programme on Immunization in 1974 @shattockContributionVaccinationImproved2024.
On the path to elimination, complex non-linear dynamics may cause an increase in the variability of annual incidence, such as the so-called "canonical path" of measles @grahamMeaslesCanonicalPath2019.
As a result, episodic outbreaks become increasingly important to the total burden of disease, and are therefore vital to detect and respond to if elimination is to be reached.
Infectious disease surveillance systems are the mechanism for this action @murrayInfectiousDiseaseSurveillance2016 @DiseaseSurveillance.

Traditional infectious disease surveillance systems are reactive in nature; suspected and laboratory confirmed cases are collated, counted, and if a pre-determined threshold is met or breached, an action is undertaken (e.g., preliminary investigation, or reactive vaccination campaign) _*[REF]*_.
However, due to the exponential trajectory of incidence often observed in the early stages of an outbreak, the reactive nature necessarily results in excess infections that cannot be prevented _*[REF]*_.
To limit the burden of disease, ideally, epidemiologists could utilize the output of a surveillance system (e.g., the trend in cases of a pathogen) to predict the risk of a future outbreak, triggering a _proactive_ action, such as a preventative vaccination campaign.

There has been growing interest in this line of reasoning, with many fields trying to identify and develop early warning signals (EWS) that are predictive of a critical transition @schefferEarlywarningSignalsCritical2009 @schefferForeseeingTippingPoints2010 @dakosSlowingEarlyWarning2008 @drakeEarlyWarningSignals2010 @boettigerQuantifyingLimitsDetection2012.
For infectious diseases, this critical transition occurs when the effective reproduction number crosses the bifurcation threshold $R_"effective" = 1$, observed during both elimination and emergence of disease.
The appeal of an alert system based upon EWS metrics is that they are model-free, only requiring the calculation of summary statistics of a time series.
Prior work has demonstrated that for infectious disease systems, computing EWS metrics on the progression of population susceptibility may be most predictive @drakeStatisticsEpidemicTransitions2019, but collecting this information is often intractable, and utilizing either the incidence or prevalence data has provided similarly useful predictions @drakeStatisticsEpidemicTransitions2019 @brettAnticipatingEpidemicTransitions2018 @southallProspectsDetectingEarly2020.
If an EWS is predictive, critical slowing down theory suggests that the EWS values will drastically change in value as a transition is approached, such as an increase in the variance.
This is the result of a slowed recovery from perturbations to the system @delecroixPotentialResilienceIndicators2023 @dakosSlowingEarlyWarning2008 @schefferEarlywarningSignalsCritical2009, e.g., an imported infection.
Prior work has demonstrated that EWS metrics are theoretically correlated with a critical transition for infectious disease systems, under emergent and extinction conditions @oreganTheoryEarlyWarning2013 @drakeStatisticsEpidemicTransitions2019 @brettAnticipatingEmergenceInfectious2017 @brettAnticipatingEpidemicTransitions2018 @southallProspectsDetectingEarly2020 @drakeMonitoringPathElimination2017.
While identifying EWS that are correlated with a transition is an important first step, there are some important limitations that arise when designing a surveillance system @dablanderOverlappingTimescalesObscure2022 @southallEarlyWarningSignals2021.
Notably, the correlation with a transition provides no information about the expected time until that transition, which is vital for any planning and preventative actions.
To address these shortcomings, various threshold-based and statistical learning based approaches have been developed, with a 2#sym.sigma threshold most commonly employed @southallHowEarlyCan2022 @drakeEarlyWarningSignals2010 @brettDynamicalFootprintsEnable2020 @clementsIncludingTraitbasedEarly2016 @obrienEarlyWarningSignal2021.
This threshold could be calculated from a single metric, or a composite of multiple EWS metrics @drakeEarlyWarningSignals2010 @obrienEarlyWarningSignal2021, with prior working demonstrating that requiring multiple consecutive flags to trigger an alert improves the accuracy in a 'noisy' system by reducing the false positive rate @southallHowEarlyCan2022 @clementsBodySizeShifts2017 @clementsEarlyWarningSignals2019.

Until now, the relatively nascent topic of EWS has only explored 'noise' in the observational process of an outbreak, such as under-reporting and aggregation of case data @brettAnticipatingEpidemicTransitions2018 @brettDetectingCriticalSlowing2020.
Our goal is to characterize the performance of EWS metrics in detecting the risk of disease emergence in a system with imperfect diagnostic tests and background infections that may be misdiagnosed as the target disease and inappropriately tested.
For diseases with non-specific symptoms, e.g., measles and rubella often co-circulate and present clinically similarly, an imperfect diagnostic test will result in false positive and negative cases.
In this paper we show the conditions under which diagnostic uncertainty overwhelms the time series used to calculate EWS summary statistics, limiting the ability to predict epidemic transitions.

= Materials & Methods
== Model Structure

To evaluate the predictive ability of EWS metrics in environments with clinically-compatible background noise that could be mistaken for the target disease and tested with an imperfect diagnostic, we constructed two stochastic compartmental non-age structured Susceptible-Exposed-Infected-Recovered (SEIR) models.
SEIR models were simulated using a Tau-leaping algorithm with a time step of 1 day, with binomial draws so that no jump resulted in negative compartment sizes @gillespieApproximateAcceleratedStochastic2001 @chatterjeeBinomialDistributionBased2005.
We assumed no seasonality in the transmission rate ($beta_t$), and set the latent and infectious periods equal to 10 days and 8 days, respectively, and an $R_0$ equal to 16, approximating measles parameters values @guerraBasicReproductionNumber2017 @gastanaduyMeasles2019.
Demographic parameters (birth and death rates) broadly reflecting those observed in Ghana were selected to evaluate the performance of EWS metrics in a setting where high, yet sub-elimination, vaccination coverage is observed, requiring ongoing vigilance @WHOImmunizationData @masreshaTrackingMeaslesRubella2024.
An initial population of 500,000 individuals was simulated, with commuter-style imports drawn from a Poisson distribution with mean proportional to the size of the population and $R_0$, to maintain a level of endemicity @keelingModelingInfectiousDiseases2008.

We generated a time series of clinically-compatible febrile rash by summing the daily test positive results from the measles incidence time series, and a noise time series.
The noise time series was the result of: independent draws of a Poisson distribution, with mean equal to a multiple (c) of the daily average measles incidence, where $c in {1, 7}$; or from an SEIR time series with rubella-like parameters with additional noise drawn from a Poisson distribution with mean equal to 15% of the daily average of the rubella incidence time series, to account for non-rubella sources of clinically-compatible febrile rash e.g., parvovirus (@tbl-model-parameters) @papadopoulosEstimatesBasicReproduction2022 @RubellaCDCYellow.
Under dynamical noise simulations, the vaccination rate at birth was selected to produce equivalent magnitudes of daily average noise incidence as observed in the Poisson-like noise simulations (10.20% and 87.34%).
Throughout the rest of the manuscript, these will be referred to as low and high Poisson/dynamical noise scenarios, accordingly.
Each day, all clinically-compatible febrile rash cases (that is, both the measles and noise time series) were tested using one of the following diagnostic tests, producing a time series of test positive cases.

- A perfect test with 100% sensitivity and specificity, with a 0-day test result delay. This was chosen to reflect the best-case scenario that the imperfect diagnostic-based alert scenarios could be compared against.
- An RDT equivalent, imperfect diagnostic with a 0-day test result delay and sensitivity and specificity equal to either 99%, 98%, 97%, 95%, 90%, or 80%.

#let import_rate = $(1.06*μ*R_0)/(√(N))$

#figure(
  table(
    columns: 4,
    [Parameters],[Measles - Emergent],[Measles - Null],[Dynamical noise],
    [R0],table.cell(colspan: 2, align: center, "16"),[5],
    [Latent period (s)],table.cell(colspan: 2, align: center, "10 days"),[7 days],
    [Infectious period (g)],table.cell(colspan: 2, align: center, "8 days"),[14 days],
    [Vaccination rate at birth \ during burn-in period (r#sub[i])],table.cell(colspan: 2, align: center, "Unif (92.69%, 100%)"),[10.20%, 83.74%],
    [Vaccination rate at birth \ after burn-in period (r#sub[e])],[Unif (60%, 80%)],[Unif (92.69%, 100%)],[10.20%, 83.74%],
    [Birth/death rate (m)],table.cell(colspan: 3, align: center, "27 per 1000 per annum"),
    [Importation rate], table.cell(colspan: 3, align: center, $(1.06*μ*R_0)/(√(N))$),
    [Population size (N)], table.cell(colspan: 3, align: center, "500,000"),
    [Initial proportion susceptible], table.cell(colspan: 3, align: center, "0.05"),
    [Initial proportion exposed], table.cell(colspan: 3, align: center, "0.0"),
    [Initial proportion infected], table.cell(colspan: 3, align: center, "0.0"),
    [Initial proportion recovered], table.cell(colspan: 3, align: center, "0.95"),
  ),
  caption: [Compartmental model parameters],
)
<tbl-model-parameters>

To evaluate the performance of the EWS metrics at predicting the approach to the critical threshold ($R_"effective" = 1$) from below, it is essential to simulate both emergent and null scenarios.
For both emergent and null scenarios, we generated 100 time series.
All measles simulation incorporated a 5 year burn-in period to produce sufficient data for calculation of the EWS metrics upon aggregation, as well as to produce greater variation in the trajectory of $R_"effective"$.
For each time series, the vaccination rate at birth during the burn-in period was sampled from a Uniform distribution between 92.69% and 100% coverage.
These bounds were selected to ensure the maximum value of $R_"effective"$ that could be reached within 10 years (twice the length of the burn-in period) was 0.9.
We simulated emergent scenarios by allowing the proportion of the population that is susceptible to grow, through lowering the vaccination rate at birth after completion of the burn-in period.
For each emergent time series, the vaccination rate at birth was independently drawn from a Uniform distribution between 60% and 80% coverage, allowing the rate of growth in $R_"effective"$ to vary between emergent time series.
For each null time series, the vaccination rate at birth was set to the coverage sampled during the burn-in period, ensuring $R_"effective"$ would not cross the critical threshold within the scope of the simulation, though it may grow slowly.
Each of the 100 emergent and 100 null time series are paired during the pre-processing steps i.e., up until the completion of the burn-in period, paired emergent and null simulations share the same vaccination rate at birth, and they are both truncated to identical lengths (the time step when $R_"effective" = 1$ in that pair's emergent simulation).

All simulations and analysis was completed in Julia version 1.10.5 @bezansonJuliaFreshApproach2017, with all code stored at #link("https://github.com/arnold-c/CSDNoise").

== Computing & Evaluating EWS
Each set of null and emergent time series are aggregated by month and numerical estimates of the EWS metrics were then calculated on the aggregated time series, de-trended using backwards-facing moving averages with bandwidth $b = 52$ weeks.
For example, the EWS metric, the mean, is given by the expectation:

$$$
hat(mu)_t &= sum_(s = t-b delta)^(s = t) X_s / b
$$$

where $X_s$ represents the aggregated incidence at time point (month) $s$, and $delta = 1$ time step (in the simulation results presented, 1 month).
At the beginning of the time series when $t < b$, $b$ is set equal to $t$.

In this paper we evaluate the performance of the following EWS metrics: the mean, variance, coefficient of variation, index of dispersion, skewness, kurtosis, autocovariance, and autocorrelation at lag-1, as they have been some evidence in the literature that they are correlated or predictive of disease emergence @brettAnticipatingEpidemicTransitions2018 @drakeStatisticsEpidemicTransitions2019 @southallEarlyWarningSignals2021 @southallEarlyWarningSignals2021 @brettDetectingCriticalSlowing2020, or commonly evaluated in critical slowing down literature at large.
Many of the EWS metrics rely on the prior computation of others e.g., the variance requires the calculation of the detrended mean, so are computed in an iterative pattern.
The full list of numerical formulas for each EWS metric can be found in @tbl-ews-formulas.

#let table_math(content) = table.cell(text(size: 16pt, content))

#figure(
  table(
    columns: 2, align: center, inset: 10pt,
    [EWS Metric],[Formula],
    [Mean ($hat(mu)_t$)],table_math[$sum_(s = t-b delta)^(s = t) X_s / b$],
    [Variance ($hat(sigma)^2_t$)], table_math[$sum_(s = t-b delta)^(s=t) (X_s - hat(mu)_s)^2 / b$],
    [Coefficient of Variation ($hat("CV")_t$) ], table_math[$hat(sigma)_t / hat(mu)_t$],
    [Index of Dispersion ($hat("IoD")_t$) ], table_math[$hat(sigma)^2_t / hat(mu)_t$],
    [Skewness ($hat("Skew")_t$)], table_math[$1/(hat(sigma)^3_t) sum_(s = t-b delta)^(s = t) (X_s - hat(mu)_s)^3 / b$],
    [Kurtosis ($hat("Kurt")_t$)], table_math[$1/(hat(sigma)^4_t) sum_(s = t-b delta)^(s = t) (X_s - hat(mu)_s)^4 / b$],
    [Autocovariance ($hat("ACov")_t$)], table_math[$sum_(s = t-b delta)^(s=t) ((X_s - hat(mu)_s)(X_(s-delta) - hat(mu)_(s-delta))) / b$],
    [Autocorrelation lag-1 ($hat("AC-1")_t$)], table_math[$hat("ACov")_t / (hat(sigma)_t hat(sigma)_(t-delta)) $],
  ),
  caption: [Numerical computations for EWS metrics, where $delta = 1$ time step, $b = 52$ weeks]
)
<tbl-ews-formulas>

Once the EWS metrics have been computed, the correlation within emergent time series is computed using Kendall's Tau-B, signifying if an EWS metric consistently increases (or decreases) in magnitude throughout the time series @kendallTREATMENTTIESRANKING1945 @knightComputerMethodCalculating1966.
Kendall's Tau is computed on two lengths of time series: from the beginning of the simulation until the critical transition is met, and from the completion of the burn-in period until the critical transition.
To evaluate the strength of the correlation, we use the area under the receiver operator curve (AUC), as described in prior papers @brettAnticipatingEpidemicTransitions2018 @southallEarlyWarningSignals2021 @southallProspectsDetectingEarly2020.
Briefly, the calculation of the AUC compares whether the distributions of Kendall's Tau differ substantially between emergent and null simulations for a given alert scenario and EWS metric.
AUC is calculated using the rank order of the EWS metrics for both emergent and null time series using the equation @flachROCAnalysis2016

$$$
(r_"null" - n_"null" (n_"null" + 1) \/ 2) / (n_"emergent" n_"null")
$$$

where $r_"null"$ equals the sum of ranks for the null time series, and $n_"null"$ and $n_"emergent"$ refer to the number of null and emergent simulations, respectively.
An AUC of 0.5 indicates the EWS is similarly correlated with both emergent and null time series, offering no benefit; values > 0.5 indicate a positive correlation with emergent time series, and < 0.5 indicates the EWS metric is negatively correlated with the emergent simulations.
AUC values are commonly transformed as $|"AUC" - 0.5|$ to highlight the strength of the correlation with emergence, with values close to 0 exhibiting poor performance, and a value of 0.5 indicating perfect correlation @brettAnticipatingEpidemicTransitions2018.

The primary mode of evaluation for the EWS metrics relies on computing and triggering an alert based upon a set of conditions.
The alert scenario is defined as the combination of diagnostic test, noise structure and magnitude, and EWS metric.
The combination of the alert scenario and EWS alert hyperparameters (the percentile threshold value of the long-running metric distribution that must be exceeded to create a flag, and the number of consecutive flags required to trigger an alert), produce distinct counts of emergent and null time series that result in an alert.
The sensitivity of the system is defined as the proportion of the emergent simulations that result in an alert, and the specificity is the proportion of the null simulations that do not result in an alert.
Taking the mean of the sensitivity and specificity produces the accuracy of the system.
For each alert scenario, a grid search over the EWS hyperparameters (percentile threshold $in [0.5, 0.99]$, consecutive flags $in [2, 30]$) is performed to identify the set of EWS hyperparameters that maximizes alert accuracy for a given alert scenario.
If multiple hyperparameter combinations produce identical alert system accuracies, the combination with the highest specificity is selected.
After the optimal EWS hyperparameters have been selected, the accuracy of each EWS metric are compared across alert scenarios, at their respective maximal values.
Finally, the speed and timing of detection relative to the critical threshold is evaluated using Kaplan-Meier survival estimates @clarkSurvivalAnalysisPart2003.

= Results
== Correlation with Emergence

The strength and direction of the raw correlation between EWS metrics and the approach to the critical threshold in emergent time series is strongly dependent upon the length of the time series evaluated than the characteristics of the diagnostic test (@tbl-tau-ranking-perfect-test).
However, when calculating AUC to normalize the correlation in the emergent time series against the correlation observed in null simulations, this affect disappears (@tbl-tau-ranking-perfect-test).
Consistent with previous studies, the autocovariance, variance, mean, and index of dispersion show the strongest correlations with emergence ($|"AUC"-0.5| = 0.2, 0.2, 0.18$, evaluated after the burn-in period, respectively) @brettDetectingCriticalSlowing2020 @brettAnticipatingEpidemicTransitions2018.

#let perfect_tau_auc_table = csv("./manuscript_files/tables/perfect-test_tau-auc.csv")

#figure(
    two_header_table(
    columns: 5,
    table.cell(rowspan: 2, align: horizon)[Rank], table.cell(colspan: 2)[Tau], table.cell(colspan: 2)[|AUC - 0.5|],
    [Full Time Series], [After Burn-In Period],
    [Full Time Series], [After Burn-In Period],
    ..perfect_tau_auc_table.slice(1).flatten()
  ),
  caption: [The ranking and mean value of Kendall's Tau computed on emergent time series, and the $|"AUC" - 0.5|$ for each metric. The values are computed on the full time series, and the subset from after the completion of the burn-in period, with a perfect test]
)
<tbl-tau-ranking-perfect-test>

With an imperfect diagnostic test, the noise structure was more influential to correlation with emergence than the noise magnitude (@tbl-auc-magnitude-ranking-rdt-comparison).
For an RDT-equivalent test with 90% sensitivity and specificity, with a 0-day turnaround time, the correlation between all EWS metrics and emergence was relatively unaffected by the magnitude of Poisson noise.
The top four metrics with a perfect diagnostics (autocovariance, variance, mean, and index of dispersion) maintained their positions as the most strongly correlated metrics.
Under low levels of Poisson noise, autocorrelation was more strongly associated with emergence when calculated on the test positive time series resulting from and RDT than from a perfect diagnostic test ($|"AUC"-0.5| = 0.17 "vs." 0.12$, respectively).
Under high levels of Poisson noise, the association of coefficient of variation with emergence also increased ($|"AUC"-0.5| = 0.17$ for an RDT vs. 0.11 for a perfect diagnostic).
When simulations included rubella-like SEIR dynamical noise, the correlation of all metrics decreased (@tbl-auc-magnitude-ranking-rdt-comparison), and was exacerbated at a higher magnitude of noise.
With low levels of dynamical noise, the autocovariance, variance, and mean maintained some correlation with emergence, albeit at a lower level than observed in Poisson noise scenarios ($|"AUC" - 0.5| = 0.16, 0.14, "and" 0.13$, respectively).
However, when the magnitude of dynamical noise was increased, these correlations disappeared, with all EWS metrics exhibiting $|"AUC"-0.5| lt.eq 0.05$.

Among the EWS metrics that were correlated with emergence, most increased in value as the tipping point neared (AUC > 0.5) when computed on the test positive time series resulting from perfect diagnostic tests or an RDT with 90% sensitivity and specificity, with the exception of: the coefficient of variation, with a perfect diagnostic test, and an RDT in low and high dynamical noise ($"AUC" = 0.39, 0.45, "and" 0.48$, respectively); kurtosis with an RDT under low Poisson noise ($"AUC" = 0.45$); and autocorrelation with an RDT under high dynamical noise ($"AUC" = 0.49$) (Supplemental Table 2).

#let auc_magnitude_comparison_table = csv("./manuscript_files/tables/auc-magnitude-comparison.csv")

#figure(
  two_header_table(
    columns: 6,
    table.cell(rowspan: 2, align: horizon)[Rank], [Perfect Test], table.cell(colspan: 4)[90% Sensitive & 90% Specific RDT],
    ..auc_magnitude_comparison_table.flatten().slice(1)
  ),
  caption: [$|"AUC" - 0.5|$ for EWS metrics, ranked in descending order of magnitude, computed on the subset of the emergent time series after the burn-in period, for a perfect test and an RDT with 90% sensitivity and 90% specificity, under high and low Poisson and dynamical noise systems]
)
<tbl-auc-magnitude-ranking-rdt-comparison>

== Predictive Ability

Each alert scenario (the combination of diagnostic test, noise structure and magnitude, and EWS metric) produced its optimal accuracy with a different combination of EWS hyperparameters (the percentile threshold of the long-running metric distribution to be exceeded to return a flag, and the number of consecutive flags required to trigger an alert) (Supplemental Figures 9-12).
At their respective maximal accuracies, the relative ranking of the EWS metrics computed with a perfect diagnostic test remained consistent to the ranking based upon $|"AUC" - 0.5|$: Mean (accuracy = 0.72), variance (0.72), autocovariance (0.7), index of dispersion (0.63), autocorrelation (0.62), skewness (0.6), kurtosis (0.58), and coefficient of variation (0.5) (Supplemental Table 3).

When EWS metrics were computed on time series generated from RDTs, each metric's accuracy generally remained constant, with a few notable exceptions (@fig-best-accuracy-line-plot, @fig-worse-accuracy-line-plot).
For the 4 most correlated metrics (autocovariance, variance, mean, and index of dispersion), the accuracy observed with imperfect diagnostic tests was primarily affected by high levels of dynamical noise e.g., the maximal accuracies achieved for autocovariance with a perfect diagnostic test and an 80% sensitive and specific RDT were 0.72 and 0.52, respectively (@fig-best-accuracy-line-plot, Supplemental Figure 12).
The primary exception to this trend is the index of dispersion, which also demonstrated higher accuracy with imperfect diagnostics under high levels of Poisson noise, increasing from 0.63 to 0.70 across the same two tests (@fig-best-accuracy-line-plot, Supplemental Figure 10).
Among the 4 least correlated EWS metrics (autocorrelation, coefficient of variation, skewness, and kurtosis), the use of RDTs generally did not change the accuracy of the alert systems, relative to a perfect test (@fig-worse-accuracy-line-plot).
However, the autocorrelation experienced a slight increase in accuracy with imperfect diagnostic tests and Poisson noise (0.67 for an 80% sensitive and specific RDT in 7x Poisson noise, vs. 0.62 with a perfect test) (Supplemental Figure 10), and slight decreases under dynamical noise (0.54 in 7x dynamical noise) (Supplemental Figure 12).
More surprisingly, the accuracy drastically improved for the coefficient of variation under high levels of Poisson noise with imperfect diagnostics (accuracy of 0.7 vs 0.5, for an 80% sensitive and specific RDT and perfect test, respectively) (@fig-worse-accuracy-line-plot, Supplemental Figure 10).

#figure(
  image("./manuscript_files/plots/accuracy-line-plot.svg"),
  caption: [The change in alert accuracy for the most correlated EWS metrics under increasing diagnostic uncertainty, and low and high levels of Poisson or dynamical noise]
)
<fig-best-accuracy-line-plot>

#figure(
  image("./supplemental_files/plots/accuracy-line-plot.svg"),
  caption: [The change in alert accuracy for less correlated EWS metrics under increasing diagnostic uncertainty, and low and high levels of Poisson or dynamical noise]
)
<fig-worse-accuracy-line-plot>

In addition to the accuracy of the EWS-based alert system, it is important to evaluate the speed and relative false positive and negative alert rates.
With a perfect diagnostic test, the 4 most correlated metrics exhibited steeper and more sustained declines in survival than the 4 least correlated metrics i.e., they alerted more frequently, and earlier in the time series (@fig-autocovariance-survival, Supplemental Figures 14-20).
For the autocovariance, the use of an RDT resulted in slightly earlier and more frequent alerts under Poisson noise, but a more specific system with fewer and later alerts under dynamical noise (@fig-autocovariance-survival).
In each of the other 3 most correlated EWS metric survival curves, the use of an RDT consistently led to more specific alert systems, regardless of noise structure or magnitude (Supplemental Figures 14-16).
For the 4 least correlated metrics, the same pattern was generally observed, with the marked exception of the coefficient of variation, which increases in alert sensitivity and discriminatory ability for Poisson noise (Supplemental Figures 17-20).

#figure(
  image("./manuscript_files/plots/survival/survival_ews-autocovariance.svg"),
  caption: [Survival curves for the autocovariance EWS metric computed on emergent and null simulations, with a perfect test and an RDT equivalent with 90% sensitivity and specificity. The histogram depicts the times when the tipping point is reached ($R_"effective" = 1$) under the emergent simulation, right-truncating the curves.]
)
<fig-autocovariance-survival>

= Discussion

Corroborating previous findings, we find that in perfectly observed systems with reporting uncertainty (monthly case aggregation), the autocovariance, variance, and mean are all well correlated with the emergence of outbreaks.
As uncertainty rises, through the use of an imperfect diagnostic tests for laboratory confirmation of clinically-compatible cases, these correlations are maintained in all scenarios, except when there is a large amount of dynamical noise.
These patterns of correlation translate to alert accuracy when incorporating a threshold-based approach, indicating the potential benefit if implemented in a real-world situation.
For example, in Botswana, the estimated incidence rates of rubella and measles for 2023 are approximately equal, suggesting EWS calculated with imperfect diagnostic tests may provide accurate and actionable results @masreshaTrackingMeaslesRubella2024.
But in Guinea Bissau, the estimated incidence rate for rubella was approximately 9 times that of measles, highlighting a region where any diagnostic test is unlikely to provide a meaningfully beneficial proactive alert system @masreshaTrackingMeaslesRubella2024.

When evaluating the ability for the EWS metrics to accurately discriminate between emergent and null simulations, it is import to contextualize the results with the system's relative speed and specificity.
Alert systems necessarily make compromises in their design: improvements to speed generally come at the cost of increased numbers of false alerts _*[REF]*_.
Depending on the context, it may be desirable to place a greater weight in preference/penalty for one of these axes; in scenarios where the expected cost to launch a preliminary investigation is low relative to the unaverted DALYs resulting from incorrect inaction in an overly specific system, higher false alert rates may be acceptable.
This analysis provides a framework to explicitly explore these trade-offs through the comparison of survival curves.
A larger separation at the end of the time series between the emergent and null simulation lines indicates higher accuracy, as there is a greater difference in the true positive and false positive rates.
Faster declines indicate a (relatively) more sensitive alert system with more advanced warning of emergence.
Under the simulation constraints placed here (i.e., equal weighting to speed and specificity to the alert system's accuracy, with the results of the more specific system being presented in cases where multiple thresholds hyperparameters provide the same accuracy), generally, the use of RDTs does not increase the speed of the warning for the EWS metrics that are predictive of emergence.
This is likely a consequence of imperfect diagnostic tests producing more false positive cases, which, without appropriate penalization, would otherwise lead to high false alert rates under the same EWS hyperparameters.
Adjusting the relative weighting of alert sensitivity and specificity to accuracy, as well as the accuracy tiebreaker preference, would allow for an exploration of alternative scenarios.

For EWS metrics to reflect the underlying dynamics of critical slowing down, careful detrending of the data is required @gamadessavreProblemDetrendingWhen2019 @dakosSlowingEarlyWarning2008 @lentonEarlyWarningClimate2012.
Our analysis utilizes a backward-facing uniform moving average to detrend the data: this was chosen as it can easily be intuited and implemented in a surveillance system.
However, it has previously been stated that insufficient detrending may lead spurious patterns that do not arise from a system's dynamical response @dakosSlowingEarlyWarning2008, and the association of an EWS with the approach to a tipping point may be sensitive to the bandwidth size selected, although the sensitivity to detrending varies by EWS metric @gamadessavreProblemDetrendingWhen2019 @lentonEarlyWarningClimate2012.
While it may be preferred to detrend using the mean over multiple realizations, this is clearly not possible in a real-world situation.
Future work could explore the effects of different detrending methods (e.g., Gaussian weighting moving average, smaller and/or larger bandwidths) on the effectiveness of EWS metrics in systems with diagnostic uncertainty.

Despite being relatively well-established in areas of study such as ecology, ecosystem collapse, and climate science @drakeEarlyWarningSignals2010 @boettigerQuantifyingLimitsDetection2012 @dakosSlowingEarlyWarning2008 @schefferEarlywarningSignalsCritical2009 @obrienEarlyWarningSignal2021 @carpenterEarlyWarningsRegime2011 @dudneyElusiveSearchTipping2020, the exploration and development of EWS for infectious disease systems is in its relative infancy.
Until recently, a large proportion of the prior work in the area has been to establish the existence of these metrics that theoretically could be used in such a system @drakeMonitoringPathElimination2017 @drakeStatisticsEpidemicTransitions2019 @oreganTheoryEarlyWarning2013.
While this is a crucial first step, for use in an proactive outbreak alert system, EWS metrics must be able to provide advance warning of the approach to the tipping point $R_"effective" = 1$.
Correlations alone are not sufficient to indicate when and what actions must be taken.
To address this, there is a growing body of work that seeks to evaluate the use of various threshold and risk-based approaches within infectious disease systems @southallHowEarlyCan2022 @southallEarlyWarningSignals2021 @brettDetectingCriticalSlowing2020 @brettDynamicalFootprintsEnable2020.
Our work expands upon these efforts, characterizing the limits of predictability for EWS metrics in systems with diagnostic uncertainty and background noise.

#pagebreak()

// #set par.line(
//   numbering: none
// )

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
