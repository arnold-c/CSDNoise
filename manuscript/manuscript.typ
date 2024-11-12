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
The appeal of an alert system based upon EWS metrics is that they are model-free, only requiring the calculation of summary statistics of a time series; in the case of infectious diseases, either the incidence or prevalence data @southallProspectsDetectingEarly2020.
If an EWS is predictive, critical slowing down theory suggests that the EWS values will drastically change in value as a transition is approached, such as an increase in the variance.
This is the result of a slowed recovery from perturbations to the system @delecroixPotentialResilienceIndicators2023 @dakosSlowingEarlyWarning2008 @schefferEarlywarningSignalsCritical2009, e.g., an imported infection.
Prior work has demonstrated that EWS metrics are theoretically correlated with a critical transition for infectious disease systems, under emergent and extinction conditions @oreganTheoryEarlyWarning2013 @drakeStatisticsEpidemicTransitions2019 @brettAnticipatingEmergenceInfectious2017 @brettAnticipatingEpidemicTransitions2018 @southallProspectsDetectingEarly2020 @drakeMonitoringPathElimination2017.
While identifying EWS that are correlated with a transition is an important first step, there are some important limitations that arise when designing a surveillance system, as noted by Southall _et al._ @southallEarlyWarningSignals2021.
Notably, the correlation with a transition provides no information about the expected time until that transition, which is vital for any planning and preventative actions.
To address these shortcomings, various threshold-based and statistical learning based approaches have been developed, with a 2#sym.sigma threshold most commonly employed @southallHowEarlyCan2022 @drakeEarlyWarningSignals2010 @brettDynamicalFootprintsEnable2020 @clementsIncludingTraitbasedEarly2016 @obrienEarlyWarningSignal2021.
This threshold could be calculated from a single metric, or a composite of multiple EWS metrics @drakeEarlyWarningSignals2010 @obrienEarlyWarningSignal2021, with prior working demonstrating that requiring multiple consecutive flags to trigger an alert improves the accuracy in a 'noisy' system by reducing the false positive rate @southallHowEarlyCan2022 @clementsBodySizeShifts2017 @clementsEarlyWarningSignals2019.

Until now, the relatively nascent topic of EWS has only explored 'noise' in the observational process of an outbreak, such as under-reporting and aggregation of case data @brettAnticipatingEpidemicTransitions2018 @brettDetectingCriticalSlowing2020.
Our goal is to characterize the performance of EWS metrics in detecting the risk of disease emergence in a system with imperfect diagnostic tests and background infections that may be misdiagnosed as the target disease and inappropriately tested.
For diseases with non-specific symptoms, e.g., measles and rubella often co-circulate and present clinically similarly, an imperfect diagnostic test will result in false positive and negative cases.
In this paper we show the conditions under which diagnostic uncertainty overwhelms the time series used to calculate EWS summary statistics, limiting the ability to predict epidemic transitions.


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
