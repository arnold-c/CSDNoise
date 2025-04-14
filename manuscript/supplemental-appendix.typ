#import "template.typ": *

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
)

#show figure.where(kind: table): set figure(supplement: [Supplemental Table])
#show figure.where(kind: image): set figure(supplement: [Supplemental Figure])

== Figures

#figure(
  image("./supplemental_files/plots/accuracy-line-plot.svg"),
  caption: [The change in alert accuracy for the least correlated EWS metrics under increasing diagnostic uncertainty, and low and high levels of static or dynamical noise. Low noise refers to simulations where the average incidence of noise is equal to the average incidence of measles. High noise refers to simulations where the average incidence of noise is equal to 7 times the average incidence of measles. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig-worst-accuracy-line-plot>

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_poisson_1.0x.svg"),
  caption: [The strength of the correlation ($|"AUC" - 0.5|$) for each EWS metric with emergence, at low levels of static noise, for diagnostic tests of varying accuracy, and was computed after the completion of the burn-in period. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-auc-mag-heatmap-poisson-1x>

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_poisson_7.0x.svg"),
  caption: [The strength of the correlation ($|"AUC" - 0.5|$) for each EWS metric with emergence, at high levels of static noise, for diagnostic tests of varying accuracy, and was computed after the completion of the burn-in period. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-auc-mag-heatmap-poisson-7x>

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_dynamical_0.8734.svg"),
  caption: [The strength of the correlation ($|"AUC" - 0.5|$) for each EWS metric with emergence, at low levels of dynamical noise, for diagnostic tests of varying accuracy, and was computed after the completion of the burn-in period. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-auc-mag-heatmap-dynamical-1x>

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_dynamical_0.102.svg"),
  caption: [The strength of the correlation ($|"AUC" - 0.5|$) for each EWS metric with emergence, at high levels of dynamical noise, for diagnostic tests of varying accuracy, and was computed after the completion of the burn-in period. The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-auc-mag-heatmap-dynamical-7x>

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_poisson_1.0x.svg"),
  caption: [The maximal alert accuracy achieved by each EWS metric under low levels of static noise. Q) refers to the long-running quantile threshold to return a flag, and C) the number of consecutive flags to trigger an alert, that in combination produce the maximal accuracy. S) refers to the resulting specificity of the alert system.  The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-accuracy-heatmap-poisson-1x>

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_poisson_7.0x.svg"),
  caption: [The maximal alert accuracy achieved by each EWS metric under high levels of static noise. Q) refers to the long-running quantile threshold to return a flag, and C) the number of consecutive flags to trigger an alert, that in combination produce the maximal accuracy. S) refers to the resulting specificity of the alert system.  The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-accuracy-heatmap-poisson-7x>

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_dynamical_0.8734.svg"),
  caption: [The maximal alert accuracy achieved by each EWS metric under low levels of dynamical noise. Q) refers to the long-running quantile threshold to return a flag, and C) the number of consecutive flags to trigger an alert, that in combination produce the maximal accuracy. S) refers to the resulting specificity of the alert system.  The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-accuracy-heatmap-dynamical-1x>

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_dynamical_0.102.svg"),
  caption: [The maximal alert accuracy achieved by each EWS metric under high levels of dynamical noise. Q) refers to the long-running quantile threshold to return a flag, and C) the number of consecutive flags to trigger an alert, that in combination produce the maximal accuracy. S) refers to the resulting specificity of the alert system.  The test sensitivity equals the test specificity for all diagnostic tests.]
)
<fig_csd-accuracy-heatmap-dynamical-7x>

== Tables

#let tau_comparison_table = csv("./manuscript_files/tables/tau-comparison.csv")
#let tau_comparison_vals = rename_noise_extract_vals(tau_comparison_table)

#figure(
  block[
    #set text(size: 11pt)
    #three_header_table(
      columns: 6,
      align: horizon,
      table.cell(rowspan: 3, align: horizon)[Rank], [Perfect Test], table.cell(colspan: 4)[90% Sensitive & Specific Imperfect Test],
      table.cell(rowspan: 2)[All Noise], table.cell(colspan: 2)[Static Noise], table.cell(colspan:2)[Dynamical Noise],
      ..tau_comparison_vals
    )
  ],
  caption: [The ranking and mean value of Kendall's Tau computed on the subset of the emergent time series after the burn-in period, for a perfect test and an imperfect test with sensitivity and specificity equal to 90%, under high and low static and dynamical noise systems]
)
<tbl-tau-ranking-rdt-comparison>

#let accuracy_comparison_table = csv("./manuscript_files/tables/accuracy-comparison.csv")
#let accuracy_comparison_vals = rename_noise_extract_vals(accuracy_comparison_table)

#figure(
  block[
    #set text(size: 11pt)
    #three_header_table(
      columns: 6,
      align: horizon,
      table.cell(rowspan: 3, align: horizon)[Rank], [Perfect Test], table.cell(colspan: 4)[90% Sensitive & Specific Imperfect Test],
      table.cell(rowspan: 2)[All Noise], table.cell(colspan: 2)[Static Noise], table.cell(colspan:2)[Dynamical Noise],
      ..accuracy_comparison_vals
    )
  ],
  caption: [The ranking and alert accuracy of the EWS-based alert system computed on the subset of the emergent time series after the burn-in period, for a perfect test and an imperfect test with sensitivity and specificity equal to 90%, under high and low static and dynamical noise systems]
)
<tbl-accuracy-ranking-rdt-comparison>

