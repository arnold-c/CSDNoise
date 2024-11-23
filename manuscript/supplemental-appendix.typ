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
)

#show figure.where(kind: table): set figure(supplement: [Supplemental Table])
#show figure.where(kind: image): set figure(supplement: [Supplemental Figure])

= Tables

#let tau_comparison_table = csv("./manuscript_files/tables/tau-comparison.csv")
#figure(
  two_header_table(
    columns: 6,
    table.cell(rowspan: 2, align: horizon)[Rank], [Perfect Test], table.cell(colspan: 4)[90% Sensitive & 90% Specific RDT],
    ..tau_comparison_table.flatten().slice(1)
  ),
  caption: [The ranking and mean value of Kendall's Tau computed on the subset of the emergent time series after the burn-in period, for a perfect test and an RDT with 90% sensitivity and 90% specificity, under high and low Poisson and dynamical noise systems]
)
<tbl-tau-ranking-rdt-comparison>


#let auc_comparison_table = csv("./manuscript_files/tables/auc-comparison.csv")

#figure(
  two_header_table(
    columns: 6,
    table.cell(rowspan: 2, align: horizon)[Rank], [Perfect Test], table.cell(colspan: 4)[90% Sensitive & 90% Specific RDT],
    ..auc_comparison_table.flatten().slice(1)
  ),
  caption: [The ranking of AUC computed on the subset of the emergent time series after the burn-in period, for a perfect test and an RDT with 90% sensitivity and 90% specificity, under high and low Poisson and dynamical noise systems]
)
<tbl-auc-ranking-rdt-comparison>

#let alert_accuracy_auc_magnitude_comparison_table = csv("./manuscript_files/tables/alert-accuracy-auc-comparison.csv")

#figure(
    two_header_table(
    columns: 7,
    table.cell(rowspan: 2, align: horizon)[Rank], table.cell(colspan:2)[ Perfect Test], table.cell(colspan: 4)[90% Sensitive & 90% Specific RDT],
    ..alert_accuracy_auc_magnitude_comparison_table.flatten().slice(1)
  ),
  caption: [The ranking and $|"AUC" - 0.5|$ for each metric computed on the emergent time series with a perfect test, and the alert accuracy with an RDT. The values are computed on the full time series, and the subset from after the completion of the burn-in period, with a perfect test]
)
<tbl-alert-accuracy-auc-rdt>


= Plots

#figure(
  image("./supplemental_files/plots/accuracy-line-plot.svg"),
  caption: [The change in alert accuracy for the least correlated EWS metrics under increasing diagnostic uncertainty, and low and high levels of Poisson or dynamical noise. Low noise refers to simulations where the average incidence of noise is equal to the average incidence of measles. High noise refers to simulations where the average incidence of noise is equal to 7 times the average incidence of measles. The tests sensitivity equals the test specificity for all diagnostic tests.]
)

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-magnitude-heatmaps/after-burnin/tau_auc-magnitude-heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x]
)

== Optimal Threshold Accuracies

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_poisson_1.0x.svg"),
  caption: [The maximal alert accuracy under 1x Poisson noise. P) refers to the long-running percentile threshold to return a flag, and C) the number of consecutive flags to trigger and alert, that in combination produce the maximal accuracy. S) refers to the specificity of the alert system]
)

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_poisson_7.0x.svg"),
  caption: [The maximal alert accuracy under 7x Poisson noise. P) refers to the long-running percentile threshold to return a flag, and C) the number of consecutive flags to trigger and alert, that in combination produce the maximal accuracy. S) refers to the specificity of the alert system]
)

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_dynamical_0.8734.svg"),
  caption: [The maximal alert accuracy under 1x Dynamical noise. P) refers to the long-running percentile threshold to return a flag, and C) the number of consecutive flags to trigger and alert, that in combination produce the maximal accuracy. S) refers to the specificity of the alert system]
)

#figure(
  image("manuscript_files/plots/optimal-threshold-heatmaps/optimal_heatmap_dynamical_0.102.svg"),
  caption: [The maximal alert accuracy under 7x Dynamical noise. P) refers to the long-running percentile threshold to return a flag, and C) the number of consecutive flags to trigger and alert, that in combination produce the maximal accuracy. S) refers to the specificity of the alert system]
)


#set bibliography(style: "elsevier-vancouver")

#bibliography("CSD.bib")
