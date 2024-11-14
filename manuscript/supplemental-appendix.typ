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

= Tables

#let tau_comparison_table = csv("./manuscript_files/tables/tau-comparison.csv")
#figure(
  two_header_table(
    columns: 6,
    table.cell(rowspan: 2, align: horizon)[Rank], [Perfect Test], table.cell(colspan: 4)[90% Sensitive & 90% Specific RDT],
    ..tau_comparison_table.flatten().slice(1)
  ),
  caption: [The ranking and mean value of Kendall's #sym.tau computed on the subset of the emergent time series after the burn-in period, for a perfect test and an RDT with 90% sensitivity and 90% specificity, under high and low Poisson and dynamical noise systems]
)
<tbl-tau-ranking-rdt-comparison>

= Plots
== AUC Heatmaps
=== Full Length

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/full-length/tau_auc-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/full-length/tau_auc-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/full-length/tau_auc-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/full-length/tau_auc-heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x]
)


=== After 5yr Burn in

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/after-burnin/tau_auc-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/after-burnin/tau_auc-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/after-burnin/tau_auc-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x]
)

#figure(
  image("./manuscript_files/plots/tau_auc-heatmaps/after-burnin/tau_auc-heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x]
)

== Tau Heatmaps
=== Full Length

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/full-length/emergent-tau-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/full-length/emergent-tau-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x noise]
)

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/full-length/emergent-tau-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/full-length/emergent-tau-heatmap_dynamical_0.102.svg"),
  caption: [Dynamical noise, 7x noise]
)

=== After 5yr Burn in

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/after-burnin/emergent-tau-heatmap_poisson_1.0x.svg"),
  caption: [Poisson noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/after-burnin/emergent-tau-heatmap_poisson_7.0x.svg"),
  caption: [Poisson noise, 7x noise]
)

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/after-burnin/emergent-tau-heatmap_dynamical_0.8734.svg"),
  caption: [Dynamical noise, 1x noise]
)

#figure(
  image("manuscript_files/plots/tau_heatmaps/emergent/after-burnin/emergent-tau-heatmap_dynamical_0.102.svg"),
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


#set bibliography(style: "elsevier-vancouver")

#bibliography("CSD.bib")
