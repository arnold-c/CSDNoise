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
)

#set bibliography(style: "elsevier-vancouver")

#bibliography("CSD.bib")
