= Referee: 1

== Comments to the Author(s)
I believe this study is scientifically relevant, for it belongs to a large body of research that aims to understand how uncertainty limits the ability of EWS metrics to anticipate critical transitions.
In this particular case, uncertainty comes from co-circulating pathogens and imperfect diagnostic tests.
The authors show that, although robust to decreasing diagnostic accuracy, outbreak detection using EWS metrics is ineffective under high dynamical noise conditions, drastically reducing their utility.
Overall, I like the way the paper is written, as well as the quality and depth of the analyses presented.
I have only few minor comments and recommend it for publication once these and other reviewer comments are addressed.

== Minor Concerns
1. More details are needed about the use of Kaplan-Meier survival estimates to evaluate the speed and timing of outbreak detection.
2. I would like to see a clearer picture of the influence of noise magnitude on alert accuracy. In other words, how is the progression from the bottom left panel to the bottom right panel in Fig. 2, for different noise magnitude values?
3. I would like the authors to expand the discussion to include indications on how this framework could be validated using historical surveillance data. What are the challenges? What are the possible ways to characterize the "true" noise structure from real-world data?

= Referee: 2

Thank you for inviting me to review “Diagnostic Uncertainty Limits the Potential of Early Warning Signals to Identify Epidemic Emergence”.
As stated in the manuscript, this work aims to “characterize the performance of EWS metrics for outbreak detection in a surveillance system with diagnostic uncertainty due to co-circulating pathogens and imperfect diagnostic tests” , a goal that, if achieved, would represent an important contribution to the field of epidemiology.
However, I’m afraid this work requires substantial improvement in order to do so. The reasons for this assessment are as follows:

== Major Concerns
1. It is difficult to follow the overall structure of the manuscript. The authors make limited use of figures to support the arguments developed in the paper. I encourage the authors to make better use of visual elements. In particular, I suggest including a schematic figure that provides an overview of the approach taken.
2. I am concerned about the generalisability of the results. If the underlying structure of the compartmental model or the parameter values were changed, would the conclusions still hold?
3. The specific findings of the study are difficult to comprehend. Overall, the paper requires significant restructuring to clearly communicate the nature and contribution of the research. Although several EWS metrics are used, the paper does not explain why some perform well and others do not under the various scenarios considered.
4. The discussion section needs improvement, as it does not follow a logical structure in its current form. For example, the second paragraph reads more like content that belongs in the introduction.
5. It is difficult to assess the importance of the relationship between the hyperparameterisation step and the performance of a specific EWS metric. More intuitive explanations or justifications would help clarify this relationship.

== Minor Concerns
1. Why was a “burn-in” period of five years used? What is its biological meaning in this context? The term "burn-in" is typically used in statistical inference (e.g. MCMC sampling) to describe the period in which the sampler converges to the target distribution. I suggest using a different term to avoid adding to the existing ambiguity in infectious disease modelling terminology.
2. Lines 94–96. The sentence beginning “These bounds…” ends abruptly and unclearly with “was 0.9.” Please check.

= Referee: 3

== Comments to the Author(s)
In this manuscript, the authors analyse the impacts of diagnostic uncertainty in the presence of different underlying case noise regimes on the effectiveness of early warning signals.
They evaluate the performance of common early warning signals on time series of febrile rash generated from an SEIR model for measles subject to noise from a static Poisson distribution and from the combination of SEIR-rubella-like infection and Poisson-distributed non-rubella sources of febrile rash.
Previous studies have investigated the impacts of imperfect case reporting on observed early warning signals using a negative binomial case reporting distribution, but haven’t combined mis-reporting with an underlying noise time series.
It is therefore a valuable contribution to the literature that helps characterise some of the observed issues with applying early warning signals to real-world time series of observed disease cases.

Although the paper is well-presented and the analysis is novel, I have serious concerns about the methodology used by the authors as well as some minor concerns about the lack of clarity in the manuscript.
My recommendation is therefore to not accept the paper in its current form and instead require major revisions.
Below I list comments and suggestions to be taken into account before resubmission.

== Major Concerns

1. A major concern is the calculation of the sensitivity and specificity used. In the accompanying codebase (e.g. on lines 535-537 of ews-hyperparam-optimization.jl) , the authors appear to have incorrectly used the ratio of true positives/negatives to the total number of simulations. In line number 148 of the manuscript, the test accuracy is then defined as the average of these two values, which is an incorrect definition (unless the authors mean the balanced accuracy), but happens to be mathematically correct due to the incorrect calculations of sensitivity and specificity. This invalidates much of the authors’ results and discussion around optimal hyperparameter choices and the alert accuracy of early warning signals under different diagnostic uncertainty and noise regimes.
2. Another major concern regards the reproducibility and applicability of the simulation methodology used. Although the simulations may take time to run due to the large population size, using only 100 time series for each noise scenario decreases reproducibility (due to large stochastic uncertainty). Looking at Figure 3 in particular, the performance of the early warning signals are unintuitively not monotonic with respect to a decreasing test specificity/sensitivity. Related to this, tests generally make a trade-off between sensitivity and specificity so it seems odd to set these to be equal (lines 84-87). Could the authors justify this choice or extend their analysis to grid-search over a combination of sensitivity + specificity combinations? Further, could the authors include graphics of the noise time series for clarity.
3. Related to the above, a more recent paper on medRxiv (https://www.medrxiv.org/content/10.1101/2025.07.11.25331369v1) exists with the same simulation methodology and noise methods (though without mutual citations). It also includes the plots of the noise time series and explanation of some of the methodology. Could the authors comment on the relationship between these papers?
4. The next major concern is the lack of sensitivity analysis for some of the parameter choices in the analysis. Early warning signals are often susceptible to changes in rolling window size, have the authors conducted a sensitivity analysis to justify the bandwidth choice of 52 weeks in line 112? Further, the authors employ a 5 year burn-in period to ensure variability in Re before calculating early warning signals. Is there a reason why random initial conditions within a certain range were not used? Is there reason to think that wouldn’t have the same impact? Finally, is there a justification for calculating early warning signals on both the full time series and just after the burn-in period? Intuitively, the differences between these results are likely captured by comparing the |AUC-0.5| for null and emergent simulations.

== Minor Concerns

1. A minor concern, throughout the paper, is that methodological choices are reported without justification. In line number 73, the authors could justify why using 1 and 7 times the daily average measles incidence is appropriate (as the range seems quite large and also seems to be quite a high noise level relative to case counts). Similarly, the parameters used for a rubella-like disease should have a more clear literature citation. In line 75 the 15% value of non-rubella febrile rash could also be justified.
2. The authors could consider reformatting for clarity. In lines 3-6, the authors change from using M and B presumably as shorthand for millions and billions to using the full word. The full words could be used for clarity. In lines 73-74 two different versions of c are used to stand for the multiple of daily average measles incidence used. I would also recommend using the Greek letter when reporting on Kendall’s Tau (as is commonly done in the literature) so that in line 158 for example it becomes more clear what Tau is referring to.
3. There are certain mathematical points that could be better phrased or further explored for improved readability. In line 122 it might be useful to explain why Kendall’s Tau-B was used (as opposed to the usual Tau-A). Explicitly stating what you mean by ‘raw’ correlation (such as in lines 158 and 167) would also help clarity. Could the authors also be explicit if normalising the correlation means calculating the |AUC-0.5| statistic? Lines 113 and 129 refer to equations that are not numbered in the submission file, the equations should also be included as part of the sentence and not referenced like with figures. I would also remind readers in the results section that |AUC-0.5| close to 0.5 is a significant result as it is a less commonly used interpretation of AUC. With reference to Figure 3, although the results are interesting, survival curves are non-standard in early warning signal analysis and so should be explained in the manuscript.


= Referee: 4

== Comments to the Author(s)
This paper shows that statistical earning warning signals of critical slowing down can predict upcoming measles outbreaks with good accuracy, despite the presence of diagnostic uncertainty and confounding effects from co-circulating rubella, as long as the dynamical noise associated with the co-circulating disease is not too large.  The paper is well-written and the topic is interesting.  As the authors point out, EWS for infectious diseases is less advanced, and the authors have addressed system-specific problems that must be faced if such methods are to be deployed to provide early warning of potential outbreaks caused by dropping vaccine coverage. As such, it moves the field forward in concrete and useful ways.  I think the paper is publishable, aside from a few things that might help improve the writing and further consolidate the findings:

== Major revision
1. I did not see any sensitivity or uncertainty analysis with respect to epidemiological parameters for the R0, latent and infectious period, although the R0 can vary from place to place and the latent/infectious periods are always just an approximation since the true periods are distributed non-exponentially. The authors should conduct some sensitivity/uncertainty analysis to assess whether the accuracy remains good under some plausible variation in these parameter values.

== Minor comments

1. Heterogeneity could significantly impact the results. For instance, diagnostic uncertainty might be highest in exactly the parts of the population where vaccine coverage is lowest, and this might make it harder to find EWS by monitoring case notifications. The discussion section could bring this up as a limitation, with reference to existing literature on spatial/network models of heterogeneity in vaccination for infectious diseases.
2. It’s possible that a data-driven dynamical systems approach might be less confounded by dynamical noise. For instance, deep learning on time series features using a training library with simulated co-circulating rubella might enable the algorithm to distinguish features associated with rubella, and discard those to generate a better signal for measles.  The authors could mention this as an opportunity for future work.

