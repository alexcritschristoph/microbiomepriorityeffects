# microbiomepriorityeffects

## Reena Debray, Robin Herbert, Alexander Jaffe, Alexander Crits-Christoph, Mary Power, Britt Koskella


## Abstract
Advances in next-generation sequencing have enabled the widespread measurement of microbiome composition across systems and over the course of microbiome assembly. Despite substantial progress in understanding deterministic drivers of community composition, the role of historical contingency has received less attention. The ability of a given strain to establish in a community can depend on the abundances of other strains at or before its arrival, a phenomenon known as a priority effect. Here, we review mechanisms of priority effects and evidence for their importance in microbial communities. We describe approaches for directly testing and predicting priority effects in complex microbial communities and illustrate with re-analysis of publicly available plant and animal microbiome datasets. Observed priority effects can recapitulate known interactions and identify novel candidates for experimental validation. Finally, we discuss shared principles that emerge across study systems, focusing on eco-evolutionary dynamics and the importance of scale. Overall, we argue that predicting when and how current community state impacts success of newly arriving microbial taxa is critical for understanding and manipulating the community assembly of microbiomes. We conclude by discussing outstanding conceptual and practical challenges faced when measuring priority effects in microbiomes.

## Funding
R.D. was supported by a graduate fellowship from the National Science Foundation.

# This repository includes:
## Arrival and persistence analysis.R
This script contains all functions needed for the analysis, and shows how they were implemented on temporal amplicon sequence data from human, mouse, and cattle intestinal microbiomes to generate Figure 3-4 in the manuscript.
## OTU tables and taxonomy
The files "humangut_otus.csv", "humangut_taxonomy.csv", "murine_otus.csv", "murine_taxonomy.csv", "rumen_otus.csv", "rumen_taxonomy.csv" were obtained from previously published work as described, and are provided here as used for the re-analysis.
## Summary tables
The files "summary_humangut_OTUs", "summary_murine_OTUs", and "summary_rumen_OTUs" contain summary information on mean arrival time, mean persistence, population prevalence, taxonomy, and dependence of persistence on arrival for each OTU that passed our initial filtering criteria. The files "humangut_prior_deseq",  "murine_prior_deseq", and "rumen_prior_deseq" contain strains identified by the negative binomial modeling package "DESeq2" to positively or negatively predict the persistence of a focal strain after arrival.
## Permutation analysis
The files "humangut_perm_prior", "murine_perm_prior", and "rumen_perm_prior" each contain a list of 1000 values, where each value represents the taxonomic overlap (proportion of pairs in the same family) of randomly sampled OTU pairs. The number of randomly sampled OTU pairs in each iteration matches the number of observed pairs identified by DESeq2 in each dataset. These permuted values were used to generate the null distributions depicted in Figure 4 as grey histograms.
