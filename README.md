# Mallard_Robertson
 
The data directory contains:
- the transition rates in High Salt: MA_lines_phenotypes.txt
- the body size measurements : body_area_means.txt
- fertility data from the Cemee: fertility_with_n.rda (from Noble et al. 2017)
- results from the competition experiment: competition.txt

The Rcode directory contains all code to reproduce the analysis.

We also included txt output files showing the different computed G/M matrices (txt, G_mat_tables, G_matrix_eigendecomposition directories) and the auto-correlation plots from the MCMCglmm models. Some of these tables have been formatted for better readability in the "formatted_tables" directory.

RData files produced by the code are also shared in the "RData" (most files) and "RData_Random"" (Randomized G matrices) directories.
Two of them, the Rdata file produced after during the Tensor analysis (see Appendix) and the one produced when comparing different priors for the computation of the ancestral G matrix are too large to be stored in github and will be stored on Dryad upon publication.
