## Summary

Here, we provide the code `chi_square.py` for identification of off-targets based on chi-square test.

Input: `design.tsv` and all tsv files

Output: `chi_square_test.csv`

Usage: `python chi_square.py`


`chi_square.py`: To calculate statistical significance of off-target editing for the ABE8e-NRCH mRNA or RNP treatments compared to control samples, we applied a Chi-square test for each of four samples (i.e., two donors, each with two replicates). The 2x2 contingency table was constructed based on the number of edited reads and the number of unedited reads in treated and control groups. FDR was calculated using the Benjamini/Hochberg method. The 54 reported significant off-targets were called based on: (1) FDR < 0.05 and (2) difference in editing frequency between treated and control > 0.5% for at least one treatment. 0.5% cutoff was determined by looking at the distribution of max_diff. (`*png`)


