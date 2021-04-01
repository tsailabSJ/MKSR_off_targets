# Quantification of base editing frequency

NOTE: Manuscript is under review. The use of any data in this repository is prohibited. 

This repository provides steps for: fastq -> processed data -> significant off-targets (Chi-square test)

## Data

See the `data` folder for raw data (fastq files have to be downloaded from NCBI SRA) and processed data.

The `raw_data` sub-folder contains inputs and codes to generate processed data

The `processed_data` sub-folder contains inputs and codes to generate the Chi-sqare test results

## Fastq -> Processed data

Raw fastq files were processed by CRISPRessoPooled (v2.0.41) with quantification_window_size 10, quantification_window_center -10, base_editor_output, conversion_nuc_from A, conversion_nuc_to G. 

`crispressoPooled_allele_edit_eff.py` reads the crispresso output and quantify editing efficiency at the "allele level". Specifically, the editing frequency for each site was calculated as the ratio between the number of reads containing the edited base (i.e., G) in a window from the 4th position to the 10th position of each protospacer and the total number of reads. Note that this window is hardcoded in the code, you can change it to other window size.

## Processed data to off-targets calling

`chi_square.py`: To calculate statistical significance of off-target editing for the ABE8e-NRCH mRNA or RNP treatments compared to control samples, we applied a Chi-square test for each of four samples (i.e., two donors, each with two replicates). The 2x2 contingency table was constructed based on the number of edited reads and the number of unedited reads in treated and control groups. FDR was calculated using the Benjamini/Hochberg method. The 54 reported significant off-targets were called based on: (1) FDR < 0.05 and (2) difference in editing frequency between treated and control > 0.5% for at least one treatment. 0.5% cutoff was determined by looking at the distribution of max_diff. 
