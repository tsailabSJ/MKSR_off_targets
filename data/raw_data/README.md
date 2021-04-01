## Please download the fastq files from NCBI SRA!

We have 3 sets of fastq files, corresponding to cas-OFFinder rhAmp-seq panel, circle-seq 1 rhAmp-seq panel, and circle-seq 2 rhAmp-seq panel (added during revision). Off-targets were examined in 2 donors with 2 replicates. Fastq file names and treatment/control info are listed as below:

| Samples           | rhAmpSeq CasOFFinder sites | rhAmpSeq CIRCLE-seq sites panel 1 | rhAmpSeq CIRCLE-seq sites panel 2 |
|-------------------|----------------------------|-----------------------------------|-----------------------------------|
| #2 mRNA_MC008     | CRL1962                    | CRL1974                           | CRL2350                           |
| #2 protein_MC008  | CRL1963                    | CRL1975                           | CRL2351                           |
| #2 UT_MC008       | CRL1964                    | CRL1976                           | CRL2352                           |
| #3 mRNA_MC008     | CRL1965                    | CRL1977                           | CRL2353                           |
| #3 protein_MC008  | CRL1966                    | CRL1978                           | CRL2354                           |
| #3 UT_MC008       | CRL1967                    | CRL1979                           | CRL2355                           |
| #2 UT_MC0010      | CRL1968                    | CRL1980                           | CRL2356                           |
| #2 mRNA_MC0010    | CRL1969                    | CRL1981                           | CRL2357                           |
| #2 protein_MC0010 | CRL1970                    | CRL1982                           | CRL2358                           |
| #3 UT_MC0010      | CRL1971                    | CRL1983                           | CRL2359                           |
| #3 mRNA_MC0010    | CRL1972                    | CRL1984                           | CRL2360                           |
| #3 protein_MC0010 | CRL1973                    | CRL1985                           | CRL2361                           |

## Summary

This folder provides the code to generate the `*.allele.edit.tsv` in the processed data folder.

Cas-OFFinder panel:

- Assay_MKSR.bed
- MKSR_CasOFF_102820.bed
- crispressoPooled_BE_clazzaro_2020-12-03_info.tsv
- crispressoPooled_BE_clazzaro_2020-12-03.input

CIRCLE-seq 1 panel:

- Assay_cirMKSR.bed
- cirMKSR.gRNA.bed
- crispressoPooled_BE_clazzaro_2020-11-23_info.tsv
- crispressoPooled_BE_clazzaro_2020-11-23.input

CIRCLE-seq 2 panel:

- Assay_MakCIRC2.bed.rename.bed
- MKSR_cir2_final_clean.bed.rename.bed
- crispressoPooled_BE_yli11_2021-03-29_info.tsv
- crispressoPooled_BE_yli11_2021-03-29.input


`MKSR_CasOFF_102820.bed` provides the genomic coordiantes for the 129 on- and off-target sites in cas-OFFinder rhAmp-seq panel.`Assay_MKSR.bed` provides the coordiantes for amplicon sequences. 

`cirMKSR.gRNA.bed` provides the genomic coordiantes for the 90 on- and off-target sites in CIRCLE-seq rhAmp-seq panel. `Assay_cirMKSR.bed`  provides the coordiantes for amplicon sequences.

`MKSR_cir2_final_clean.bed.rename.bed` provides the genomic coordiantes for the 480 off-target sites in CIRCLE-seq panel. This is our second circle-seq rhAmp-seq panel. `Assay_MakCIRC2.bed.rename.bed` provides the coordiantes for amplicon sequences. 

`*_info.tsv` follows the format of `AMPLICONS_FILE.txt` from https://github.com/pinellolab/CRISPResso2

`*.input` is a tsv file providing the R1, R2 fastq files, sample label, and amplicon_info. Used for `BaseE.lsf`, for the LSF job array specification. If you don't have the LSF job system, you can run the example following the usage below. 

## Usage

This usage provides steps for you to run each fastq files and calculate editing efficiency.

Run the crispresso2 for all fastq files, here `$COL1` is R1 read, `$COL2` is R2 read, `$COL3` is file label, `$COL4` is a tsv file containing amplicon sequence and gRNA sequence.

```
CRISPRessoPooled -r1 ${COL1} -r2 ${COL2} --name ${COL3} -o results/${COL3}_results -f ${COL4} --quantification_window_size 10 --quantification_window_center -10 --base_editor_output --conversion_nuc_from A --conversion_nuc_to G

crispressoPooled_allele_edit_eff.py ${COL3} ${COL4} A G 

```

The output is `${COL3}.allele.edit.tsv`

## Usage on stjude HPC

https://hemtools.readthedocs.io/en/latest/content/CRISPR/crispressoPooled_BE.html

`crispressoPooled_BE.py` will generate a job script, `BaseE.lsf`, which is then submitted to the LSF job system using `bsub < BaseE.lsf`. Note that `crispressoPooled_BE.py` is specifically design to fit the stjude HPC, may not work outside stjude. 

```

crispressoPooled_BE.py -a Assay_MKSR.bed -gRNA MKSR_CasOFF_102820.bed -f fastq.tsv -g hg38 --ref A --alt G

crispressoPooled_BE.py -a Assay_MakCIRC2.bed.rename.bed -gRNA MKSR_cir2_final_clean.bed.rename.bed -f fastq.tsv -g hg38 --ref A --alt G

crispressoPooled_BE.py -a Assay_cirMKSR.bed -gRNA cirMKSR.gRNA.bed -f fastq.tsv -g hg38 --ref A --alt G

```
