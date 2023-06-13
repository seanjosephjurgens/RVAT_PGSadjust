# RVAT PGS adjustment

## Publication
Some example scripts for rare variant analysis within the UK Biobank DNAnexus RAP, where we apply PGS adjustment to improve power in rare variant association discovery. Analyses reported in the manuscript ```Jurgens, Pirruccello, et al. (2023). Adjusting for common variant polygenic scores improves yield in rare variant association analyses. Nature Genetics.``` https://www.nature.com/articles/s41588-023-01342-w.

## Scripts
The scripts are made for use on the UKB DNAnexus RAP -- Cloud-based platform (with some adjustments for file paths and names). Please read the hashes in these files as they may refer to other files that need to be made to run the analyses (such as phenotype files).

* The ```PRScs``` directory contains general bash scripts to run the PRScs algorithm for construction of genome-wide polygenic scores, as applied in our study.
* The ```genotype_array_file_forGRM``` file contains scripts used for pruning and merging of the genotyping array data for use in the genetic relatedness matrices in the mixed effects models.
* The ```SAIGE_GENE+``` directory contains bash scripts used for submission of SAIGE-GENE+ analyses on the UKB DNAnexus RAP, utilizing speed-optimized sparse mixed-effects models.
* The ```BOLT-LMM``` directory contains bash scripts used for submission of BOLT-LMM analyses on the UKB DNAnexus RAP.

## Summary statistcs
The summary statistics, including the PRS weights used in our study (and summary stats for all main GWAS and RVAS), have been submitted to the Cardiovascular Disease Knowledge Portal (https://cvd.hugeamp.org/downloads.html). We recommend trying out the out-of-sample PRS for some traits using SAIGE-GENE+; good traits to test include Height, LDL cholesterol and HDL cholesterol.

## Final remarks
NB 1: The scripts in this repos are mainly focussed on SAIGE-GENE+ and BOLT-LMM (as these are nicely implementable in the UKB RAP). We have also included a general script for running PRScs that was used in our study, although investigators may use various methods to construct their PGS. Of note, we have submitted all PGS weights used in our study to the Cardiovascular Disease Knowledge Portal for download (https://cvd.hugeamp.org/downloads.html) so people can replicate analyses with our scores! 

NB 2: The current scripts are optimized for analysis of rare and ultra-rare variation, i.e. MAF<0.1% and MAF<0.01%, in gene-based testing. We recommend trying out Leave-One-Chromosome-Out (LOCO) PGS when applying our method to low-frequency variants (e.g. MAF<1% or MAF<5%) if investigators want to avoid linkage disequillibrium (LD) between PGS-variants and the tested variants. For rare variation (MAF<0.1%), however, we find that LOCO scores perform similarly to all chromosome PGS (given very minimal LD between such rare variants and common variants).

## Citation
If you use our approach our code, please cite:
```
Jurgens, S.J., Pirruccello, J.P., Choi, S.H. et al. Adjusting for common variant polygenic scores improves yield in rare variant association analyses. Nat Genet 55, 544–548 (2023). https://doi.org/10.1038/s41588-023-01342-w
```
