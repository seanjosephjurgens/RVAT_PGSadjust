Some example scripts for rare variant analysis within the UK Biobank DNAnexus RAP, where we apply PGS adjustment to improve power in rare variant association discovery. Analyses reported in the manuscript: Jurgens, Pirruccello, et al. (2022) Adjusting for common variant polygenic scores improves power in rare variant association analyses. (under consideration). 

The scripts are made for use on the UKB DNAnexus RAP (with some adjustments). Please read the hashes in these files as they may refer to other files that need to be made to run the analyses (such as phenotype files).


####################################################################################


The _genotype_array_file_forGRM_ file contains scripts used for pruning and merging of the genotyping array data for use in the genetic relatedness matrices in the mixed effects models.


The _SAIGE_GENE+_ directory contains bash scripts used for submission of SAIGE-GENE+ analyses on the UKB DNAnexus RAP, utilizing speed-optimized sparse mixed-effects models.


The _BOLT-LMM_ directory contains bash scripts used for submission of BOLT-LMM analyses on the UKB DNAnexus RAP.


####################################################################################

NB:
The current scripts are optimized for analysis of rare and ultra-rare variation, e.i. MAF<0.1% and MAF<0.01%, in gene-based testing. The scripts in this repos are restricted to association analyses in SAIGE-GENE+ and BOLT-LMM, and do not include those for PGS construction (as investigators may use various methods to construct their PGS). Of note, we recommend trying out Leave-One-Chromosome-Out (LOCO) PGS when applying our method to low-frequency variants (e.g. MAF<1% or MAF<5%) if investigators want to avoid linkage between PGS-variants and the tested variants. For rare variation (MAF<0.1%), however, we find that LOCO scores perform similarly to all chromosome PGS (given very minimal linkage between rare and common variants).
