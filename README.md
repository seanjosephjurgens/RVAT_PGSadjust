Some example scripts for rare variant analysis within the UK Biobank DNAnexus RAP, where we apply PGS adjustment to improve power in rare variant association discovery. Analyses reported in the manuscript: Jurgens, Pirruccello, et al. (2022) Adjusting for common variant polygenic scores improves power in rare variant association analyses. (under consideration). 

The scripts are made for use on the UKB DNAnexus RAP (with some adjustments). Please read the hashes in these files as they may refer to other files that need to be made to run the analyses (such as phenotype files).

################################################################################################################################################

The _genotype_array_file_forGRM_ file contains scripts used for pruning and merging of the genotyping array data for use in the genetic relatedness matrices in the mixed effects models.

The _BOLT-LMM_ directory contains bash scripts used for submission of BOLT-LMM analyses on the UKB DNAnexus RAP.
