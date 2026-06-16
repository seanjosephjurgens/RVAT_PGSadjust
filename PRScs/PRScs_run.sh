#!/bin/bash

# This script expects that PRScs (https://github.com/getian107/PRScs) has been
# installed at /install/PRScs/PRScs.py and it requires the following variables
# to be set in the environment:

#PRSCS_LD_REF: The path to the ldblk_* file precomputed by the PRScs authors (or
#a custom LD reference panel). See
#https://github.com/getian107/PRScs#getting-started for precomputed files that
#can be downloaded for this purpose. The `ldblk_ukbb_eur.tar.gz` data was used
#for the present manuscript.

#SUMSTATS: The path to a BOLT-LMM format GWAS summary statistics file that will
#be used to create the polygenic score.

#PRSCS_OUT: The output path where PRScs output will be placed.

#chrom: The chromosome on which PRScs is being computed. Run separately per
#chromosome.

#APPROXIMATE_N_GWAS: The sample size of the GWAS.

#PRSCS_PHI: The Phi parameter from PRScs. If not set, PRScs will be run in
#'auto' mode, which is the mode used for this manuscript.

# # Script runs below here

# Use the same seed every time so that the program creates reproducible output
NOTHING_UP_MY_SLEEVE=31337

# # INPUT VARS
PRSCS_LD_REF_DIR="$(dirname "${PRSCS_LD_REF}")"

# # OUTPUT VARS
PRSCS_OUT_DIR="$(dirname "${PRSCS_OUT}")"

# PHI
PRSCS_PHI_LINE="--phi=${PRSCS_PHI}"
if [[ $PRSCS_PHI = 'auto' ]] || [[ $PRSCS_PHI = '' ]]; then
  PRSCS_PHI_LINE=""
fi

# Produce a PRScs bim file from a BOLT-LMM formatted summary statistics file.
echo "Producing BIM"
cat ${SUMSTATS} | tail -n +2 | awk 'BEGIN {FS="\t"; OFS=" "} {print $2 OFS $1 OFS "0" OFS $3 OFS $5 OFS $6}' > converted.bim
OPTIONAL_BIM="${PWD}/converted.bim"
# Strip out the final .bim to get to the path+name format expected by prscs
BIM_PREFIX=${OPTIONAL_BIM%.*}

# Produce a PRScs input-summary-stats-formatted file from a BOLT-LMM formatted
# summary statistics file.
echo "Converting sumstats to PRScs format"
echo "SNP	A1	A2	BETA	P" > converted.sumstats
cat ${SUMSTATS} | tail -n +2 | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1 OFS $5 OFS $6 OFS $11 OFS $16}' >> converted.sumstats

CONVERTED_SUMSTATS="${PWD}/converted.sumstats"

echo "Running PRScs with --phi=${PRSCS_PHI}"
/install/PRScs/PRScs.py \
    --ref_dir=${PRSCS_LD_REF_DIR} \
    --bim_prefix=${BIM_PREFIX} \
    --sst_file=${CONVERTED_SUMSTATS} \
    --n_gwas=${APPROXIMATE_N_GWAS} \
    --chrom=${chrom} \
    --seed=${NOTHING_UP_MY_SLEEVE} \
    ${PRSCS_PHI_LINE} \
    --out_dir=${PRSCS_OUT_DIR}/chr${chrom}
