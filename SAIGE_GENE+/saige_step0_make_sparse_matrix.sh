#! /bin/bash

# Run on Swiss Army Knife on UKB RAP

mkdir tmp/
cd tmp/

### This script requires loading of the SAIGE docker by specifiying in dx run input!
### Input files for bed/bim/fam need to be specified in dx run input!

LD_pruned_PLINK_file="ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01"
OUTNAME=../saige_sparse_matrix

ls -la
#Try
docker run wzhou88/saige:1.0.9 createSparseGRM.R \
    --plinkFile=${LD_pruned_PLINK_file} \
    --nThreads=94  \
    --outputPrefix=${OUTNAME}       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05

#Try
Rscript createSparseGRM.R   \
    --plinkFile=${LD_pruned_PLINK_file} \
    --nThreads=94  \
    --outputPrefix=${OUTNAME}       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05

cd ..
rm -rf tmp/
