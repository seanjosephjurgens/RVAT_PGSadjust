#! /bin/bash

# Analyses run on Swiss Army Knife on the UKB RAP

PROJECT="RAP_PROJECT_NAME"
TRAIT_NUM=NUM # the integer trait number of 65 traits assessed in the current study
mkdir tmp/
cd tmp/

echo "READING IN INPUT FILES"
# Import genotyping array data previously pruned to ~280k array variants (see the 'genotype_array_file_forGRM' file from this repos)
dx download ${PROJECT}:/path/to/pruned_array_data/ukbb200k_array_chrall.*

# Import phenotype and covar file for QCd samples from the initial 200k freeze of UKB, which includes the out-of-sample derived PGS for all traits
dx download ${PROJECT}:/path/to/phenofile/phenofile_wes_total_quantitative_forBOLT.txt

# Input names
PLINK_FILE=ukbb200k_array_chrall
PHENO_FILENAME=phenofile_wes_total_quantitative_forBOLT.txt

# Specify WES plink dataset path
WES_PLINK_PREFIX=/mnt/project/path/to/plink/file_chromosome
WES_PLINK_SUFFIX=_genotype_variant_sample_QCd

# REGENIE annotation files; see REGENIE documentation for more info
ANNO_FILE=/mnt/project/path/to/REGENIE/groupingfile_chr${chr}.annotationfile.tsv
SETLIST_FILE=/mnt/project/path/to/REGENIE/groupingfile_chr${chr}.setlistfile.tsv
AAF_FILE=/mnt/project/path/to/REGENIE/groupingfile_chr${chr}.AAFfile.tsv
SETINCLUSION_FILE=/mnt/project/path/to/REGENIE/groupingfile_chr${chr}.setinclusionfile.tsv
MASK_FILE=/mnt/project/path/to/REGENIE/groupingfile_chr${chr}.maskdefinitionfile.tsv

# Set up REGENIE
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.2/regenie_v3.2.2.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.2.2.gz_x86_64_Linux_mkl.zip

# Traits analyzed in the current study
traits=(
invnorm_whitebloodcellleukocytecount
invnorm_redbloodcellerythrocytecount
invnorm_haemoglobinconcentration
invnorm_haematocritpercentage
invnorm_meancorpuscularvolume
invnorm_meancorpuscularhaemoglobin
invnorm_meancorpuscularhaemoglobinconcentration
invnorm_redbloodcellerythrocytedistributionwidth
invnorm_plateletcount
invnorm_plateletcrit
invnorm_meanplateletthrombocytevolume
invnorm_plateletdistributionwidth
invnorm_lymphocytecount
invnorm_monocytecount
invnorm_neutrophillcount
invnorm_eosinophillcount
invnorm_basophillcount
invnorm_nucleatedredbloodcellcount
invnorm_lymphocytepercentage
invnorm_monocytepercentage
invnorm_neutrophillpercentage
invnorm_eosinophillpercentage
invnorm_basophillpercentage
invnorm_nucleatedredbloodcellpercentage
invnorm_reticulocytepercentage
invnorm_reticulocytecount
invnorm_meanreticulocytevolume
invnorm_meanspheredcellvolume
invnorm_immaturereticulocytefraction
invnorm_highlightscatterreticulocytepercentage
invnorm_highlightscatterreticulocytecount
invnorm_alanineaminotransferase
invnorm_albumin
invnorm_alkalinephosphatase
invnorm_apolipoproteina
invnorm_apolipoproteinb
invnorm_aspartateaminotransferase
invnorm_creactiveprotein
invnorm_calcium
invnorm_cholesterol
invnorm_creatinine
invnorm_cystatinc
invnorm_directbilirubin
invnorm_gammaglutamyltransferase
invnorm_glucose
invnorm_glycatedhaemoglobinhba1c
invnorm_hdlcholesterol
invnorm_igf1
invnorm_ldldirect
invnorm_lipoproteina
invnorm_phosphate
invnorm_shbg
invnorm_testosterone
invnorm_totalbilirubin
invnorm_totalprotein
invnorm_triglycerides
invnorm_urate
invnorm_urea
invnorm_vitamind
invnorm_bmi
invnorm_sbp
invnorm_dbp
invnorm_pulse_rate
invnorm_weight_kg
invnorm_height_cm
)
#phenoList=$(printf ",%s" "${traits[@]}")
#phenoList=${phenoList:1}
TRAIT=${traits[TRAIT_NUM-1]}

echo "STARTING ANALYSIS"

#### Outputs
STEP1_OUTPUT=STEP1_out

echo "STARTING REGENIE STEP1"

# Run REGENIE STEP 1
./regenie_v3.2.2.gz_x86_64_Linux_mkl \
    --step 1 \
    --qt \
    --bed ${PLINK_FILE} \
    --covarFile ${PHENO_FILENAME} \
    --covarColList sequencing_batch,genotyping_array,Inferred_Gender,enroll_age,enroll_age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
    --catCovarList sequencing_batch,genotyping_array,Inferred_Gender \
    --phenoFile ${PHENO_FILENAME} \
    --phenoColList ${TRAIT} \
    --bsize 1000 \
    --out ${STEP1_OUTPUT}

echo "STARTING REGENIE STEP2"

for chr in {1..22};
do
echo "   BUSY WITH CHROMOSOME ${chr}"
# Run REGENIE STEP 2
./regenie_v3.1.3.gz_x86_64_Linux_mkl \
    --step 2 \
    --pgen ${WES_PLINK_PREFIX}${chr}${WES_PLINK_SUFFIX} \
    --covarFile ${PHENO_FILENAME} \
    --covarColList sequencing_batch,genotyping_array,Inferred_Gender,enroll_age,enroll_age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
    --catCovarList sequencing_batch,genotyping_array,Inferred_Gender \
    --phenoFile ${PHENO_FILENAME} \
    --phenoColList ${TRAIT} \
    --pred ${STEP1_OUTPUT}_pred.list \
    --anno-file ${ANNO_FILE} \
    --set-list ${SETLIST_FILE} \
    --aaf-file ${AAF_FILE} \
    --extract-sets ${SETINCLUSION_FILE} \
    --mask-def ${MASK_FILE} \
    --aaf-bins 0.0004 \
    --minMAC 20 \
    --print-pheno \
    --bsize 200 \
    --out ../REGENIE_STEP2_${TRAIT}_noPRS
done

echo "DONE WITH ANALYSIS."

cd ..
rm -rf tmp

echo "MOVING FILES TO PROJECT"
# When using Swiss Army knife, the files will be automatically uploaded to the destination directory
