#! /bin/bash

# Analyses run on Swiss Army Knife on the UKB RAP

PROJECT="RAP_PROJECT_NAME"
TRAIT_NUM="TRAIT_NUMBER" # the integer trait number of 65 traits assessed in the current study
mkdir tmp/
cd tmp/

echo "READING IN INPUT FILES"
# Import BOLT files that come with BOLT for running analyses
dx download ${PROJECT}:path/to/BOLT_files/genetic_map_hg19_withX.txt.gz
dx download ${PROJECT}:path/to/BOLT_files/LDSCORE.1000G_EUR.tab.gz

# Import genotyping array data previously pruned to ~280k array variants (see the 'genotype_array_file_forGRM' file from this repos)
dx download ${PROJECT}:path/to/pruned_array_data/ukbb200k_array_chrall.*

# Import the gene-burdens as saved in bgen format (so each gene collaps/burden is saved as thought being a variant).
dx download ${PROJECT}:path/to/gene_burdens/in_bgen_format/ukbb200k_wesQC_chrall_carrier.bgen
dx download ${PROJECT}:path/to/gene_burdens/in_bgen_format/ukbb200k_wesQC_chrall_carrier.sample

# Import phenotype and covar file for QCd samples from the initial 200k freeze of UKB, which includes the out-of-sample derived PGS for all traits
##### The column for the out-of-sample PGS based on PRScs-auto is named according to the format 'PRSprscsauto_${TRAIT}' in this file
dx download ${PROJECT}:path/to/phenofile/phenofile_wes_total_quantitative_forBOLT.txt

#### Inputs
PLINK_FILE=ukbb200k_array_chrall
PHENO_FILENAME=phenofile_wes_total_quantitative_forBOLT.txt
SAMPLE_FILE=ukbb200k_wesQC_chrall_carrier.sample
BGEN_FILE=ukbb200k_wesQC_chrall_carrier.bgen
LDSCORES_FILE=LDSCORE.1000G_EUR.tab.gz
GENETIC_MAP_FILE=genetic_map_hg19_withX.txt.gz

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
TRAIT=${traits[TRAIT_NUM-1]}
echo "STARTING ANALYSIS"
echo "RUNING TRAIT: "${TRAIT}

#### Outputs
IMPUTED_OUTPUT="../ukbb200k_PGSadjust_results_${TRAIT}_PRSprscsauto_bgenoutput"
GT_OUTPUT="../ukbb200k_PGSadjust_results_${TRAIT}_PRSprscsauto_statsoutput"

echo "STARTING BOLT"
bolt \
    --bfile=${PLINK_FILE} \
    --phenoFile=${PHENO_FILENAME} \
    --phenoCol=${TRAIT} \
    --covarFile=${PHENO_FILENAME} \
    --sampleFile=${SAMPLE_FILE} \
    --bgenFile=${BGEN_FILE} \
    --bgenMinMAF=0 \
    --bgenMinINFO=0 \
    --statsFileBgenSnps=${IMPUTED_OUTPUT} \
    --covarMaxLevels=1000 \
    --covarCol=sequencing_batch \
    --covarCol=genotyping_array \
    --covarCol=Inferred_Gender \
    --qCovarCol=PRSprscsauto_${TRAIT} \
    --qCovarCol=enroll_age \
    --qCovarCol=enroll_age2 \
    --qCovarCol=PC{1:20} \
    --LDscoresFile=${LDSCORES_FILE} \
    --geneticMapFile=${GENETIC_MAP_FILE} \
    --LDscoresMatchBp \
    --lmm \
    --numThreads=$(nproc) \
    --maxMissingPerSnp=0.01 \
    --maxMissingPerIndiv=0.25 \
    --statsFile=${GT_OUTPUT} \
    --lmmForceNonInf \
    --verboseStats

# Note that --maxMissingPerIndiv is now set to a very high number so that no
# additional filtering will be done beyone what you did to create pheno.tsv

echo "DONE WITH ANALYSIS."

cd ..
rm -rf tmp

echo "MOVING FILES TO PROJECT"
# When using Swiss Army Knife files will be automatically moved to the destination directory
