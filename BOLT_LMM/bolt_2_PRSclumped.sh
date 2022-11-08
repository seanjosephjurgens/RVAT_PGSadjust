#! /bin/bash

TRAIT_NUM=$1
mkdir tmp/
cd tmp/

echo "READING IN INPUT FILES"
# Import BOLT files
dx download exome-seq:sjj/short_projects/PRSadjust/R2/data/BOLT_files/genetic_map_hg19_withX.txt.gz
dx download exome-seq:sjj/short_projects/PRSadjust/R2/data/BOLT_files/LDSCORE.1000G_EUR.tab.gz

# Import genotyping array; 89k version
dx download exome-seq:sjj/short_projects/PRSadjust/R2/data/step1_genetic/ukbb200k_array_chrall.*

# Import bgen format for gene-burdens
dx download exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/ukbb200k_wesQC_chrall_carrier.bgen
dx download exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/ukbb200k_wesQC_chrall_carrier.sample

# Import pheno/covar file
dx download exome-seq:sjj/short_projects/PRSadjust/R2/data/pheno/phenofile_wes_total_quantitative_forBOLT.txt

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
IMPUTED_OUTPUT="../ukbb200k_PGSadjust_results_${TRAIT}_PRSclumped_bgenoutput"
GT_OUTPUT="../ukbb200k_PGSadjust_results_${TRAIT}_PRSclumped_statsoutput"

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
    --qCovarCol=PRSclumped_${TRAIT} \
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
