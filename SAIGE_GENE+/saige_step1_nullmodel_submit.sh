## Submission script to submit SAIGE step1 nullmodel fitting using Swiss Army knife

PROJECT="RAP_PROJECT_NAME"
TRAIT_NUM=TRAIT_NUM # integer trait number you want to analyze from the 65 quantitative traits assessed in our study

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
  height_cm
)

TRAIT=${traits[TRAIT_NUM-1]}

#####################
## No PGS added
#####################

### Submit noPRS nullmodel
instance_type="mem3_ssd1_v2_x4"
priority=normal
jobname=saige_step1_${TRAIT}_noPRS 
nThreads=3
destination=${PROJECT}:path/to/SAIGE_nullmodels/
traitType=quantitative
isCateVarianceRatio=TRUE
invNormalize=FALSE
useSparseGRMtoFitNULL=TRUE
relatednessCutoff=0.05
covariatesList=sex,enroll_age,enroll_age2,sequencing_batch,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20
qCovarColList=sex,sequencing_batch,genotyping_array
sampleIDCol=IID
PLINK_for_vr=ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01
PLINK_for_vr_path=${PROJECT}:path/to/array_data_for_saige/
sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
sparseGRM_path=${PROJECT}:path/to/SAIGE_step0_output_files/
phenoFile=phenofile_wes_total_quantitative_forBOLT.txt
phenoFile_path=${PROJECT}:path/to/phenofile/
phenoCol=${TRAIT}
outputPrefix=saige_step1_nullmodel_${TRAIT}_noPRS
tauInit=0,0 

dx run swiss-army-knife \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bed" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bim" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.fam" \
-iin="${phenoFile_path}${phenoFile}" \
-iin="${sparseGRM_path}${sparseGRMFile}" \
-iin="${sparseGRM_path}${sparseGRMSampleIDFile}" \
-icmd="step1_fitNULLGLMM.R --plinkFile=${PLINK_for_vr} --nThreads=${nThreads} --outputPrefix=${outputPrefix} --useSparseGRMtoFitNULL=${useSparseGRMtoFitNULL} --relatednessCutoff=${relatednessCutoff} --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --isCateVarianceRatio=${isCateVarianceRatio} --phenoFile=${phenoFile} --phenoCol=${phenoCol} --covarColList=${covariatesList} --qCovarColList=${qCovarColList} --sampleIDColinphenoFile=${sampleIDCol} --traitType=${traitType} --invNormalize=${invNormalize} --tauInit=${tauInit}" \
-iimage_file="${PROJECT}:path/to/docker/saige_1.0.9.tar.gz" \
--name ${jobname} \
--instance-type ${instance_type} \
--priority ${priority}  \
--yes \
--destination ${destination}


#####################
## Out-of-sample scores
#####################

### PRSclumped
instance_type="mem3_ssd1_v2_x4"
priority=normal
jobname=saige_step1_${TRAIT}_PRSclumped
nThreads=3
destination=exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/
traitType=quantitative
isCateVarianceRatio=TRUE
invNormalize=FALSE
useSparseGRMtoFitNULL=TRUE
relatednessCutoff=0.05
phenoCol=${TRAIT}
covariatesList=PRSclumped_${TRAIT},sex,enroll_age,enroll_age2,sequencing_batch,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20
qCovarColList=sex,sequencing_batch,genotyping_array
sampleIDCol=IID
PLINK_for_vr=ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01
PLINK_for_vr_path=${PROJECT}:path/to/array_data_for_saige/
sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
sparseGRM_path=${PROJECT}:path/to/SAIGE_step0_output_files/
phenoFile=phenofile_wes_total_quantitative_forBOLT.txt
phenoFile_path=${PROJECT}:path/to/phenofile/
phenoCol=${TRAIT}
outputPrefix=saige_step1_nullmodel_${TRAIT}_PRSclumped
tauInit=0,0

dx run swiss-army-knife \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bed" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bim" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.fam" \
-iin="${phenoFile_path}${phenoFile}" \
-iin="${sparseGRM_path}${sparseGRMFile}" \
-iin="${sparseGRM_path}${sparseGRMSampleIDFile}" \
-icmd="step1_fitNULLGLMM.R --plinkFile=${PLINK_for_vr} --nThreads=${nThreads} --outputPrefix=${outputPrefix} --useSparseGRMtoFitNULL=${useSparseGRMtoFitNULL} --relatednessCutoff=${relatednessCutoff} --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --isCateVarianceRatio=${isCateVarianceRatio} --phenoFile=${phenoFile} --phenoCol=${phenoCol} --covarColList=${covariatesList} --qCovarColList=${qCovarColList} --sampleIDColinphenoFile=${sampleIDCol} --traitType=${traitType} --invNormalize=${invNormalize} --tauInit=${tauInit}" \
-iimage_file="${PROJECT}:path/to/docker/saige_1.0.9.tar.gz" \
--name ${jobname} \
--instance-type ${instance_type} \
--priority ${priority}  \
--yes \
--destination ${destination}

### PRSprscsauto
instance_type="mem3_ssd1_v2_x4"
priority=normal
jobname=saige_step1_${TRAIT}_PRSprscsauto
nThreads=3
destination=exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/
traitType=quantitative
isCateVarianceRatio=TRUE
invNormalize=FALSE
useSparseGRMtoFitNULL=TRUE
relatednessCutoff=0.05
phenoCol=${TRAIT}
covariatesList=PRSprscsauto_${TRAIT},sex,enroll_age,enroll_age2,sequencing_batch,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20
qCovarColList=sex,sequencing_batch,genotyping_array
sampleIDCol=IID
PLINK_for_vr=ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01
PLINK_for_vr_path=${PROJECT}:path/to/array_data_for_saige/
sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
sparseGRM_path=${PROJECT}:path/to/SAIGE_step0_output_files/
phenoFile=phenofile_wes_total_quantitative_forBOLT.txt
phenoFile_path=${PROJECT}:path/to/phenofile/
phenoCol=${TRAIT}
outputPrefix=saige_step1_nullmodel_${TRAIT}_PRSprscsauto
tauInit=0,0 

dx run swiss-army-knife \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bed" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bim" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.fam" \
-iin="${phenoFile_path}${phenoFile}" \
-iin="${sparseGRM_path}${sparseGRMFile}" \
-iin="${sparseGRM_path}${sparseGRMSampleIDFile}" \
-icmd="step1_fitNULLGLMM.R --plinkFile=${PLINK_for_vr} --nThreads=${nThreads} --outputPrefix=${outputPrefix} --useSparseGRMtoFitNULL=${useSparseGRMtoFitNULL} --relatednessCutoff=${relatednessCutoff} --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --isCateVarianceRatio=${isCateVarianceRatio} --phenoFile=${phenoFile} --phenoCol=${phenoCol} --covarColList=${covariatesList} --qCovarColList=${qCovarColList} --sampleIDColinphenoFile=${sampleIDCol} --traitType=${traitType} --invNormalize=${invNormalize} --tauInit=${tauInit}" \
-iimage_file="${PROJECT}:path/to/docker/saige_1.0.9.tar.gz" \
--name ${jobname} \
--instance-type ${instance_type} \
--priority ${priority}  \
--yes \
--destination ${destination}


#####################
## In-sample scores
#####################

### PRSclumped_insample
instance_type="mem3_ssd1_v2_x4"
priority=normal
jobname=saige_step1_${TRAIT}_PRSclumped_insample
nThreads=3
destination=exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/
traitType=quantitative
isCateVarianceRatio=TRUE
invNormalize=FALSE
useSparseGRMtoFitNULL=TRUE
relatednessCutoff=0.05
phenoCol=${TRAIT}
covariatesList=PRSclumped_${TRAIT},sex,enroll_age,enroll_age2,sequencing_batch,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20
qCovarColList=sex,sequencing_batch,genotyping_array
sampleIDCol=IID
PLINK_for_vr=ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01
PLINK_for_vr_path=${PROJECT}:path/to/array_data_for_saige/
sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
sparseGRM_path=${PROJECT}:path/to/SAIGE_step0_output_files/
phenoFile=phenofile_wes_total_quantitative_forBOLT.txt
phenoFile_path=${PROJECT}:path/to/phenofile/
phenoCol=${TRAIT}
outputPrefix=saige_step1_nullmodel_${TRAIT}_PRSclumped_insample
tauInit=0,0

dx run swiss-army-knife \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bed" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bim" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.fam" \
-iin="${phenoFile_path}${phenoFile}" \
-iin="${sparseGRM_path}${sparseGRMFile}" \
-iin="${sparseGRM_path}${sparseGRMSampleIDFile}" \
-icmd="step1_fitNULLGLMM.R --plinkFile=${PLINK_for_vr} --nThreads=${nThreads} --outputPrefix=${outputPrefix} --useSparseGRMtoFitNULL=${useSparseGRMtoFitNULL} --relatednessCutoff=${relatednessCutoff} --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --isCateVarianceRatio=${isCateVarianceRatio} --phenoFile=${phenoFile} --phenoCol=${phenoCol} --covarColList=${covariatesList} --qCovarColList=${qCovarColList} --sampleIDColinphenoFile=${sampleIDCol} --traitType=${traitType} --invNormalize=${invNormalize} --tauInit=${tauInit}" \
-iimage_file="${PROJECT}:path/to/docker/saige_1.0.9.tar.gz" \
--name ${jobname} \
--instance-type ${instance_type} \
--priority ${priority}  \
--yes \
--destination ${destination}

### PRSprscsauto_insample
instance_type="mem3_ssd1_v2_x4"
priority=normal
jobname=saige_step1_${TRAIT}_PRSprscsauto_insample
nThreads=3
destination=exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/
traitType=quantitative
isCateVarianceRatio=TRUE
invNormalize=FALSE
useSparseGRMtoFitNULL=TRUE
relatednessCutoff=0.05
phenoCol=${TRAIT}
covariatesList=PRSprscsauto_${TRAIT},sex,enroll_age,enroll_age2,sequencing_batch,genotyping_array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20
qCovarColList=sex,sequencing_batch,genotyping_array
sampleIDCol=IID
PLINK_for_vr=ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01
PLINK_for_vr_path=${PROJECT}:path/to/array_data_for_saige/
sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
sparseGRM_path=${PROJECT}:path/to/SAIGE_step0_output_files/
phenoFile=phenofile_wes_total_quantitative_forBOLT.txt
phenoFile_path=${PROJECT}:path/to/phenofile/
phenoCol=${TRAIT}
outputPrefix=saige_step1_nullmodel_${TRAIT}_PRSprscsauto_insample
tauInit=0,0 

dx run swiss-army-knife \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bed" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.bim" \
-iin="${PLINK_for_vr_path}${PLINK_for_vr}.fam" \
-iin="${phenoFile_path}${phenoFile}" \
-iin="${sparseGRM_path}${sparseGRMFile}" \
-iin="${sparseGRM_path}${sparseGRMSampleIDFile}" \
-icmd="step1_fitNULLGLMM.R --plinkFile=${PLINK_for_vr} --nThreads=${nThreads} --outputPrefix=${outputPrefix} --useSparseGRMtoFitNULL=${useSparseGRMtoFitNULL} --relatednessCutoff=${relatednessCutoff} --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --isCateVarianceRatio=${isCateVarianceRatio} --phenoFile=${phenoFile} --phenoCol=${phenoCol} --covarColList=${covariatesList} --qCovarColList=${qCovarColList} --sampleIDColinphenoFile=${sampleIDCol} --traitType=${traitType} --invNormalize=${invNormalize} --tauInit=${tauInit}" \
-iimage_file="${PROJECT}:path/to/docker/saige_1.0.9.tar.gz" \
--name ${jobname} \
--instance-type ${instance_type} \
--priority ${priority}  \
--yes \
--destination ${destination}
