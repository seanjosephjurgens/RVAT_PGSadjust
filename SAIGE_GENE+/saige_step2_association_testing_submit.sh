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
  
for chr in {1..22};
do
  
  #####################
  ## No PGS added
  #####################

  ### noPRS
  instance_type="mem3_ssd1_v2_x32"
  priority=normal
  jobname=saige_step2_${TRAIT}_${chr}_noPRS
  nThreads=30
  destination=exome-seq:sjj/short_projects/PRSadjust/R2/results/SAIGE/chr${chr}
  bedFile=${PROJECT}:path/to/gene_burdens_file/plink_format/ukbb200k_wesQC_chr${chr}_bed.bed
  bimFile=${PROJECT}:path/to/gene_burdens_file/plink_format/ukbb200k_wesQC_chr${chr}_bed.bim
  famFile=${PROJECT}:path/to/gene_burdens_file/plink_format/ukbb200k_wesQC_chr${chr}_bed.fam
  AlleleOrder=ref-first
  SAIGEOutputFile=saige_step2_association_${TRAIT}_chr${chr}_noPRS
  sampleFile=${PROJECT}:path/to/QC_passingIDs/ids_afterQC.txt
  nullmod_prefix=${PROJECT}:path/to/SAIGE_nullmodels/saige_step1_nullmodel_${TRAIT}_noPRS
  sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
  sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
  groupFile=SAIGE_GENE_groupingfile_chr${chr}_SAIGE_GENEplus_bed.txt
  maxMAC_in_groupTest=40
  maxMAF_in_groupTest=0.0001
  #minGroupMAC_in_BurdenTest=20
  LOCO=FALSE
  is_fastTest=TRUE
  is_output_markerList_in_groupTest=FALSE
  is_single_in_groupTest=FALSE
  is_no_weight_in_groupTest=FALSE
  annotation_in_groupTest=LOFmiss
  
  ~/dx/dx-toolkit/bin/dx run swiss-army-knife \
  -iin="${bedFile}" \
  -iin="${bimFile}" \
  -iin="${famFile}" \
  -iin="${sampleFile}" \
  -iin="${nullmod_prefix}.rda" \
  -iin="${nullmod_prefix}.varianceRatio.txt" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMSampleIDFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/grouping_files/${groupFile}" \
  -icmd="step2_SPAtests.R --bedFile=${bedFile} --bimFile=${bimFile} --famFile=${famFile} --AlleleOrder=${AlleleOrder} --chr=${chr} --SAIGEOutputFile=${SAIGEOutputFile} --minMAF=0 --minMAC=0.5 --sampleFile=${sampleFile} --GMMATmodelFile=${nullmod_prefix}.rda --varianceRatioFile=${nullmod_prefix}.varianceRatio.txt --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --groupFile=${groupFile} --maxMAC_in_groupTest=${maxMAC_in_groupTest} --maxMAF_in_groupTest=${maxMAF_in_groupTest}  --is_output_markerList_in_groupTest=${is_output_markerList_in_groupTest} --LOCO=${LOCO} --is_fastTest=${is_fastTest} --is_single_in_groupTest=${is_single_in_groupTest} --is_no_weight_in_groupTest=${is_no_weight_in_groupTest} --annotation_in_groupTest=${annotation_in_groupTest}" \
  -iimage_file="exome-seq:sjj/docker/saige_1.0.9.tar.gz" \
  --name ${jobname} \
  --instance-type ${instance_type} \
  --priority ${priority}  \
  --yes \
  --destination ${destination}
  
  ### PRSclumped
  instance_type="mem3_ssd1_v2_x32"
  priority=normal
  jobname=saige_step2_${TRAIT}_${chr}_PRSclumped
  nThreads=30
  destination=exome-seq:sjj/short_projects/PRSadjust/R2/results/SAIGE/chr${chr}
  bedFile=ukbb200k_wesQC_chr${chr}_bed.bed
  bimFile=ukbb200k_wesQC_chr${chr}_bed.bim
  famFile=ukbb200k_wesQC_chr${chr}_bed.fam
  AlleleOrder=ref-first
  SAIGEOutputFile=saige_step2_association_${TRAIT}_chr${chr}_PRSclumped
  sampleFile=ids_afterQC_steveIDs.txt
  nullmod_prefix=saige_step1_nullmodel_${TRAIT}_PRSclumped
  sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
  sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
  groupFile=SAIGE_GENE_groupingfile_chr${chr}_SAIGE_GENEplus_bed.txt
  maxMAC_in_groupTest=40
  maxMAF_in_groupTest=0.0001
  #minGroupMAC_in_BurdenTest=20
  LOCO=FALSE
  is_fastTest=TRUE
  is_output_markerList_in_groupTest=FALSE
  is_single_in_groupTest=FALSE
  is_no_weight_in_groupTest=FALSE
  annotation_in_groupTest=LOFmiss
  
  ~/dx/dx-toolkit/bin/dx run swiss-army-knife \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bedFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bimFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${famFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${sampleFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.rda" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.varianceRatio.txt" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMSampleIDFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/grouping_files/${groupFile}" \
  -icmd="step2_SPAtests.R --bedFile=${bedFile} --bimFile=${bimFile} --famFile=${famFile} --AlleleOrder=${AlleleOrder} --chr=${chr} --SAIGEOutputFile=${SAIGEOutputFile} --minMAF=0 --minMAC=0.5 --sampleFile=${sampleFile} --GMMATmodelFile=${nullmod_prefix}.rda --varianceRatioFile=${nullmod_prefix}.varianceRatio.txt --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --groupFile=${groupFile} --maxMAC_in_groupTest=${maxMAC_in_groupTest} --maxMAF_in_groupTest=${maxMAF_in_groupTest}  --is_output_markerList_in_groupTest=${is_output_markerList_in_groupTest} --LOCO=${LOCO} --is_fastTest=${is_fastTest} --is_single_in_groupTest=${is_single_in_groupTest} --is_no_weight_in_groupTest=${is_no_weight_in_groupTest} --annotation_in_groupTest=${annotation_in_groupTest}" \
  -iimage_file="exome-seq:sjj/docker/saige_1.0.9.tar.gz" \
  --name ${jobname} \
  --instance-type ${instance_type} \
  --priority ${priority}  \
  --yes \
  --destination ${destination}
  
  ### PRSprscsauto
  instance_type="mem3_ssd1_v2_x32"
  priority=normal
  jobname=saige_step2_${TRAIT}_${chr}_PRSprscsauto
  nThreads=30
  destination=exome-seq:sjj/short_projects/PRSadjust/R2/results/SAIGE/chr${chr}
  bedFile=ukbb200k_wesQC_chr${chr}_bed.bed
  bimFile=ukbb200k_wesQC_chr${chr}_bed.bim
  famFile=ukbb200k_wesQC_chr${chr}_bed.fam
  AlleleOrder=ref-first
  SAIGEOutputFile=saige_step2_association_${TRAIT}_chr${chr}_PRSprscsauto
  sampleFile=ids_afterQC_steveIDs.txt
  nullmod_prefix=saige_step1_nullmodel_${TRAIT}_PRSprscsauto
  sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
  sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
  groupFile=SAIGE_GENE_groupingfile_chr${chr}_SAIGE_GENEplus_bed.txt
  maxMAC_in_groupTest=40
  maxMAF_in_groupTest=0.0001
  #minGroupMAC_in_BurdenTest=20
  LOCO=FALSE
  is_fastTest=TRUE
  is_output_markerList_in_groupTest=FALSE
  is_single_in_groupTest=FALSE
  is_no_weight_in_groupTest=FALSE
  annotation_in_groupTest=LOFmiss
  
  ~/dx/dx-toolkit/bin/dx run swiss-army-knife \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bedFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bimFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${famFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${sampleFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.rda" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.varianceRatio.txt" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMSampleIDFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/grouping_files/${groupFile}" \
  -icmd="step2_SPAtests.R --bedFile=${bedFile} --bimFile=${bimFile} --famFile=${famFile} --AlleleOrder=${AlleleOrder} --chr=${chr} --SAIGEOutputFile=${SAIGEOutputFile} --minMAF=0 --minMAC=0.5 --sampleFile=${sampleFile} --GMMATmodelFile=${nullmod_prefix}.rda --varianceRatioFile=${nullmod_prefix}.varianceRatio.txt --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --groupFile=${groupFile} --maxMAC_in_groupTest=${maxMAC_in_groupTest} --maxMAF_in_groupTest=${maxMAF_in_groupTest}  --is_output_markerList_in_groupTest=${is_output_markerList_in_groupTest} --LOCO=${LOCO} --is_fastTest=${is_fastTest} --is_single_in_groupTest=${is_single_in_groupTest} --is_no_weight_in_groupTest=${is_no_weight_in_groupTest} --annotation_in_groupTest=${annotation_in_groupTest}" \
  -iimage_file="exome-seq:sjj/docker/saige_1.0.9.tar.gz" \
  --name ${jobname} \
  --instance-type ${instance_type} \
  --priority ${priority}  \
  --yes \
  --destination ${destination}

  ### PRSclumped_insample
  instance_type="mem3_ssd1_v2_x32"
  priority=normal
  jobname=saige_step2_${TRAIT}_${chr}_PRSclumped_insample
  nThreads=30
  destination=exome-seq:sjj/short_projects/PRSadjust/R2/results/SAIGE/chr${chr}
  bedFile=ukbb200k_wesQC_chr${chr}_bed.bed
  bimFile=ukbb200k_wesQC_chr${chr}_bed.bim
  famFile=ukbb200k_wesQC_chr${chr}_bed.fam
  AlleleOrder=ref-first
  SAIGEOutputFile=saige_step2_association_${TRAIT}_chr${chr}_PRSclumped_insample
  sampleFile=ids_afterQC_steveIDs.txt
  nullmod_prefix=saige_step1_nullmodel_${TRAIT}_PRSclumped_insample
  sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
  sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
  groupFile=SAIGE_GENE_groupingfile_chr${chr}_SAIGE_GENEplus_bed.txt
  maxMAC_in_groupTest=40
  maxMAF_in_groupTest=0.0001
  #minGroupMAC_in_BurdenTest=20
  LOCO=FALSE
  is_fastTest=TRUE
  is_output_markerList_in_groupTest=FALSE
  is_single_in_groupTest=FALSE
  is_no_weight_in_groupTest=FALSE
  annotation_in_groupTest=LOFmiss
  
  ~/dx/dx-toolkit/bin/dx run swiss-army-knife \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bedFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bimFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${famFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${sampleFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.rda" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.varianceRatio.txt" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMSampleIDFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/grouping_files/${groupFile}" \
  -icmd="step2_SPAtests.R --bedFile=${bedFile} --bimFile=${bimFile} --famFile=${famFile} --AlleleOrder=${AlleleOrder} --chr=${chr} --SAIGEOutputFile=${SAIGEOutputFile} --minMAF=0 --minMAC=0.5 --sampleFile=${sampleFile} --GMMATmodelFile=${nullmod_prefix}.rda --varianceRatioFile=${nullmod_prefix}.varianceRatio.txt --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --groupFile=${groupFile} --maxMAC_in_groupTest=${maxMAC_in_groupTest} --maxMAF_in_groupTest=${maxMAF_in_groupTest}  --is_output_markerList_in_groupTest=${is_output_markerList_in_groupTest} --LOCO=${LOCO} --is_fastTest=${is_fastTest} --is_single_in_groupTest=${is_single_in_groupTest} --is_no_weight_in_groupTest=${is_no_weight_in_groupTest} --annotation_in_groupTest=${annotation_in_groupTest}" \
  -iimage_file="exome-seq:sjj/docker/saige_1.0.9.tar.gz" \
  --name ${jobname} \
  --instance-type ${instance_type} \
  --priority ${priority}  \
  --yes \
  --destination ${destination}
  
  ### PRSprscsauto
  instance_type="mem3_ssd1_v2_x32"
  priority=normal
  jobname=saige_step2_${TRAIT}_${chr}_PRSprscsauto_insample
  nThreads=30
  destination=exome-seq:sjj/short_projects/PRSadjust/R2/results/SAIGE/chr${chr}
  bedFile=ukbb200k_wesQC_chr${chr}_bed.bed
  bimFile=ukbb200k_wesQC_chr${chr}_bed.bim
  famFile=ukbb200k_wesQC_chr${chr}_bed.fam
  AlleleOrder=ref-first
  SAIGEOutputFile=saige_step2_association_${TRAIT}_chr${chr}_PRSprscsauto_insample
  sampleFile=ids_afterQC_steveIDs.txt
  nullmod_prefix=saige_step1_nullmodel_${TRAIT}_PRSprscsauto_insample
  sparseGRMFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
  sparseGRMSampleIDFile=saige_sparse_matrix_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt
  groupFile=SAIGE_GENE_groupingfile_chr${chr}_SAIGE_GENEplus_bed.txt
  maxMAC_in_groupTest=40
  maxMAF_in_groupTest=0.0001
  #minGroupMAC_in_BurdenTest=20
  LOCO=FALSE
  is_fastTest=TRUE
  is_output_markerList_in_groupTest=FALSE
  is_single_in_groupTest=FALSE
  is_no_weight_in_groupTest=FALSE
  annotation_in_groupTest=LOFmiss
  
  ~/dx/dx-toolkit/bin/dx run swiss-army-knife \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bedFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${bimFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${famFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/step2_genetic/carrier_files/${sampleFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.rda" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/nullmodels/${nullmod_prefix}.varianceRatio.txt" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/${sparseGRMSampleIDFile}" \
  -iin="exome-seq:sjj/short_projects/PRSadjust/R2/data/SAIGE_files/grouping_files/${groupFile}" \
  -icmd="step2_SPAtests.R --bedFile=${bedFile} --bimFile=${bimFile} --famFile=${famFile} --AlleleOrder=${AlleleOrder} --chr=${chr} --SAIGEOutputFile=${SAIGEOutputFile} --minMAF=0 --minMAC=0.5 --sampleFile=${sampleFile} --GMMATmodelFile=${nullmod_prefix}.rda --varianceRatioFile=${nullmod_prefix}.varianceRatio.txt --sparseGRMFile=${sparseGRMFile} --sparseGRMSampleIDFile=${sparseGRMSampleIDFile} --groupFile=${groupFile} --maxMAC_in_groupTest=${maxMAC_in_groupTest} --maxMAF_in_groupTest=${maxMAF_in_groupTest} --is_output_markerList_in_groupTest=${is_output_markerList_in_groupTest} --LOCO=${LOCO} --is_fastTest=${is_fastTest} --is_single_in_groupTest=${is_single_in_groupTest} --is_no_weight_in_groupTest=${is_no_weight_in_groupTest} --annotation_in_groupTest=${annotation_in_groupTest}" \
  -iimage_file="exome-seq:sjj/docker/saige_1.0.9.tar.gz" \
  --name ${jobname} \
  --instance-type ${instance_type} \
  --priority ${priority}  \
  --yes \
  --destination ${destination}
  
  done
