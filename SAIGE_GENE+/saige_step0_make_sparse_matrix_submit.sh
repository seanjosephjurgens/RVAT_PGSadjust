## Submission script to submit SAIGE step0 analysis using Swiss Army knife

PROJECT="RAP_PROJECT_NAME"

~/dx/dx-toolkit/bin/dx run swiss-army-knife \
-iin="${PROJECT}:path/to/array_data_for_saige/ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01.bed" \
-iin="${PROJECT}:path/to/array_data_for_saige/ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01.bim" \
-iin="${PROJECT}:path/to/array_data_for_saige/ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01.fam" \
-icmd="createSparseGRM.R --plinkFile=ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22_geno0.01 --nThreads=94  --outputPrefix=saige_sparse_matrix --numRandomMarkerforSparseKin=5000 --relatednessCutoff=0.05" \
-iimage_file="${RAP_PROJECT_NAME}:path/to/docker/saige_1.0.9.tar.gz" \
--name saige_step0 \
--instance-type mem2_ssd1_v2_x96 \
--priority normal \
--yes \
--destination ${RAP_PROJECT_NAME}:path/to/SAIGE_step0_output_files/
