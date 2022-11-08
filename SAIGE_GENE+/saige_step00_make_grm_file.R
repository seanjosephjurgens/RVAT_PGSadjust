#! Rscript

### SAIGE-GENE+ requires a PLINK file for step0 which includes very rare variants for fitting of the nullmodel. 
### Here, we will use rare variants from the exome sequencing data to add some random rare variation to the PLINK file.

# This script was run interactively using a Jupyter Notebook (R kernell) on the UKB RAP.

PROJECT="RAP_PROJECT_NAME"

# Import genotyping array data previously pruned to ~280k array variants (see the 'genotype_array_file_forGRM' file from this repos)
system("dx download ${PROJECT}:path/to/pruned_array_data/ukbb200k_array_chrall.*")

# Import phenotype and covar file for QCd samples from the initial 200k freeze of UKB, which includes the out-of-sample derived PGS for all traits
system("dx download ${PROJECT}:path/to/phenofile/phenofile_wes_total_quantitative_forBOLT.txt")

# Import QCd exome sequencing plink data for the samples from the initial 200k release of UKB exomes; we will use chr22 data
system("dx download ${PROJECT}:path/to/200k_exome_data/ukb23156_c22_v1.genotype_variant_sample_QCed.*")

# Import double column ID file (FID and IID) with IDs for samples passing the QC for the genotyping array and exome sequencing for the current analyses
system("dx download ${PROJECT}:path/to/ids_double_afterQC.txt

# Change variant names in bim file
bimm <- fread('ukbb200k_array_chrall.bim', stringsAsFactors=F, data.table=F)
bimm$V2 <- paste0(bimm$V1, ":", bimm$V4, ":", bimm$V5, ":", bimm$V6)
write.table(bimm, file='ukbb200k_array_chrall.bim', col.names=F, row.names=F, quote=F, sep='\t')


###Take independent rare variant subset for addition into GRM file, we will use only chr22 as this comput less expensive
# extract rare variants maf<0.001; for variance ratio estimation we will take additional QC measures: Variant missingnes <1%, hwe p < 1e-6, chr22
write.table(paste0('22 20000000 2000000000 1'), file='range.tsv', col.names=F, row.names=F, quote=F)
system(paste0(' plink2 ',
              ' --bfile ukb23156_c22_v1.genotype_variant_sample_QCed ',
              ' --max-mac 1 ',
              ' --make-bed ',
              ' --exclude range range.tsv ',
              ' --keep ids_double_afterQC.txt ',
              ' --write-snplist --out ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac1_include '
))
write.table(rbind(paste0('22 1 20000000'),paste0('22 36000000 2000000000')), file='range.tsv', col.names=F, row.names=F, quote=F)
system(paste0(' /medpop/afib/chaffin/PLINK2/07.19.2018/plink2 ',
              ' --bfile /medpop/afib/data/ukb9389/exome/200K/genotype_variant_sample_QCed/plink/ukb23156_c22_v1.genotype_variant_sample_QCed ',
              ' --max-mac 2 --mac 2 ',
              ' --make-bed ',
              ' --exclude range range.tsv ',
              ' --keep ids_double_afterQC_steveIDs.txt ',
              ' --write-snplist --out ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac2_include '
))
write.table(paste0('22 1 36000000'), file='range.tsv', col.names=F, row.names=F, quote=F)
system(paste0(' /medpop/afib/chaffin/PLINK2/07.19.2018/plink2 ',
              ' --bfile /medpop/afib/data/ukb9389/exome/200K/genotype_variant_sample_QCed/plink/ukb23156_c22_v1.genotype_variant_sample_QCed ',
              ' --max-mac 20 --mac 3 ',
              ' --make-bed ',
              ' --exclude range range.tsv ',
              ' --keep ids_double_afterQC_steveIDs.txt ',
              ' --write-snplist --out ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac3_20_include '
))
f1 <- fread('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac1_include.snplist', stringsAsFactors=F, data.table=F, header=F)
f2 <- fread('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac2_include.snplist', stringsAsFactors=F, data.table=F, header=F)
f3 <- fread('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac3_20_include.snplist', stringsAsFactors=F, data.table=F, header=F)
f <- rbind(f1, f2, f3)
write.table(f, file='keep_list.txt', col.names=F, row.names=F, quote=F)
system('rm *include*')

