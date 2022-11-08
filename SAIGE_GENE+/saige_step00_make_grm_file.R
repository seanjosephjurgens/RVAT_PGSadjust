#! Rscript

### SAIGE-GENE+ requires a PLINK file for step0 which includes very rare variants for fitting of the nullmodel. 
### Here, we will use rare variants from the exome sequencing data to add some random rare variation to the PLINK file.

# This script was run interactively using a Jupyter Notebook (R kernell) on the UKB RAP.

PROJECT="RAP_PROJECT_NAME"

system("mkdir tmp")
system("cd tmp")


# download plink and plink2
system("wget https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip")
system("unzip plink2_linux_avx2.zip")
system("wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip")
system("unzip plink_linux_x86_64_20220402.zip")


# Import genotyping array data previously pruned to ~280k array variants (see the 'genotype_array_file_forGRM' file from this repos)
system("dx download ${PROJECT}:path/to/pruned_array_data/ukbb200k_array_chrall.*")

# Import phenotype and covar file for QCd samples from the initial 200k freeze of UKB, which includes the out-of-sample derived PGS for all traits
system("dx download ${PROJECT}:path/to/phenofile/phenofile_wes_total_quantitative_forBOLT.txt")

# Import QCd exome sequencing plink data for the samples from the initial 200k release of UKB exomes; we will use chr22 data
system("dx download ${PROJECT}:path/to/200k_exome_data/ukb23156_c22_v1.genotype_variant_sample_QCed.*")

# Import double column ID file (FID and IID) with IDs for samples passing the QC for the genotyping array and exome sequencing for the current analyses
system("dx download ${PROJECT}:path/to/QC_passingIDs/ids_double_afterQC.txt")


#############################################
# Make adjusted GRM PLINK file for SAIGE-GENE
#############################################

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
system(paste0(' plink2 ',
              ' --bfile ukb23156_c22_v1.genotype_variant_sample_QCed ',
              ' --max-mac 2 --mac 2 ',
              ' --make-bed ',
              ' --exclude range range.tsv ',
              ' --keep ids_double_afterQC_steveIDs.txt ',
              ' --write-snplist --out ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_mac2_include '
))
write.table(paste0('22 1 36000000'), file='range.tsv', col.names=F, row.names=F, quote=F)
system(paste0(' plink2 ',
              ' --bfile ukb23156_c22_v1.genotype_variant_sample_QCed ',
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

system(paste0(' plink2 ',
              ' --bfile ukb23156_c22_v1.genotype_variant_sample_QCed ',
              ' --make-bed ',
              ' --extract keep_list.txt ',
              ' --keep ids_double_afterQC_steveIDs.txt ',
              ' --out ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter '
))
# LD prune
system(paste0(' plink2 ',
              ' --bfile ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter ',
              ' --indep-pairwise 500 50 0.2  ',
              ' --out plink2_wes '
))
system(paste0(' plink2 ',
              ' --bfile ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter ',
              ' --make-bed ',
              ' --exclude plink2_wes.prune.out ',
              ' --out ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter2 '
))


#Manually adjust variant positions so PLINK merge can handle multiallelic variants without realizing it
bim <- data.table::fread('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter2.bim', stringsAsFactors=F, header=F)
bim <- as.data.frame(bim)
bim2 <- data.table::fread('ukbb200k_array_chrall.bim', stringsAsFactors=F, header=F)
bim2 <- as.data.frame(bim2)
options(scipen=999)
for(chr in unique(bim2$V1)){
        cat('Busy with chromosome',chr,'...\n')
        rowlength <- nrow(bim2[bim2$V1 == chr,])

        bim2[bim2$V1 == chr, 'V4'] <- seq(1, rowlength, 1)
        if(chr %in% paste0(c(1:22))){
                rowlength2 <- nrow(bim[bim$V1 == chr,])
                if(rowlength2>0){
                        start <- rowlength+1
                        end <- rowlength + rowlength2
                        bim[bim$V1 == chr, 'V4'] <- seq(start, end, 1)
                }
        }
}
bim$key <- bim$V2
bim$V2 <- paste0("WES", seq(1, nrow(bim), 1))
bim2$key <- bim2$V2
bim2$V2 <- paste0("ARRAY", seq(1, nrow(bim2), 1))

write.table(bim[,c(1:6)], file='ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter2_fakepositions.bim', col.names=F, row.names=F, quote=F)
write.table(bim2[,c(1:6)], file='ukbb200k_array_chrall_fakepositions.bim', col.names=F, row.names=F, quote=F)
write.table(as.data.frame(rbind(bim[,c(7,2)],bim2[,c(7,2)])), file='fakepositions_key.tsv', sep='\t', col.names=F, row.names=F, quote=F)

### Merge the files
mergefile <- NULL
bed <- paste0('ukbb200k_array_chrall.bed')
bim <- paste0('ukbb200k_array_chrall_fakepositions.bim')
fam <- paste0('ukbb200k_array_chrall.fam')
line <- c(bed, bim, fam)
mergefile <- rbind(mergefile, line)
bed <- paste0('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter2.bed')
bim <- paste0('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter2_fakepositions.bim')
fam <- paste0('ukbb200k_wes_rarevariantsforGRM_QC_all_chr22_inter2.fam')
line <- c(bed, bim, fam)
mergefile <- rbind(mergefile, line)

write.table(mergefile, file='mergefile2.txt', col.names=F, row.names=F, quote=F)

#Merge requires a lot of memory
system(paste0(' plink ',
              ' --merge-list mergefile2.txt ',
              ' --out ../ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22 '
))

#Fix positions and variant names
bim <- data.table::fread('../ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22.bim', stringsAsFactors=F, header=F)
key <- data.table::fread('fakepositions_key.tsv', stringsAsFactors=F, header=F)
colnames(key)[1] <- "variant"
bim <- plyr::join(bim, key, by='V2', type='inner')
bim$V2 <- bim$variant
bim$V4 <- bim$variant
bim$V4 <- gsub('A', '', bim$V4)
bim$V4 <- gsub('T', '', bim$V4)
bim$V4 <- gsub('C', '', bim$V4)
bim$V4 <- gsub('G', '', bim$V4)
bim$V4 <- gsub('::', '', bim$V4)
bim$V4 <- gsub('.*:', '', bim$V4)

write.table(bim[,c(1:6)], file='../ukbb200k_all_commonvariantsfromarray_chrall_and_rarevariantsfromWES_chr22.bim', col.names=F, row.names=F, quote=F)
system('rm *inter* *prune.* range.tsv *fakepositions* *.log ')

system("cd ..")
system("rm -rf tmp/")
