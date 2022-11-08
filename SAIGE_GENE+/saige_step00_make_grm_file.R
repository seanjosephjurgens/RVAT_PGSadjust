#! Rscript

### SAIGE-GENE+ requires a PLINK file for step0 which includes very rare variants for fitting of the nullmodel. 
### Here, we will use rare variants from the exome sequencing data to add some random rare variation to the PLINK file.



# Change variant names in bim file
system(paste0(' rm *array_chr*_inter*bed *array_chr*_inter*bim *array_chr*_inter*fam *_inter*log '))
bimm <- fread('ukbb200k_array_chrall.bim', stringsAsFactors=F, data.table=F)
bimm$V2 <- paste0(bimm$V1, ":", bimm$V4, ":", bimm$V5, ":", bimm$V6)
write.table(bimm, file='ukbb200k_array_chrall.bim', col.names=F, row.names=F, quote=F, sep='\t')
# Change IDs in fam file
famm <- fread('ukbb200k_array_chrall.fam', stringsAsFactors=F, data.table=F)
nrow(famm)
pheno <- fread('/medpop/afib/projects/ukbb_200k_exome/cardiometabolic/V2/PRS_adjust/V2/data/pheno/EUR/out_sample/phenofile_wes_total_quantitative.tsv', stringsAsFactors=F)
pheno <- pheno[,c(1:2)]
all(famm$V1 == pheno$IID)
#TRUE
famm$V1 <- famm$V2 <- pheno$sample_id
write.table(famm, file='ukbb200k_array_chrall.fam', col.names=F, row.names=F, quote=F)

###### Make file with further pruning for experimenting in BOLT#####
system(paste0(' /medpop/afib/chaffin/PLINK2/07.19.2018/plink2 ',
         ' --bfile ukbb200k_array_chrall ',
         ' --indep-pairwise 250 100 0.05 ',
         ' --out plink2_bolt '
))
system(paste0(' /medpop/afib/chaffin/PLINK2/_10.09.2017/plink2 ',
              ' --bfile ukbb200k_array_chrall ',
              ' --make-bed ',
              ' --maf 0.01 --geno 0.01 ',
              ' --exclude plink2_bolt.prune.out ',
              ' --out ukbb200k_array_secondpruneround_chrall '
))
