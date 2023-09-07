################################################
#                                              #
#  title: 'Maternal miR-eQTL analyses in BiSC' #
#  date: '2023-07-23'                          #
#  Mariona Bustamante / Dolors Pelegrí /       #
#            Antonio Gómez                     #
#                                              #
################################################

# Sample selection
## Restrict analyses to European
## Non-imputed
## SNPs: call rate 0.97, MAF 5%, HWE 1E-06

################
# miRNA levels #
################

# Read miRNA expression (sample x miRNA) 

## Load matrix expression miRNA:
ex <- readRDS("Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20220803/denoising/BISC_miRNA_TMM_5cpm_10p_c.lib.rds")

## Load data.frame of correspondences ID_lab PLACENTA to ID child (id_bisc, XXXXXX31) and covariates
cov <- read.table("Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20220803/pheno_ddbb/pheno389_complete_rin.txt", sep = "\t", header = TRUE)
rownames(ex) <- cov$idparticipant[match(rownames(ex), cov$Id_lab)]

# Convert ID child (id_bisc, XXXXXX31) to ID Mother (XXXXXX11)
row.names(ex) <- sub("31$", "11", row.names(ex))

# Load data.frame of correspondences ID Mother (SubjectID) to LabID blood (BISC_XXXX)
bisc_names <- read.table(file = "Y:/analyses/BiSC_21/006_PlamiQTLs_AG/script/BISC_names.txt", header = TRUE)
rownames(ex) <- bisc_names$GWASID[match(rownames(ex), bisc_names$SubjectID)]

table(is.na(rownames(ex)))

# FALSE  TRUE 
#   383     6 

dim(ex)
#[1] 389 417 

ex <- ex[!is.na(rownames(ex)),]
dim(ex)
#[1] 383 417 <- 383 samples x 417 mirRNA

ex <- t(ex)


#	Read miRNA annotation (miRNA x annotation) 

annot <- read.table(file = "Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20230313/DBs/mirbase/miRNA_annotation_v22_GRCh37.txt")

library(stringr)
annot$mir_names <- str_extract(annot$V9, "Name=([^;]+)")
annot$mir_names <- gsub("Name=", "", annot$mir_names)
annot$mir_names <- gsub("-", ".", annot$mir_names)

rownames(ex)[!(rownames(ex) %in% annot$mir_names)]
#  [1] "hsa.let.7f.2.3p"  "hsa.miR.1299"     "hsa.miR.141.3p"   "hsa.miR.141.5p"   "hsa.miR.200c.3p" 
#  [6] "hsa.miR.200c.5p"  "hsa.miR.221.3p"   "hsa.miR.221.5p"   "hsa.miR.29a.3p"   "hsa.miR.335.3p"  
# [11] "hsa.miR.335.5p"   "hsa.miR.362.5p"   "hsa.miR.3664.3p"  "hsa.miR.500a.5p"  "hsa.miR.501.3p"  
# [16] "hsa.miR.501.5p"   "hsa.miR.506.3p"   "hsa.miR.508.3p"   "hsa.miR.509.3.5p" "hsa.miR.509.3p"  
# [21] "hsa.miR.509.5p"   "hsa.miR.511.3p"   "hsa.miR.511.5p"   "hsa.miR.514a.3p"  "hsa.miR.532.3p"  
# [26] "hsa.miR.532.5p"   "hsa.miR.590.3p"   "hsa.miR.660.5p"   "hsa.miR.98.5p"   

ex <- ex[(which(rownames(ex) %in% annot$mir_names)),]
dim(ex) #388 miRNA <- 417-388 = 29 miRNA in ex but not annoted (because we have lost them when converting GRCh38 to GRCh37 with Assembly Converter). 383 samples


# EXPRESION BED FILE

## chr, start, end, phenotype_id, samples (same ID as genotype input)
## start = end-1

ex_df <- as.data.frame(ex)

annot_fin <- annot[,c(1,4,10)]

annot_fin$stop <- annot_fin$V4+1
names(annot_fin) <- c("chr", "start","phenotype_id", "end")

ex_df_fin <- merge(ex_df, annot_fin, by.x="row.names", by.y="phenotype_id")

ex.bed <- ex_df_fin[,c(385:387, 1, 2:384)]
colnames(ex.bed)[1]<-"#chr"
colnames(ex.bed)[4]<-"phenotype_id"

## Paste miRNA id with chromosomic locations to obtain a unique gene (miRNA) identifier
ex.bed$phenotype_id <- paste(ex.bed$phenotype_id, ex.bed$"#chr", ex.bed$start, ex.bed$end, sep="_")


write.table(ex.bed, "BISC_miRNA.bed", 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
## ATENTION!!! The bed file needs to be gzipped (i.e., the file in TensorQTL should be of the type expression.bed.gz)
## gzip -c BISC_miRNA.bed > BISC_miRNA.bed.gz

##################
# COVARIATE FILE #
##################

# PCA only European subsample

pcs <- read.table("Y:/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/PCA/bisc_m_snpqc_maf5_EUR_inclsex_PC.eigenvec")
pcs <- pcs[,2:7]


# Merge covariates

cov$Id_child <- sub("31$", "11", cov$idparticipant)
cov.bisc <- merge(bisc_names, cov, by.x = "SubjectID", by.y = "Id_child")

covTQTL <- cov.bisc[ , c("GWASID", "b_sexo", "gestage_0y_c_days", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC")] #, "Syncytiotrophoblast"

covTQTL <- merge(covTQTL, pcs, by.x = "GWASID", by.y = "V2")
dim(covTQTL)
# [1] 271  13 <- EUR subsample

names(covTQTL) <-  c("GWASID", "b_sexo", "gestage_0y_c_days", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "PC1", "PC2", "PC3", "PC4", "PC5")


# Recode qualitative covariates with integers

covTQTL$b_sexo[covTQTL$b_sexo=="Boy"] <- 2
covTQTL$b_sexo[covTQTL$b_sexo=="Girl"] <- 1

covTQTL_t <- t(covTQTL)

write.table(covTQTL_t, "BISC_Cells_covariates.txt", 
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)


# Here keep samples included if needed to match with samples in GWAS

write.table(covTQTL_t[1,],"samples_to_filter_Cells.txt", 
            sep=",",col.names = FALSE, row.names=FALSE)





