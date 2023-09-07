#####################################
#                                   #
# title: 'Graphics'                 #
# date: '2023-08-10'                #
#                                   #
#####################################

#############
# Histogram #
#############

#Load data
cov <- read.table("Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20220803/pheno_ddbb/pheno389_complete_rin.txt", sep = "\t", header = TRUE)

# Variables
variables <- c("Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "Syncytiotrophoblast")


# Number of rows and columns to organise histograms
rows <- 2
cols <- 3

# Set the grid layout
par(mfrow = c(rows, cols))

# Letter size
font_size <- 1.5


# Creating histograms individually and organising them
for (var in variables) {
  hist(cov[[var]], main = var, xlab = "", ylab = "Frecuency",
       col = "steelblue", border = "black", cex.main = font_size, cex.lab = 1.1, cex.axis = 1)
  
  # Change the text on the x-axis
  axis(side = 1, labels = FALSE)  # Suppresses the original x-axis labels
  mtext("Cell Proportions", side = 1, line = 2, cex = 0.8)  # Adds the new text on the x-axis
}

# Restore the grid layout to its original state
par(mfrow = c(1, 1))

# Other
pdf("./sample_qc/bisc_m_snpqc97_SNPcallrate.pdf")
hist(l$cr, breaks=50, main="Histogram of F_MISS Values - snp call rate", xlab="F_MISS")
dev.off()


###########
# Corplot #
###########

# Load data
library(readxl)
pheno <- read_excel("Y:/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/db/BISC_db_gwas_20230506_codebook.xlsx", 
                    sheet = "db_miQTL_m_EUR", col_types = c("text", 
                                                            "text", "text", "numeric", "numeric", 
                                                            "text", "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric"))
# Select variables
covTQTL <- pheno[ , c("gestational_age_birth", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC", "Syncytiotrophoblast")] 

# Load corrplot package
library(corrplot)

# Crear una matriz de correlación de ejemplo
matriz_cor <- cor(covTQTL[, c("gestational_age_birth", "Trophoblasts", "Stromal", "Hofbauer", "Endothelial", "nRBC","Syncytiotrophoblast")], method = "spearman")

colnames(matriz_cor)[1] <- "Gestational Age"
rownames(matriz_cor)[1] <- "Gestational Age"

# Generate the corplot
corrplot(matriz_cor, method = "circle", tl.cex = 0.7, type = "upper")

corrplot.mixed(matriz_cor, upper = "circle", diag = "n", tl.pos = "d", tl.cex = 0.7, number.cex = 1.8)

library(PerformanceAnalytics)
chart.Correlation(covTQTL, histogram = TRUE, method = "pearson")
chart.Correlation(covTQTL, histogram = TRUE, method = "spearman")


###########
# Heatmap #
###########

V <- apply(ex, 1, var)

selectedGenes <- names(V[order(V, decreasing = T)])

#BiocManager::install("pheatmap")
library(pheatmap)
pheatmap(ex[selectedGenes,], scale = 'row', show_rownames = FALSE)


################
# Volcano plot #
################

library(readr)
data_v <- read_delim("Y:/analyses/BiSC_21/006_PlamiQTLs_AG/results/mat_miQTL_TensorQTL/bisc_m_TensorQTL-bim_pvalue1E-05.txt", 
                     delim = "\t", escape_double = FALSE)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(data_v,
                lab = data_v$rsID,
                labSize = 3.5,
                x = 'Beta',
                y = 'pval',
                drawConnectors = TRUE)

#In server

# Generate the graphic in PNG format
CairoPNG("VolcanoPlot.png", width = 1600, height = 1200, dpi = 300)

EnhancedVolcano(merged_df,
                lab = merged_df$rsID,
                labSize = 3.5,
                x = 'Beta',
                y = 'pval',
                drawConnectors = TRUE)
dev.off()


##################
# Manhattan plot #
##################

library(qqman)  
data_v$P <- data_v$pval
data_v$SNP <- data_v$rsID

manhattan(data_v, main = "", 
          col = c("blue4", "orange3"), 
          suggestiveline = -log10(1e-05), 
          genomewideline = -log10(5e-08), 
          annotateTop = TRUE, 
          annotatePval = -log10(5e-08), 
          chrlabs = as.character(sort(unique(data_v$CHR))))


# In server
library(Cairo)

# Generate the graphic in PNG format
CairoPNG("ManhattanPlot.png", width = 1600, height = 1200, dpi = 300)

# Create the Manhattan chart
manhattan(merged_df, main = "",
          chr="CHR", bp="BP", snp="rsID", p="pval",
          col = c("blue4", "orange3"),
          suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08),
          annotateTop = TRUE,
          annotatePval = -log10(5e-08),
          chrlabs = as.character(sort(unique(merged_df$CHR))))

# Close the Cairo PNG device
dev.off()

# Other Methods
library(tidyverse)
library(ggrepel)

#source("https://raw.githubusercontent.com/isglobal-brge/book_omic_association/master/R/manhattanPlot.R")

setwd("Y:/analyses/BiSC_21/006_PlamiQTLs_AG/results/mat_miQTL_TensorQTL") #to modify

library(data.table)
res<-fread("bisc_m_TensorQTL_pvalue1E-05_Rformated.txt")
dim(res)
names(res)

pvals <- data.frame(SNP=res$rsID, 
                    CHR=res$CHR,
                    BP=res$BP,
                    P=res$pval)
# missing data is not allowed
pvals <- subset(pvals, !is.na(CHR) & !is.na(P)) 

#install.packages("qqman")
library(qqman)

bonf.sig <- 1e-5 #5e-8

rsIDs <- c("rs2072474", "rs1397755", "GSA-rs6985501", "GSA-rs12375547", "rs61992671")
snpsOfInterest <- res[res$pval<5E-8, ]

plt <- manhattan(pvals, chr="CHR", bp="BP", snp="SNP", p="P", highlight = snpsOfInterest$rsID, suggestiveline = -log10(bonf.sig), main = "Figure 8. Manhattan plot") 


###################
# Locus Zoom plot #
###################

#Create a table including the SNPs with a p-value <1e-4, its chromosome, genomic position, minor allele, minor allele frequency, the annotated gene and the OR under dominant, recessive and additive model. 
#(HINT: use odds function in SNPassoc package). NOTE: If there are no SNPs that pass this significant level use another threshold.

#BiocManager::install("biomaRt")
library(biomaRt)
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host="https://www.ensembl.org")

snpInfo <- getBM(attributes=c("refsnp_id", "chr_name", "chrom_start", "allele", "minor_allele",  "minor_allele_freq", "ensembl_gene_stable_id"),
                 filters = c("snp_filter"), 
                 values = snpsOfInterest$rsID, 
                 mart = snpmart)
head(snpInfo)

tableSNPs <- merge(snpInfo, topPvals, by.x = "refsnp_id", by.y = "SNP")

tableSNPs <- cbind(tableSNPs, odds(WGa, model=c("dominant")), odds(WGa, model=c("recessive")), odds(WGa, model=c("log-additive")))

knitr::kable(tableSNPs, caption = "Table 7. significant SNPs info") 

#Create a Locus Zoom plot for those SNPs that are significantly associated after Bonferroni correction (top 5 as maximum).
#install.packages("topr")
library(topr)
# CHROM,POS,P 
names(tableSNPs)[8] = "CHROM"
names(tableSNPs)[9] = "POS"
names(tableSNPs)[10] = "P"

R <- read_excel("BISC_pla_miQTL_res_20230712.xlsx",sheet = "snps")
R <- R[,c("rsID","r2")]
res_R <- merge(res, R, by = "rsID")
names(res_R)[2] = "CHROM"
names(res_R)[3] = "POS"
names(res_R)[10] = "P"
names(res_R)[14] = "R2"

res_R_2 <- res_R[res_R$CHR==2,]

locuszoom(res_R_2)

################
# Density plot #
################

# The affy library has a density plotting function
ex <- read.table("Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20220803/QC_MB_20230701/db/BISC_pla_miRNA_qc.s.m_norm_lcpm_20230701.txt", header = TRUE)

#BiocManager::install("affy")
library(affy)

# Average of the genes
ex_mean <- rowMeans(ex)

## Create a list of 4 colors to use which are the same used throughout this chapter 
library(scales)
myColors <- hue_pal()(4)


## Plot the log2-transformed data
plotDensity(as.matrix(ex_mean), col='red',
            xlab='log2(normalised CPM)',
            main='Expression Distribution',
            lwd=4)

## Plot the log2-transformed data
plotDensity(ex, col=rep(myColors, each=3),
            lty=c(1:ncol(ex)),
            xlab='log2(normalised CPM)',
            main='Expression Distribution')

## Add a legend and vertical line
legend('topright', rownames(ex),
       col=rep(myColors, each=3))
abline(v=-1.5, lwd=1, col='red', lty=2)


############
# Box plot #
############

library(ggplot2)

ex <- read.table("Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20220803/QC_MB_20230701/db/BISC_pla_miRNA_qc.s.m_norm_lcpm_20230701.txt", header = TRUE)

# Suppose expresion_genes is your expression matrix
# You can replace this with your own data

# 1. Calculate the expression average for each gene
average_expression <- rowMeans(expresion_genes)

# 2. Sort the genes based on their average expression
sorted_genes <- order(average_expression)

# 3. Select the top 10 and bottom 10 most expressed genes
top_10_genes <- sorted_genes[(length(sorted_genes) - 9):length(sorted_genes)]
bottom_10_genes <- sorted_genes[1:10]

# 4. Prepare the data for the boxplot
boxplot_data <- data.frame(
  Gene = rep(c("Top 10", "Bottom 10"), each = 10),
  Expression = c(average_expression[top_10_genes], average_expression[bottom_10_genes])
)

# 5. Create the boxplot using ggplot2
ggplot(boxplot_data, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  labs(title = "Expression of the top 10 and bottom 10 expressed genes",
       x = "Genes",
       y = "Average Expression")

################
# Scatter plot #
################

# Check sex
pdf("bisc_m_snpqc_sexcheck.pdf")
cc <- ifelse(pheno$V5==1, "red", "blue")
plot(info.x$F, col=cc, 
     pch=16, xlab="Individuals", 
     ylab="Heterozygosity in chromosome X",
     main = "Figure 4. Sex Discrepancies")
legend("topleft", c("male", "female"), col=c("red", "blue"),
       pch=16)
abline(h = 0.2, col = "red")
abline(h = 0.8, col = "red")
dev.off()

# Heterozygosity
pdf("bisc_m_snpqc_het.pdf")
plot(info.het$F, pch = 16, xlab = "Index", ylab = "Heterozygosity", main = "Figure 5. Heterozigosity")

mean_value <- mean(info.het$F)
sd_value <- sd(info.het$F)
threshold1 <- mean_value - 4 * sd_value
threshold2 <- mean_value + 4 * sd_value

abline(h = threshold1, col = "red")
abline(h = threshold2, col = "red")

dev.off()

# PCA

library(ggplot2)

setwd("/PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/PCA/")

# read metadata
dd <- read.table("/PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/db/BISC_db_gwas_20230510_filter_ancestry.txt", header=T, sep="\t")
dim(dd)
head(dd)
dd1 <- dd[,c(1,2,7,8,11)]

rm(pcs)
pcs <- read.table("bisc_m_snpqc_maf1_inclsex_PC.eigenvec")
names(pcs) <- c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
p <- merge(pcs,dd1, by="IID")
dim(p) #949

pdf("ALL_inclsex_PCA-plink_PC1vsPC2_cohort.pdf")
qplot(PC1, PC2, colour = hospital_parto, 
      data = p)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC3vsPC4_cohort.pdf")
qplot(PC3, PC4, colour = hospital_parto, 
      data = p)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC1vsPC2_ancestry.pdf")
qplot(PC1, PC2, colour = ancestry_Grafpop,
      data = p)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC3vsPC4_ancestry.pdf")
qplot(PC3, PC4, colour = ancestry_Grafpop,
      data = p)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC1vsPC2_ethn_quest.pdf")
qplot(PC1, PC2, colour = ancestry_quest,
      data = p)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC3vsPC4_ethn_quest.pdf")
qplot(PC3, PC4, colour = ancestry_quest,
      data = p)
dev.off()



pdf("ALL_inclsex_PCA-plink_PC1vsPC2_ancestry_label.pdf")
qplot(PC1, PC2, colour = ancestry_Grafpop,
      data = p) +
  geom_text(aes(label = IID), size = 1, vjust = -1)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC3vsPC4_ancestry_label.pdf")
qplot(PC3, PC4, colour = ancestry_Grafpop,
      data = p)+
  geom_text(aes(label = IID), size = 1, vjust = -1)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC1vsPC2_ethn_quest_label.pdf")
qplot(PC1, PC2, colour = ancestry_quest,
      data = p) +
  geom_text(aes(label = IID), size = 1, vjust = -1)
dev.off()

pdf("ALL_inclsex_PCA-plink_PC3vsPC4_ethn_quest_label.pdf")
qplot(PC3, PC4, colour = ancestry_quest,
      data = p) +
  geom_text(aes(label = IID), size = 1, vjust = -1)
dev.off()
