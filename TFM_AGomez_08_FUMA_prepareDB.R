setwd("Y:/analyses/BiSC_21/006_PlamiQTLs_AG/results/mat_miQTL_TensorQTL") #to modify

library(data.table)
res<-fread("bisc_m_TensorQTL_pvalue1E-05_formated.csv")
dim(res)

#### Headers ####
# Column names are automatically detected based on the following headers (case insensitive).
# 
# SNP | snpid | markername | rsID: rsID
# CHR | chromosome | chrom: chromosome
# BP | pos | position: genomic position (hg19)
# A1 | effect_allele | allele1 | alleleB: affected allele
# A2 | non_effect_allele | allele2 | alleleA: another allele
# P | pvalue | p-value | p_value | pval: P-value (Mandatory)
# OR: Odds Ratio
# Beta | be: Beta
# SE: Standard error
# If your input file has alternative names, these can be entered in the respective input boxes when specifying the input file. Note that any columns with the name listed above but with different element need to be avoided. For example, when the column name is "SNP" but the actual element is an id such as "chr:position" rather than rsID will cause an error.
# Extra columns will be ignored.
# Rows that start with "#" will be ignored.
# Column "N" is described in the Parameters section.
# Be careful with the alleles header in which A1 is defined as effect allele by default. Please specify both effect and non-effect allele column to avoid mislabeling.
# If wrong labels are provided for alleles, it does not affect any annotation and prioritization results. It does however affect eQTLs results (alignment of risk increasing allele of GWAS and tested allele of eQTLs). Be aware of that when you interpret results.

#### Delimiter ####
# Delimiter can be any of white space including single space, multiple space and tab. Because of this, each element including column names must not include any space.

names(res)
# [1] "rsID"       "CHR"        "BP"         "miRNA"      "chr_miRNA"  "pos1_miRNA" "pos2_miRNA" "A1"         "A2"         "pval"       "Beta"      
# [12] "SE"         "af_A1"  

#res1<-res[,c(2,3,1,8,9,11,12,10)]
#names(res1)<-c("CHR", "BP", "SNP", "REF", "ALT", "BETA", "SE", "P")

table(res1$CHR_SNP)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  25 
# 95  81 113 114 105 260  94  80  54  57  55  81  58  50  24  46  28  60  35  24   8  28 100   4

#chr:pos_ref/alt
res$marker <- paste(res$chr_snp, ":", res$pos_snp, "_", res$ref_A2, "/", res$alt_A1, sep = "")


res$chr_snp<-as.numeric(res$chr_snp)
res <- res[order(res$chr_snp,res$pos_snp), ]
write.table(res, "bisc_m_TensorQTL_pvalue1E-05_Rformated.txt", quote=F, sep="\t", row.names=F)
