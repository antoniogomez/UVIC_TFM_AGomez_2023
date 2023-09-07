
# To the TensorQTL output file, you would have to paste the info from the bim file: chr, pos, A1, A2 from the SNP. 

trans_df <- read.csv("Y:/analyses/BiSC_21/006_PlamiQTLs_AG/results/mat_miQTL_TensorQTL/bisc_m_TensorQTL_pvalue1E-05.csv")

bim <- read.delim("Y:/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf5_EUR.bim", header=FALSE)

bim <- bim[,-3]

colnames(bim) <- c("CHR","rsID", "BP", "A1", "A2")
colnames(trans_df) <- c("rsID","phenotype_id", "pval", "Beta", "SE", "af_A1")

merged_df <- merge(trans_df, bim, by.x = "rsID", all.x = TRUE)
write.table(merged_df, "bisc_m_TensorQTL-bim_pvalue1E-05.txt", quote=F, sep="\t", row.names=F)
  
