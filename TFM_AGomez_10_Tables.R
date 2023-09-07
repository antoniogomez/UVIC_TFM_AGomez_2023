#####################################
#                                   #
# title: 'Tables'                   #
# date: '2023-08-10'                #
#                                   #
#####################################

# Load data
library(readxl)
data <- read_excel("Y:/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/db/BISC_db_gwas_20230506_codebook.xlsx", 
                    sheet = "db_miQTL_m_EUR", col_types = c("text", 
                                                            "text", "text", "numeric", "numeric", 
                                                            "text", "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric", "numeric", 
                                                            "numeric", "numeric"))

#install.packages("gtsummary")
library(gtsummary)
tbl_summary <- data %>% 
  select(hospital_parto, sex, gestational_age_birth, birth_weight, ancestry_child, Trophoblasts, Stromal, Hofbauer, Endothelial, nRBC, Syncytiotrophoblast) %>% # keep only columns of interest
  tbl_summary(     
    #by = InfantSex,   # stratify entire table by outcome
    statistic = list(all_continuous() ~ "{mean} ({sd})", # stats and format for continuous columns
                     all_categorical() ~ "{n} ({p}%)"), # stats and format for categorical columns
    digits = c(Trophoblasts, Stromal, Hofbauer, Endothelial, nRBC, Syncytiotrophoblast) ~ 4, # rounding for continuous columns
    type   = all_categorical() ~ "categorical",                 # force all categorical levels to 
    missing_text = "Missing"                                   # how missing values should display
  )  %>%
  #   add_p(test = list(all_continuous() ~ "t.test", all_categorical() ~ "chisq.test"), group = InfantSex) %>%
  #   #modify_header(bold_labels()) %>% # update the column header
  bold_labels() #%>%
#   #modify_spanning_header(all_stat_cols() ~ "**M_RN.T1**")%>%
#   bold_p()

#tbl_merge <- tbl_merge(list(tbl_summary_sesgo_GENEIDA, tbl_summary_sesgo_1, tbl_summary_sesgo_2, tbl_summary_sesgo_3)) 
# modify_spanning_header(
#   list(
#     all_stat_cols() ~ "**Sesgo Analysis**",
#     starts_with("p.value") ~ "**p-values**"
#   )
# )
tbl_summary



ex <- readRDS("Y:/data_qc/transS/child/0y/placenta_RNAseq_20220803/QC_seqClusterBuster_20220803/denoising/BISC_miRNA_TMM_5cpm_10p_c.lib.rds")

data <- as.data.frame(ex[1])

table1 <- data %>%
  tbl_summary(     
    type = all_continuous() ~ "continuous2",       # indicate that you want to print multiple statistics 
    statistic = all_continuous() ~ c(
      "{mean}",                             # line 1: mean 
      "{geometric.mean}",
      "{min}",
      "{median} ({p50}, {p75}, {p95})",                   # line 2: median and IQR
      "{max}"),                              # line 3: min and max
    digits = all_continuous() ~ 2,                              # rounding for continuous columns
    missing_text = "Missing"                                  # how missing values should display
  )

table1
