#####################################################
# script to prepare db for QC of maternal GWAS BISC
# 20230506
#####################################################

#####################################################
# libraries

library(readstata13)
library(xlsx)
library(openxlsx)
library(sjmisc)
library(tidyverse)

#####################################################
# IDs

# FamilyID: family ID (in plink FID)
# SubjectID: invidudal ID (in plink IID)
# LabID: ID used in the isglobal.lab - one ID for each DNA - abbreviated format (_bis)
# GWASID: ID used in the cegen.lab - one ID for each genotyped sample (_dup)

#####################################################
# 1052 samples - 22 batches 
# between pilot + main mothers

#####################################################
# read file with info on gwas

### read pilot

setwd("//fs01.isglobal.lan/hpc_bisc_omics/data_raw/gwas/mother_child/GSA_pilot_20220803/IDAT_pilot_20220628")

gp<-read.csv("sampleSheet_pilot.csv", skip=8)
dim(gp)#48
head(gp)

# add round
gp$Round<-"pilot"

### read main mothers

setwd("//fs01.isglobal.lan/hpc_bisc_omics/data_raw/gwas/mother_child/GSA_mothers_20230227/sampleSheet_MBustamante_987m_GSA")

gm<-read.csv("sampleSheet_MBustamante_987m_GSA_concatenados.csv", skip=8)
dim(gm)#1008
head(gm)

# add round
gm$Round<-"main_mother"

### bind rounds

gpm<-rbind(gp, gm)
names(gpm)
gpmf<-gpm[,c(1,2,3,4,5,8,9,13)]
dim(gpmf)#1056 8

### identify controls

gpmf$control<-"no"
gpmf$control[substr(gpmf$Sample_ID, 1, 2) == 'NA']<-"yes"
table(gpmf$control)
#no  yes 
#1035   21

### check duplicates

table(duplicated(gpmf$Sample_ID))
#FALSE  TRUE 
#1055     1
which(duplicated(gpmf$Sample_ID))
gpmf[1056,]#it is a HapMap present in pilot + main_mother round
gpmf[gpmf$Sample_ID=="NA11920",]
#Sample_ID SentrixBarcode_A SentrixPosition_A Sample_Plate Sample_Well Sample_Name Replicate       Round control  Lab_ID
#648    NA11920     207207180019            R12C02       231578         H03     NA11920           main_mother     yes NA11920
#1056   NA11920     207207180067            R12C02       231582         H06     NA11920           main_mother     yes NA11920
gpmf[gpmf$Sample_ID=="NA11920_dup_dup",]
#Sample_ID SentrixBarcode_A SentrixPosition_A Sample_Plate Sample_Well Sample_Name Replicate       Round control    Lab_ID
#696 NA11920_dup_dup     207207180034            R12C02       231578         H09 NA11920_dup   NA11920 main_mother     yes NA11920_d

gpmf$Sample_ID[1056]<-"NA11920_dup"
gpmf$Replicate[1056]<-"NA11920"

### identify replicates

table(gpmf$Replicate)
# 11 HapMaps
# 4 BISC

### create LabID (eliminate "dup")

gpmf$LabID<-str_remove(gpmf$Sample_ID, "_dup")
table(gpmf$LabID%in%gpmf$Sample_ID)
#FALSE  TRUE 
#9   1047
table(duplicated(gpmf$LabID))
#FALSE  TRUE 
#1051     5
which(duplicated(gpmf$LabID))
gpmf[c(118,140,558,594,1056),]
#Sample_ID
#BISC_0674_dup     
#BISC_0675_dup    
#BISC_0681_dup    
#BISC_0682_dup   
#NA11920_dup

### create GWASID

gpmf$GWASID<-gpmf$Sample_ID
table(duplicated(gpmf$GWASID))#0


#####################################################
# read file with correspondence BiSC ID to GWAS ID (stock plates - pilot also plated)
# first mother LabID is 0673

### read info on maternal DNAs

setwd("//fs01/bisc/GWAS/DNA_plates/BISC_XX_11_Preg_DNA_plates/StockPlates")

id<-read.xlsx("BISC_XX_11_Preg_D01_Plates_20230313.xlsx", sheet="list")
dim(id)#1022   14
head(id)
names(id)

id<-id[,c(1,2,4,9,13)]
head(id)

### merge

m<-merge(gpmf, id, by.x="LabID", by.y="LabID", all.x=T)
dim(m)#1056

### if missing in SubjectID, then no maternal DNAs

table(is.na(m$SubjectID))
#FALSE  TRUE 
#1009    47
# there are 21 controls, the rest  might be cord blood samples of the pilot study

### create variable to identify maternal DNAs

m$BISC_mother<-"yes"
m$BISC_mother[is.na(m$SubjectID)]<-"no"
table(m$BISC_mother)
#no  yes 
#47 1009

table(duplicated(m$LabID))
#FALSE  TRUE 
#1051     5

#####################################################
# read file with info on main BiSC variables

### read

setwd("//fs01.isglobal.lan/hpc_bisc_data/analyses/info_sample_selection")

fu<-read.xlsx("BISC_followup_codebook_20230330.xlsx", sheet="list")
dim(fu)
head(fu)

### merge

m1<-merge(m, fu, by.x="SubjectID", by.y="id_mother", all.x=T)
dim(m1)#1056
names(m1)
m1$FamilyID<-substr(m1$SubjectID,1,6)

### add twins

which(duplicated(m1$FamilyID)&m1$BISC_mother=="yes"&!duplicated(m1$SubjectID))
m1[c(27, 137, 619, 650),]
#SubjectID     
#27   10002812       
#137  10013912 
#619  12001813 
#650  12004912

m1$twin<-rep(NA,nrow(m1))
m1$twin[m1$FamilyID=="100028"]<-"yes"
m1$twin[m1$FamilyID=="100139"]<-"yes"
m1$twin[m1$FamilyID=="120018"]<-"yes"
m1$twin[m1$FamilyID=="120049"]<-"yes"
table(m1$twin)
#yes 
#8

### clean variables

names(m1)
m2<-m1[,c(55,1,2,12,16,11,9,10,56,4:7,14,17:54)]



### descriptives for all

table(m2$abortion)
table(m2$stillbirth)
table(m2$lost_follow_up)

### descriptive only in mothers

m3<-m2[m2$BISC_mother=="yes",]
dim(m3)#1009 

table(m3$abortion) #0
table(m3$stillbirth) #1 (we exclude it)
table(m3$lost_follow_up) #98 (we keep them)
table(m3$ci_gral_mare) #2 (we keep them)
table(m3$ci_genet_mare) #33 (we exclude them)
table(m3$twin)#8 (we exclue them at the end)
table(duplicated(m3$SubjectID))#4 (we exclude them at the end)


#####################################################
# save file

setwd("//fs01.isglobal.lan/hpc_bisc_omics/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/db")

write.table(m2, "BISC_db_gwas_20230506.txt", quote=F, sep="\t", row.names=F)
# manually  modified to excel file

####################################################

#####################################################
# Filter
# 20230510
#####################################################

### Read file

setwd("//fs01.isglobal.lan/hpc_bisc_omics/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/db")

data <- read.table("BISC_db_gwas_20230506.txt", header=TRUE, sep="\t")
dim(data) #1056   52

### Filter

dim(subset(data, BISC_mother != "yes" & control != "yes"))[1] #26
dim(subset(data, stillbirth != "no" ))[1] #1
dim(subset(data, ci_genet_mare == "0" ))[1] #33


data_f <- subset(data, control == "yes" | (BISC_mother != "no" & stillbirth != "yes" & ci_genet_mare != "0"))
dim(data_f)[1] #996
table(data_f$control) #no 975  yes 21

#####################################################
# save file

setwd("//fs01.isglobal.lan/hpc_bisc_omics/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/db")

write.table(data_f, "BISC_db_gwas_20230510_filter.txt", quote=F, sep="\t", row.names=F)
