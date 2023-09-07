#####################################
#                                   #
# title: 'BISC - GWAS mothers - QC' #
# date: '2023-05-10'                #
#                                   #
#####################################

# The present project shows how QC of BISC GWAS mothers data was performed.


###########
# SERVER  #
###########

## Open interactive session on the server. See 'TFM_AGomez_0_ServerCommands.sh'


################
# BINARIE DATA #
################

# Load module PLINK

module load bio/PLINK/1.9b_6.21-x86_64 #There is other versions: 2.00a3.6-GCC-11.3.0 and 2.00a2.3_x86_6

# Prior to performing any tasks, some other plink files have to be generated in order to perform the analysis. In the following web page can be found the code for generating those files.
# https://www.cog-genomics.org/plink/1.9/data
# https://www.cog-genomics.org/plink/2.0/data

# PILOT MOTHERS

## Creating bed, bim and fam files selected

plink --bfile /PROJECTES/BISC_OMICS/data_raw/gwas/mother_child/GSA_pilot_20220803/PLINKbinno0_pilot_20220628/MBustamante_GSA_48m_piloto_bin --keep /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_pilot_sel_20230516/sample_list.txt --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_pilot_sel_20230516/gwas_pilot

# MAIN MOTHERS

## Creating bed, bim and fam files

plink --file /PROJECTES/BISC_OMICS//data_raw/gwas/mother_child/GSA_mothers_20230227/PLINK_230223_0420/MBustamante_987m_GSA --make-bed --out /PROJECTES/BISC_OMICS//data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_mothers_20230509/gwas_mothers

## Creating bed, bim and fam files selected

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_mothers_20230509/gwas_mothers --keep /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_mothers_sel_20230515/sample_list.txt --make-bed --out ./qc/intermediate_qc/PLINKbinno0_mothers_sel_20230515/gwas_mothers

# ALL MOTHERS Merge file

## Creating bed, bim and fam files merged Pilot + Main Mothers

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_pilot_sel_20230516/gwas_pilot --bmerge /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_mothers_sel_20230515/gwas_mothers --make-bed --out ./qc/intermediate_qc/PLINKbinno0_allmothers_merged_20230516/gwas_allmothers

## Add family ID

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_20230516/gwas_allmothers --update-ids /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/family_list.txt --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/bisc_m

## Add sex 

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/bisc_m --update-sex /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/sex_list.txt 3 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/bisc_m


###################
# QUALITY CONTROL #
###################

# The present protocol details the QC steps for the Genome wide data coming from the BiSC project. We will begin by performing sample quality control, including identification of individuals with outlying missing genotype or heterozygosity rates, identification of individuals with discordant sex information, identification of duplicated or related individuals and identification of individuals of divergent ancestry. We will then perform genotype quality control including calculation of call rates, analysis of minor allele frequency (MAF) and deviation from Hardy-Weinberg equilibrium (HWE).  The quality control analysis is performed with the PLINK v1.90b6.2 (12 June 2018) software (Purcell et al. 2007).

################
# 1. Sample QC #
################

# 1.1. Sample Call Rate

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/bisc_m --missing --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m # This command creates the files bisc_m.imiss (Samples) and .lmiss (SNPs) files.

## Filter SNPs
  
plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/bisc_m --geno 0.03 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97 # After discussing about SNPs to discard, we decided to discard those with a F_MISS>0.03 (97%)

head /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m.imiss
  
## Filter Samples
  
### Make list Samples with SNP call rate <97%

awk '$6 > 0.03 {print $1, $2}' /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc97.imiss > /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/exclude_sample_less97callrate.txt # There are 4 samples that <97% sample call rate and they should be eliminated before relatedness because otherwise they come out as related to all of them and it is an error

### Remove with PLINK

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97 --remove /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/exclude_sample_less97callrate.txt --make-bed --out ./qc/intermediate_qc/bisc_m_snpqc97_sampleqc97


# 1.2. Check sex

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97_sampleqc97 --check-sex --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc # We need .map and .ped files

## Visualize chart: See 'TFM_AGomez_Graphics.R'

# 1.3. Heterozygosity

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97_sampleqc97 --het --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc

head /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc.het # ± 4 standard deviations from the mean are excluded

## Visualize chart: See 'TFM_AGomez_Graphics.R'


##############
## 2. SNP QC #
##############

# 1. Allele frequency

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97_sampleqc97 --freq --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/snp_qc/bisc_m_snpqc

# 2. Hardy-Weinberg Equilibrium (HWE)

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97_sampleqc97 --hardy --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/snp_qc/bisc_m_snpqc

## Filter HWE

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97_sampleqc97 --hwe 0.000001 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc


## Filter Minor Allele Frequency (MAF)

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc --maf 0.05 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf5 # MAF>5%

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc --maf 0.01 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf1 # MAF>1%

###################
# 3. RELATEDNESS  #
###################

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc97_sampleqc97 --genome --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc

# Duplicates

awk 'NR==1 || $10 >= 0.9' /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc.genome > /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/duplicates.txt

# Another relatednees

awk 'NR==1 || ($10 >= 0.4 && $10 < 0.9)' /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/bisc_m_snpqc.genome > /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/siblings.txt


# Filter samples
plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_sex_family_20230516/bisc_m --keep /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_FinalSelection_20230529/Sample_Selection.txt --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_FinalSelection_20230529/bisc_m_snpqc


###############
# 4. ANCESTRY #
###############

# Filter HWE <1E-6

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/PLINKbinno0_allmothers_merged_FinalSelection_20230529/bisc_m_snpqc --hwe 0.000001 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc


# Filter MAF

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc --maf 0.05 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf5

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc --maf 0.01 --make-bed --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf1


# Grafpop 1.0

## Intall Grafpop

### mkdir grafpop #Create the folder
### cd grafpop
### wget https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetZip.cgi?zip_name=GrafPop1.0.tar.gz #download
### tar -zxvf GetZip.cgi?zip_name=GrafPop1.0.tar.gz #Unzip software

## Analysis

./grafpop/grafpop /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf5.bed  /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_v2.txt

## Load modules

module load lang/Perl/5.30.0-GCCcore-8.3.0
module load lib/libgd/2.2.5-GCCcore-7.3.0

## Install CPAN to be able to install Perl modules. Accept all default options

cpan App::cpanminus

## Load shell of CPAN
perl -e shell -MCPAN

## We install the necessary modules to make the grafpop plots
install CGI
install GD
install GD::Text
install GD::Graph

## Close shell of CPAN
exit 

## Execute the code to generate the plot

## Plot
perl ./grafpop/PlotGrafPopResults.pl /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores.txt /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores.png

## SaveSample
perl ./grafpop/SaveSamples.pl /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores.txt /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_savesample.txt


# Grafpop 2.4

## Intall Grafpop

mkdir grafpop2 #Create the folder
cd grafpop2
wget http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetZip.cgi?zip_name=GRAF_files.zip #download
tar -zxvf GetZip.cgi?zip_name=GRAF_files.zip #Unzip software

## Analysis

./grafpop2/graf -plink /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf5  -pop /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_v1.txt

perl ./grafpop2/PlotPopulations.pl /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_v1.txt /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_popv1.txt

# Plots
perl ./grafpop2/PlotPopulations.pl /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_v1.txt /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_gd2.png -spf /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/SbjSuperPop.txt -dot 5 -cutoff 1 -pops 2,3,4,5,6,7,8,9 -gw 1200 -ymin 1.1 -ymax 1.5 -xmin 1.2 -xmax 2

perl ./grafpop2/PlotPopulations.pl /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_v1.txt /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/bisc_m_snpqc_maf5_pop_scores_gd4.png -spf /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/ancestry/SbjSuperPop.txt -dot 5 -cutoff 1 -pops 2,3,4,5,6,7,8,9 -gw 1200 -gd4 1


#####################################
# 5. Principal Componentes Analysis #
#####################################

# Prunning SNPs r2 0.2

plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf1 --indep-pairwise 50 5 0.2 --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf1_inclsex_prunned

# PCA
plink --bfile /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf1 --extract /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/intermediate_qc/bisc_m_snpqc_maf1_inclsex_prunned.prune.in --pca --out /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/qc/sample_qc/PCA/bisc_m_snpqc_maf1_inclsex_PC

# Plot: See 'TFM_AGomez_Graphics.R'
