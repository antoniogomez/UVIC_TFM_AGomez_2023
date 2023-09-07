#!/bin/bash

#SBATCH --job-name=bisc_m_imputation # Set a Job name
#SBATCH --partition=long # Set a partition or queue to submit to (pot ser short, long o no_limits (però per no_limits cal demanar permís a SRI)
#SBATCH --account=generic # Set a user account to submit as (always generic, unless there is a special request to SRI for the unlimited one)
#SBATCH --mail-type=END,FAIL,BEGIN # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=antonio.gomez@isglobal.org # Where to send mail events
#SBATCH --nodes=1 # Run all processes on a single node (you can change this)
#SBATCH --ntasks=1 # Run single tasks
#SBATCH --cpus-per-task=24 # Use 4 cpus for each task (you can change this)
#SBATCH --mem=128gb # Job memory request (you can change this)
#SBATCH --output=/home/isglobal.lan/agomez/logs/bisc_m_imputation1.log # Standard output and error log

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load R modules to execute script:
source ~/.bashrc
module load bio/PLINK/2.00a3.6-GCC-11.3.0   # plink2
module load bio/BCFtools/1.10.2-GCC-9.3.0   # bcftools
module load lang/Java/15.0.1                # Java
module load lang/R/4.2.1-foss-2022a         # R

# Load conda environment
source activate genomics # Conda environment with samtools (tabix)

# And finally run the job
cd  /PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/scripts/

source_file=/PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/imp_1000g/input/bisc_m_snpqc_maf5_nohapmap.bed
dest_path=/PROJECTES/BISC_OMICS/data_qc/gwas/mother_child/gsa_qc_mothers_20230502/imp_1000g/output/
ref_panel="1000G"
ncores=24

source TFM_AGomez_GWAS_3_Imputation.sh $source_file $dest_path $ref_panel $ncores
