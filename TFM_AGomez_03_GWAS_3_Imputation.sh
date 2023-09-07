#!/bin/bash

exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3 RETURN
mkdir -p logs
exec 1>logs/"`date +"%Y%m%d_%H:%M"`"_Imputation.out 2>&1

# ################################################################
# IMPUTE chromosomes
#
# Impute chromosomes
#   chromosomes:  ALL
#   Imputation Software: minimac4
#
#   Usage:
#       sh 01_ImputeData.sh <source_file> <dest_path> <reference_panel> <ncores>
#
#       Example: 
#           sh 01_ImputeData.sh /PROJECTES/SPORTOMICS/data_raw/futfem_20230331_GSA.ped /PROJECTES/SPORTOMICS/data_imputed/ 1000G 16
#
#   Reference Panel: 1000G Phase 3 v5
#       - /PROJECTES/PUBLICDATA/REFERENCES/reference_panels/1000GP_Phase3_v5/*.m3vcf
#
#   Input: 
#       - Source file:
#           A file with data to be imputed file formats allowed: 
#               - .vcf files
#               - .ped files
#               - .pgen
#               - .bgen
#
#       - Destination path:
#           Path to be used to impute data and save the final imputed file
#           IMPORTANT!!: There must be enough space to save the intermediate and final imputation results
#
#       - Reference panel:
#           One of the panels defined in  /PROJECTES/SPORTOMICS/data_qc/gwas/futfem/GSA_20230331/scripts/refpanels.txt,
#           use only the variable name
#               - 1000G
#               - JAPAN
#           NOTE: New reference panels can be added in refpanels.txt
#
#       - ncores: 
#           Number of cores to be used during the process
# 
#   Output:  
#       - Imputed data in vcf format
#
#   Needed modules (HPC):
#        module load bio/PLINK/2.00a3.6-GCC-11.3.0
#        module load bio/BCFtools/1.10.2-GCC-9.3.0
#        module load lang/Java/15.0.1                # Java
#        module load lang/R/4.2.1-foss-2022a
#   Also needs:
#       Samtools 
#           Samtools installed in my conda environment genomics: 
#               source activate "genomics"
#
#  Author: BRGE - Dolors Pelegrí
#  Contact: dolors.pelegri@isglobal.org
# ################################################################

toImpute=$1      # Path with plink files
dest_path=$2     # Path to upload imputed files
ref_p=$3
cores=$4


# Harmonize with reference panel before imputation
function harmonize() {
    fchr=$1
    fftoImpute=$2
    fref_panel=$3
    fref_panel_sufix=$4
    fcores=$5
    
    if ! [[ -f "${fref_panel}/vcf/ALL.chr${fchr}.$fref_panel_sufix.vcf.gz.tbi" ]]; then
        tabix -p vcf $fref_panel/vcf/ALL.chr$fchr.$fref_panel_sufix.vcf.gz # Create tabix index
    fi
    java -jar -Xms16g -Xmx16g -jar /PROJECTES/PUBLICDATA/software/GenotypeHarmonizer-1.4.25/GenotypeHarmonizer.jar  --input $fftoImpute"_"chr$fchr.vcf.gz --output $fftoImpute"_"chr$fchr"_harm" --keep --ref $fref_panel/vcf/ALL.chr$fchr.$fref_panel_sufix

    # ped to vcf
    plink2  --threads $fcores \
            --bfile $fftoImpute"_"chr$fchr"_harm" \
            --snps-only just-acgt \
            --export vcf vcf-dosage=DS id-delim=+ \
            --out $fftoImpute"_"chr$fchr"_harm"
    bgzip -c $fftoImpute"_"chr$fchr"_harm".vcf > $fftoImpute"_"chr$fchr"_harm".vcf.gz
    bcftools index $fftoImpute"_"chr$fchr"_harm".vcf.gz
}


# #. DEBUG.# 
# dest_path="/PROJECTES/SPORTOMICS/data_qc/gwas/futfem/GSA_20230331/ImputedNotPreQC/"
# toImpute="/PROJECTES/SPORTOMICS/data_raw/gwas/futfem/GSA_20230331/PLINK_140423_1220/futfem_20230331_GSA.ped"
# cores=8
# var1="1000G"

minimac4="/PROJECTES/PUBLICDATA/software/minimac4-1.0.2-Linux/bin/minimac4"
ftoImpute=""
basef="$(basename -- $toImpute)"
filename=$(echo $basef | grep -Po '.*(?=\.)')
extension="${basef##*.}"

# Create destination path
 mkdir -p $desti_path

# Set ref_panel
 ref_panel=$(sed -n 's/^'$ref_p'=\(.*\)/\1/p' < "/PROJECTES/PUBLICDATA/REFERENCES/reference_panels/refpanels.txt")
 if [ -z "${ref_panel}" ];
 then
     echo "Unknown Reference panel \"$ref_p\""
     return 0
 fi

# Convert data to VCF
 if [[ -f "${toImpute}" ]];
 then

     echo "Input format: $extension"

     if  [[ "$extension" != "vcf" ]]; 
     then
         if  [[ "$extension" == "ped" ]]; 
         then
             echo "Converting \"$filename.$extension\" to .vcf format"
             # First, we need to convert file to pgen sorting chr
             plink2  --threads $cores \
                     --pedmap $(echo $toImpute | grep -Po '.*(?=\.)') \
                     --make-pgen --sort-vars --chr 1-22, X, Y \
                     --out $dest_path$filename"_tmp_plink"

             plink2  --threads $cores \
                     --pfile $dest_path$filename"_tmp_plink" \
                     --make-pgen --rm-dup error list \
                     --out $dest_path$filename

             rm $dest_path$filenam*_tmp_plink*
             # Now from pgen to vcf
             specialFormat="--pfile $dest_path$filename "

         elif [[ "$extension" == "bed" ]]; then
             echo "Converting \"$filename.$extension\" to .vcf format"
             specialFormat="--bfile $(echo $toImpute | grep -Po '.*(?=\.)') "

         elif [[ "$extension" == "gen" ]]; then
             echo "Converting \"$filename.$extension\" to .vcf format"
             specialFormat="--bgen $(echo $toImpute | grep -Po '.*(?=\.)') "

         else 
             echo "Unknown format \"$extension\""
             exit 1
         fi

         plink_command="plink2 --threads $cores $specialFormat --chr 1-22, X, Y --snps-only just-acgt --export vcf vcf-dosage=DS id-delim=+ --out $dest_path$filename"
         eval " $plink_command"
     fi

     if  [[ "$extension" == "vcf" ]]; then 
         cp $toImpute $dest_path
     fi

     ftoImpute=$dest_path$filename
     #..# ftoImpute=$dest_path$filename.vcf
     bgzip -c $ftoImpute.vcf > $dest_path$filename".vcf.gz"
     bcftools index $ftoImpute".vcf.gz"

    # Remove duplicates
     bcftools norm \
          -d none $dest_path$filename".vcf.gz" \
          --output-type z \
          -o $ftoImpute"_nodupli.vcf.gz"
     bcftools index $ftoImpute"_nodupli.vcf.gz"

 else    
     echo "Error: File \"$toImpute\" does not exists "
     exit 1
 fi


# IMPUTE DATA - Autosomal Chromosomes
#   splits data by chr before imputation
 chrs=($(seq 1 1 22))
 for chr in ${chrs[@]};
 do
     if [[ -f "${ftoImpute}.vcf.gz" ]]; then        
         bcftools view -Oz -r $chr $ftoImpute"_nodupli.vcf.gz" > $ftoImpute"_"chr$chr.vcf.gz
         tabix -p vcf $ftoImpute"_"chr$chr.vcf.gz

         harmonize $chr $ftoImpute $ref_panel phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes $cores

        # Impute data
         eval ${minimac4} --refHaps ${ref_panel}${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
             --haps $ftoImpute"_"chr$chr"_harm".vcf.gz \
             --prefix $ftoImpute"_"chr$chr"_harm_imp" \
             --format GT,DS,GP \
             --ignoreDuplicates --rsid --log --cpus $cores

         bcftools index $ftoImpute"_"chr$chr"_harm_imp.dose.vcf.gz" #Not tested

     else    
         echo "Error: File \"$toImpute\" does not exists "
     fi
 done

# IMPUTE DATA - Chromosome X
 chrs=("X")
 for chr in ${chrs[@]}; do

      if [[ -f "${ftoImpute}.vcf.gz" ]]; then        
         bcftools view -Oz -r $chr $ftoImpute"_nodupli.vcf.gz" > $ftoImpute"_"chr$chr.vcf.gz
         tabix -p vcf $ftoImpute"_"chr$chr.vcf.gz

         harmonize $chr $ftoImpute $ref_panel Non.Pseudo.Auto.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes $cores

        # Impute data
         eval ${minimac4} --refHaps ${ref_panel}${chr}.Non.Pseudo.Auto.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
             --haps $ftoImpute"_"chr$chr"_harm".vcf.gz \
             --prefix $ftoImpute"_"chr$chr"_harm_imp" \
             --format GT,DS,GP \
             --ignoreDuplicates --log --cpus $cores

         bcftools index $ftoImpute"_"chr$chr"_harm_imp.dose.vcf.gz" #Not tested
     else    
         echo "Error: File \"$toImpute\" does not exists "
     fi

 done

bcftools concat ${dest_path}*_harm_imp*.vcf.gz \
         --output-type z \
         -o ${dest_path}${filename}_harm_imp.vcf.gz
bcftools index ${dest_path}${filename}_harm_imp.vcf.gz


# ANNOTATION (3-Steps)
#   0. Remove all ID annotations from .vcf file
#   1. Annotate SNPs using initial file (from Lab)
#   2. Anootate not annotated SNPs with dbSNP database

# Remove all rsID from imputed file previous annotation steps
bcftools annotate -x ID ${dest_path}${filename}_harm_imp.vcf.gz \
         --output-type z \
         -o ${dest_path}${filename}_harm_imp_rsID_stp0.vcf.gz
bcftools index ${dest_path}${filename}_harm_imp_rsID_stp0.vcf.gz

# Annotate rsID after imputation (database: dbSNP) and Index file
bcftools annotate \
         -a /PROJECTES/PUBLICDATA/REFERENCES/GRCh37/dbSNP/GRCh37_b151/All_20180423.vcf.gz \
         -c ID ${dest_path}${filename}_harm_imp_rsID_stp0.vcf.gz \
         --output-type z \
         -o ${dest_path}${filename}_harm_imp_rsID_stp1.vcf.gz
bcftools index ${dest_path}${filename}_harm_imp_rsID_stp1.vcf.gz

# Annotate rsID after imputation with dbSNP database - 
# ID annotated only if ID is present in source file and missing in target file
bcftools annotate \
         -a /PROJECTES/PUBLICDATA/REFERENCES/GRCh37/dbSNP/GRCh37_b151/All_20180423.vcf.gz \
         -c +ID ${dest_path}${filename}_harm_imp_rsID_stp1.vcf.gz \
         --output-type z \
         -o ${dest_path}${filename}_harm_imp_rsID.vcf.gz

#..# bgzip -c ${dest_path}${filename}_imp_rsID.vcf > ${dest_path}${filename}_imp_rsID.vcf.gz
bcftools index ${dest_path}${filename}_harm_imp_rsID.vcf.gz
