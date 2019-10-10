#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"
path_dock="$path_user_cellar/Data/dock"
path_dock_gtex_8="$path_dock/gtex-8"
path_relation="$path_dock_gtex_8/relation"

path_gtex_8="$path_user_nrnb/gtex-8"
path_genotype_vcf="$path_gtex_8/gtex-8_genotype.vcf.gz"
path_genotype_bed="$path_gtex_8/gtex-8_genotype"

path_plink_1="$path_user_cellar/plink"
path_plink_2="$path_user_cellar/plink2"
path_gcta="$path_user_cellar/gcta_1.92.4beta/gcta64"

#rm -r $path_relation
#mkdir $path_relation

# Temporary test directories...
path_relation_gcta="$path_relation/gcta"
path_relation_plink_1="$path_relation/plink_1"
path_relation_plink_2="$path_relation/plink_2"

rm -r $path_relation_gcta
rm -r $path_relation_plink_1
rm -r $path_relation_plink_2
mkdir $path_relation_gcta
mkdir $path_relation_plink_1
mkdir $path_relation_plink_2

# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "The script is to initialize genetic relationship matrix (GRM) for"
echo "heritability analysis in GCTA."
echo "----------"
echo "--------------------------------------------------"

# Echo each command to console.
set -x

##########
# Convert data format.
# Possible to use PLINK to filter by person and minimal allelic frequency.
#$path_plink --vcf $path_genotype --keep $path_persons --maf 0.01 --make-pgen --out $path_dock/gtex-8_genotype
# However, only use PLINK to convert and filter in GCTA.
# GCTA requires PLINK 1 format files, .bed, .bim, and .fam.
#$path_plink --vcf $path_genotype_vcf --no-fid --make-bed --out $path_genotype_ped
# I think that the "no-fid" flag does not change anything when importing a VCF.
#$path_plink_2 --vcf $path_genotype_vcf --make-bed --out $path_genotype_bed --threads 10


##########
# Generate GRM for all autosomal chromosomes.
#$path_gcta --bfile $path_genotype_bed --autosome --maf 0.01 --make-grm --out $path_relation/autosome_common --threads 10
#$path_gcta --bfile $path_genotype_bed --autosome --make-grm --out $path_relation/autosome_rare-common --threads 10

##########
# Calculate principal components.
# Preserve rare variants in this analysis (do not filter by minor allele
# frequency).
#$path_gcta --grm $path_relation/autosome_rare-common --pca 10 --out $path_relation_gcta/autosome_rare-common
#$path_gcta --grm $path_relation/autosome_common --pca 10 --out $path_relation_gcta/autosome_common

$path_plink_1 --bfile $path_genotype_bed --autosome --maf 0.01 --pca 10 "header" "tabs" --out $path_relation_plink_1/autosome_common
#$path_plink_2 --bfile $path_genotype_bed --autosome --maf 0.01 --pca 10 "header" "tabs" --out $path_relation_plink_2/autosome_common
