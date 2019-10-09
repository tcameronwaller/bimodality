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
path_genotype_ped="$path_gtex_8/gtex-8_genotype"

path_gcta="$path_user_cellar/gcta_1.92.3beta3/gcta64"
path_plink="$path_user_cellar/plink2"

rm -r $path_relation
mkdir $path_relation

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
$path_plink --vcf $path_genotype_vcf --make-bed --out $path_genotype_ped --threads 10

##########
# Generate GRM for all autosomal chromosomes.
$path_gcta --bfile $path_genotype_ped --autosome --maf 0.01 --make-grm --out $path_relation/autosome --threads 10
