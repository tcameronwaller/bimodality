#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"

path_plink_1="$path_user_cellar/plink"
path_plink_2="$path_user_cellar/plink2"
path_gcta="$path_user_cellar/gcta_1.93.0beta/gcta64"

path_gtex_8="$path_user_nrnb/gtex-8"
path_genotype_vcf="$path_gtex_8/gtex-8_genotype.vcf.gz"
path_gtex_bed_bim_fam="$path_gtex_8/bed_bim_fam"
path_gtex_pgen_pvar_psam="$path_gtex_8/pgen_pvar_psam"
path_genotype_bed_bim_fam="$path_gtex_bed_bim_fam/gtex-8_genotype"
path_genotype_pgen_pvar_psam="$path_gtex_pgen_pvar_psam/gtex-8_genotype"

path_access_private="$path_user_cellar/Data/dock/access_private"
path_relation_gcta="$path_access_private/relation/gcta"
path_relation_gcta_bed_bim_fam="$path_relation_gcta/bed_bim_fam"
path_relation_gcta_pgen_pvar_psam="$path_relation_gcta/pgen_pvar_psam"
path_relation_plink="$path_access_private/relation/plink"

rm -r $path_gtex_bed_bim_fam
mkdir -p $path_gtex_bed_bim_fam
rm -r $path_gtex_pgen_pvar_psam
mkdir -p $path_gtex_pgen_pvar_psam

rm -r $path_relation_gcta
mkdir -p $path_relation_gcta_bed_bim_fam
mkdir -p $path_relation_gcta_pgen_pvar_psam
rm -r $path_relation_plink
mkdir -p $path_relation_plink

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
# PLINK 1 binary files: .bed, .bim, .fam
$path_plink_2 --vcf $path_genotype_vcf --make-bed --out $path_genotype_bed_bim_fam --threads 10
# PLINK 2 binary files: .pgen, .pvar, .psam
$path_plink_2 --vcf $path_genotype_vcf --make-pgen --out $path_genotype_pgen_pvar_psam --threads 10

##########
# Generate GRM for all autosomal chromosomes.
$path_gcta --bfile $path_genotype_bed_bim_fam --autosome --maf 0.01 --make-grm --out $path_relation_gcta_bed_bim_fam/autosome_common --threads 10
$path_gcta --pfile $path_genotype_pgen_pvar_psam --autosome --maf 0.01 --make-grm --out $path_relation_gcta_pgen_pvar_psam/autosome_common --threads 10
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --autosome --maf 0.01 --make-rel --out $path_relation_plink/autosome_common --threads 10

##########
# Calculate principal components.
# Plink2 does not support principal components yet.
# Plink1 gives memory error when trying to calculate principal components.
#$path_plink_1 --bfile $path_genotype_bed --autosome --maf 0.01 --pca 10 "header" "tabs" --out $path_relation_plink_1/autosome_common
# Use the Genetic Relationship Matrix (GRM) that does not include rare variants
# (minor allele frequence < 0.01).
# Principal components on rare variants produces Eigenvalues near zero and
# missing values across Eigenvectors.
$path_gcta --grm $path_relation_gcta_pgen_pvar_psam/autosome_common --pca 25 --out $path_relation_gcta_pgen_pvar_psam/components
$path_gcta --grm $path_relation_gcta_bed_bim_fam/autosome_common --pca 25 --out $path_relation_gcta_bed_bim_fam/components
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --autosome --maf 0.01 --pca 25 --out $path_relation_plink/components
