#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"

path_gcta="$path_user_cellar/gcta_1.92.3beta3/gcta64"
path_plink="$path_user_cellar/plink2"
path_dock="$path_user_nrnb/dock"

path_persons="$path_user_cellar/Data/heritability/families_persons.tsv"
path_genotype_vcf="$path_user_nrnb/gtex_genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"
path_genotype_ped="$path_dock/gtex-8_genotype"
path_distribution="$path_user_cellar/Data/heritability/families_persons_signals.tsv"
path_relation="$path_dock/relation"
path_result="$path_dock/result"


# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "The script is for heritability analysis."
echo "----------"
echo "--------------------------------------------------"

# Echo each command to console.
set -x

rm -r $path_dock
mkdir $path_dock
rm -r $path_relation
mkdir $path_relation

# Only convert data format.
#$path_bin/plink2 --vcf $path_genotype --out $path_dock/gtex-8_genotype

# Convert data format.
# Possible to use PLINK to filter by person and minimal allelic frequency.
#$path_plink --vcf $path_genotype --keep $path_persons --maf 0.01 --make-pgen --out $path_dock/gtex-8_genotype
# However, only use PLINK to convert and filter in GCTA.
# GCTA requires PLINK 1 format files, .bed, .bim, and .fam.
#$path_plink --vcf $path_genotype_vcf --no-fid --make-bed --out $path_genotype_ped
# I think that the "no-fid" flag does not change anything when importing a VCF.
$path_plink --vcf $path_genotype_vcf --make-bed --out $path_genotype_ped

# Calculate genetic relationship matrix (GRM).
# For cis heritability, GRM must be specific to each gene's chromosome.
# Filter by persons and minimal allelic frequence (MAF).
# GCTA's format requirement for list of persons is text with tab delimiters.
#$path_gcta --bfile $path_dock/gtex-8_genotype --autosome --maf 0.01 --make-grm --out $path_dock/gtex-8_grm_autosomes --threads 10
#$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 6 --maf 0.01 --make-grm --out $path_relation/chromosome_6 --threads 10

$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 1 --maf 0.01 --make-grm --out $path_relation/chromosome_1 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 2 --maf 0.01 --make-grm --out $path_relation/chromosome_2 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 3 --maf 0.01 --make-grm --out $path_relation/chromosome_3 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 4 --maf 0.01 --make-grm --out $path_relation/chromosome_4 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 5 --maf 0.01 --make-grm --out $path_relation/chromosome_5 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 6 --maf 0.01 --make-grm --out $path_relation/chromosome_6 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 7 --maf 0.01 --make-grm --out $path_relation/chromosome_7 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 8 --maf 0.01 --make-grm --out $path_relation/chromosome_8 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 9 --maf 0.01 --make-grm --out $path_relation/chromosome_9 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 10 --maf 0.01 --make-grm --out $path_relation/chromosome_10 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 11 --maf 0.01 --make-grm --out $path_relation/chromosome_11 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 12 --maf 0.01 --make-grm --out $path_relation/chromosome_12 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 13 --maf 0.01 --make-grm --out $path_relation/chromosome_13 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 14 --maf 0.01 --make-grm --out $path_relation/chromosome_14 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 15 --maf 0.01 --make-grm --out $path_relation/chromosome_15 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 16 --maf 0.01 --make-grm --out $path_relation/chromosome_16 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 17 --maf 0.01 --make-grm --out $path_relation/chromosome_17 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 18 --maf 0.01 --make-grm --out $path_relation/chromosome_18 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 19 --maf 0.01 --make-grm --out $path_relation/chromosome_19 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 20 --maf 0.01 --make-grm --out $path_relation/chromosome_20 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 21 --maf 0.01 --make-grm --out $path_relation/chromosome_21 --threads 20
$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 22 --maf 0.01 --make-grm --out $path_relation/chromosome_22 --threads 20

$path_gcta --bfile $path_genotype_ped --keep $path_persons --autosome --maf 0.01 --make-grm --out $path_relation/autosome --threads 20


# Analysis
#$path_gcta --grm $path_relation/chromosome_6 --pheno $path_distribution --reml --out $path_result
