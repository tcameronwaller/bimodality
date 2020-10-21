#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"

path_qtl_tools="$path_user_cellar/plink"

path_gtex_8="$path_user_nrnb/gtex-8"
path_genotype_vcf="$path_gtex_8/gtex-8_genotype.vcf.gz"

path_dock="$path_user_cellar/Data/dock"
path_selection_trait="$path_dock/selection/tight/trait"
path_covariates="$path_selection_trait/data_persons_variables.tsv"
path_distribution_trait="$path_dock/distribution/collection/trait"
path_phenotypes="$path_distribution_trait/data_signals_genes_persons_trait.tsv"
path_phenotypes_compress="$path_phenotypes.zip"


# TODO: I need a list of persons for the analysis... ie those with valid genotypes...

path_persons="$path_access_private/families_persons.tsv"


# TODO: I do need to set up a destination directory...

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

# Create index for genotype file.
tabix -p vcf $path_genotype_vcf

# Create index for phenotype file.
# might need to sort table...
# tabix requires the input table to be compressed by gzip...
# https://bedtools.readthedocs.io/en/latest/content/tools/sort.html
sort -k 1,1 -k 2,2n $path_phenotypes | bgzip > $path_phenotypes_compress
tabix -p bed $path_phenotypes_compress

#



# TODO: list of persons (samples) to include
# TODO: list of genes (phenotypes) to include
# TODO: use --include-samples followed by a list of person identifiers with genotypes
# TODO: use --include-phenotypes followed by a list of multimodal genes






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
$path_gcta --bfile $path_genotype_bed_bim_fam --keep $path_persons --autosome --maf 0.01 --make-grm --out $path_relation_gcta_bed_bim_fam/autosome_common --threads 10
$path_gcta --pfile $path_genotype_pgen_pvar_psam --keep $path_persons --autosome --maf 0.01 --make-grm --out $path_relation_gcta_pgen_pvar_psam/autosome_common --threads 10
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons --autosome --maf 0.01 --make-rel --out $path_relation_plink/autosome_common --threads 10

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
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons --autosome --maf 0.01 --pca 25 --out $path_relation_plink/components
