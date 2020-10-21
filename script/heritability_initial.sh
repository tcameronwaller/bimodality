#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"

path_plink_1="$path_user_cellar/plink"
path_plink_2="$path_user_cellar/plink2"
path_gcta="$path_user_cellar/gcta_1.93.2beta/gcta64"

# Source files.

path_gtex="/nrnb/data/controlled/2020_dbGaP_GTEx/75875/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU"
path_genotype_vcf="$path_gtex/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"

path_persons_properties="$path_user_cellar/Data/dock/selection/tight/persons_properties"
path_persons_selection="$path_persons_properties/selection/persons.txt"
path_families_selection="$path_persons_properties/selection/heritability/families_persons.tsv"
path_persons_respiration="$path_persons_properties/respiration/persons.txt"
path_families_respiration="$path_persons_properties/respiration/heritability/families_persons.tsv"
path_persons_ventilation="$path_persons_properties/ventilation/persons.txt"
path_families_ventilation="$path_persons_properties/ventilation/heritability/families_persons.tsv"

# Product files.

path_gtex_bed_bim_fam="$path_user_nrnb/gtex_8/bed_bim_fam"
path_gtex_pgen_pvar_psam="$path_user_nrnb/gtex_8/pgen_pvar_psam"
path_genotype_bed_bim_fam="$path_gtex_bed_bim_fam/gtex_8_genotype"
path_genotype_pgen_pvar_psam="$path_gtex_pgen_pvar_psam/gtex_8_genotype"
rm -r $path_gtex_bed_bim_fam
mkdir -p $path_gtex_bed_bim_fam
rm -r $path_gtex_pgen_pvar_psam
mkdir -p $path_gtex_pgen_pvar_psam

path_access_private="$path_user_cellar/Data/dock/access_private"

path_relation_selection="$path_access_private/relation/selection"
path_selection_gcta_bed_bim_fam="$path_relation_selection/gcta/bed_bim_fam"
path_selection_gcta_pgen_pvar_psam="$path_relation_selection/gcta/pgen_pvar_psam"
path_selection_plink="$path_relation_selection/plink"
rm -r $path_relation_selection
mkdir -p $path_selection_gcta_bed_bim_fam
mkdir -p $path_selection_gcta_pgen_pvar_psam
mkdir -p $path_selection_plink
cp $path_persons_selection "$path_relation_selection/persons.txt"
cp $path_families_selection "$path_relation_selection/families_persons.tsv"

path_relation_respiration="$path_access_private/relation/respiration"
path_respiration_gcta_bed_bim_fam="$path_relation_respiration/gcta/bed_bim_fam"
path_respiration_gcta_pgen_pvar_psam="$path_relation_respiration/gcta/pgen_pvar_psam"
path_respiration_plink="$path_relation_respiration/plink"
rm -r $path_relation_respiration
mkdir -p $path_respiration_gcta_bed_bim_fam
mkdir -p $path_respiration_gcta_pgen_pvar_psam
mkdir -p $path_respiration_plink
cp $path_persons_respiration "$path_relation_respiration/persons.txt"
cp $path_families_respiration "$path_relation_respiration/families_persons.tsv"

path_relation_ventilation="$path_access_private/relation/ventilation"
path_ventilation_gcta_bed_bim_fam="$path_relation_ventilation/gcta/bed_bim_fam"
path_ventilation_gcta_pgen_pvar_psam="$path_relation_ventilation/gcta/pgen_pvar_psam"
path_ventilation_plink="$path_relation_ventilation/plink"
rm -r $path_relation_ventilation
mkdir -p $path_ventilation_gcta_bed_bim_fam
mkdir -p $path_ventilation_gcta_pgen_pvar_psam
mkdir -p $path_ventilation_plink
cp $path_persons_ventilation "$path_relation_ventilation/persons.txt"
cp $path_families_ventilation "$path_relation_ventilation/families_persons.tsv"

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

# GCTA does not seem to work well with .pgen, .pvar, .psam format files.

# Selection cohort.
$path_gcta --bfile $path_genotype_bed_bim_fam --keep $path_families_selection --autosome --maf 0.01 --make-grm --out $path_selection_gcta_bed_bim_fam/autosome_common --threads 10
#$path_gcta --pfile $path_genotype_pgen_pvar_psam --keep $path_families_selection --autosome --maf 0.01 --make-grm --out $path_selection_gcta_pgen_pvar_psam/autosome_common --threads 10
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons_selection --autosome --maf 0.01 --make-rel --out $path_selection_plink/autosome_common --threads 10

# Respiration cohort.
$path_gcta --bfile $path_genotype_bed_bim_fam --keep $path_families_respiration --autosome --maf 0.01 --make-grm --out $path_respiration_gcta_bed_bim_fam/autosome_common --threads 10
#$path_gcta --pfile $path_genotype_pgen_pvar_psam --keep $path_families_respiration --autosome --maf 0.01 --make-grm --out $path_respiration_gcta_pgen_pvar_psam/autosome_common --threads 10
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons_respiration --autosome --maf 0.01 --make-rel --out $path_respiration_plink/autosome_common --threads 10

# Ventilation cohort.
$path_gcta --bfile $path_genotype_bed_bim_fam --keep $path_families_ventilation --autosome --maf 0.01 --make-grm --out $path_ventilation_gcta_bed_bim_fam/autosome_common --threads 10
#$path_gcta --pfile $path_genotype_pgen_pvar_psam --keep $path_families_ventilation --autosome --maf 0.01 --make-grm --out $path_ventilation_gcta_pgen_pvar_psam/autosome_common --threads 10
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons_ventilation --autosome --maf 0.01 --make-rel --out $path_ventilation_plink/autosome_common --threads 10

##########
# Calculate principal components.
# Plink2 does not support principal components yet.
# Plink1 gives memory error when trying to calculate principal components.
#$path_plink_1 --bfile $path_genotype_bed --autosome --maf 0.01 --pca 10 "header" "tabs" --out $path_relation_plink_1/autosome_common
# Use the Genetic Relationship Matrix (GRM) that does not include rare variants
# (minor allele frequency < 0.01).
# Principal components on rare variants produces Eigenvalues near zero and
# missing values across Eigenvectors.

# Selection cohort.
#$path_gcta --grm $path_selection_gcta_pgen_pvar_psam/autosome_common --pca 25 --out $path_selection_gcta_pgen_pvar_psam/components
$path_gcta --grm $path_selection_gcta_bed_bim_fam/autosome_common --pca 25 --out $path_selection_gcta_bed_bim_fam/components
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons_selection --autosome --maf 0.01 --pca 25 --out $path_selection_plink/components

# Respiration cohort.
#$path_gcta --grm $path_respiration_gcta_pgen_pvar_psam/autosome_common --pca 25 --out $path_respiration_gcta_pgen_pvar_psam/components
$path_gcta --grm $path_respiration_gcta_bed_bim_fam/autosome_common --pca 25 --out $path_respiration_gcta_bed_bim_fam/components
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons_respiration --autosome --maf 0.01 --pca 25 --out $path_respiration_plink/components

# Ventilation cohort.
#$path_gcta --grm $path_ventilation_gcta_pgen_pvar_psam/autosome_common --pca 25 --out $path_ventilation_gcta_pgen_pvar_psam/components
$path_gcta --grm $path_ventilation_gcta_bed_bim_fam/autosome_common --pca 25 --out $path_ventilation_gcta_bed_bim_fam/components
$path_plink_2 --pfile $path_genotype_pgen_pvar_psam --keep $path_persons_ventilation --autosome --maf 0.01 --pca 25 --out $path_ventilation_plink/components
