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

#rm -r $path_genotype_ped
#mkdir $path_genotype_ped

rm -r $path_relation
mkdir $path_relation

##########
# Only convert data format.
#$path_bin/plink2 --vcf $path_genotype --out $path_dock/gtex-8_genotype

# Convert data format.
# Possible to use PLINK to filter by person and minimal allelic frequency.
#$path_plink --vcf $path_genotype --keep $path_persons --maf 0.01 --make-pgen --out $path_dock/gtex-8_genotype
# However, only use PLINK to convert and filter in GCTA.
# GCTA requires PLINK 1 format files, .bed, .bim, and .fam.
#$path_plink --vcf $path_genotype_vcf --no-fid --make-bed --out $path_genotype_ped
# I think that the "no-fid" flag does not change anything when importing a VCF.
#$path_plink --vcf $path_genotype_vcf --make-bed --out $path_genotype_ped/genotype --threads 20

# Calculate genetic relationship matrix (GRM).
# For cis heritability, GRM must be specific to each gene's chromosome.
# Filter by persons and minimal allelic frequence (MAF).
# GCTA's format requirement for list of persons is text with tab delimiters.
#$path_gcta --bfile $path_dock/gtex-8_genotype --autosome --maf 0.01 --make-grm --out $path_dock/gtex-8_grm_autosomes --threads 10
#$path_gcta --bfile $path_genotype_ped --keep $path_persons --chr 6 --maf 0.01 --make-grm --out $path_relation/chromosome_6 --threads 10


# TODO: create cis and trans GRMs for each chromosome
# chromosome 6 --> 6_cis (chrom 6 only), 6_trans (all autosomes other than 6)

##########
# Generate GRM for all autosomal chromosomes.
#$path_gcta --bfile $path_genotype_ped/genotype --keep $path_persons --autosome --maf 0.01 --make-grm --out $path_relation/autosome --threads 20

##########
# Generate GRMs for individual chromosomes.
# These GRMs are directly useful as cis GRMs.
path_cis="$path_heritability/cis"
rm -r $path_cis
mkdir $path_cis
# Iterate across autosomal chromosomes.
for cis in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    $path_gcta --bfile $path_genotype_ped/genotype --keep $path_persons --chr $cis --maf 0.01 --make-grm --out $path_cis/$cis --threads 20
done

##########
# Define combinations of chromosomes for unification of trans GRMs.
path_combinations="$path_heritability/combinations"
rm -r $path_combinations
mkdir $path_combinations
# Iterate across cis autosomal chromosomes.
for cis in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    # Initialize a combination list of trans chromosomes corresponding to the
    # cis chromosome.
    path_combination="$path_combinations/$cis.txt"
    # Iterate across autosomal chromosomes.
    for chromosome in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
    do
        # Determine whether chromosome is cis or trans.
        if [ $chromosome != $cis ]
        then
            # Chromosome is trans.
            # Define complete path to chromosome's GRM.
            path_chromosome="$path_cis/$chromosome"
            # Include chromosome in trans collection.
            echo $path_chromosome >> $path_combination
        fi
    done
done

##########
# Unify GRMs for combinations of trans chromosomes
# multi_grm.txt file should include complete paths to all individual grms to unify...

path_trans="$path_heritability/trans"
rm -r $path_trans
mkdir $path_trans
# Iterate across cis autosomal chromosomes.
for cis in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    # Define path to list of chromosomes for trans combination.
    path_combination="$path_combinations/$cis.txt"
    # Define path to union GRM.
    path_union="$path_trans/$cis"
    $path_gcta --mgrm $path_combination --unify-grm --out $path_union
done

# Analysis
#$path_gcta --grm $path_relation/chromosome_6 --pheno $path_distribution --reml --out $path_result
