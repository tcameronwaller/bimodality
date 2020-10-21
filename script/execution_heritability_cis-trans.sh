#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"
path_gcta="$path_user_cellar/gcta_1.92.3beta3/gcta64"
path_plink="$path_user_cellar/plink2"

path_genotype="$path_user_nrnb/gtex_genotype"
path_genotype_vcf="$path_genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"
path_genotype_ped="$path_genotype/gtex-8_genotype"

path_dock="$path_user_cellar/Data/dock"
path_selection="$path_dock/selection"
path_persons="$path_selection/families.tsv"
path_heritability="$path_dock/heritability"
path_relation="$path_heritability/relation"

rm -r $path_heritability
mkdir $path_heritability
rm -r $path_relation
mkdir $path_relation

# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "The script is for heritability analysis."
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
#-->   $path_plink --vcf $path_genotype_vcf --make-bed --out $path_genotype_ped --threads 10

##########
# Generate GRM for all autosomal chromosomes.
#$path_gcta --bfile $path_genotype_ped/genotype --keep $path_persons --autosome --maf 0.01 --make-grm --out $path_relation/autosome --threads 20

# TODO: 2019-10-04 I decided to remove the --threads 10 from the end of the --make-grm and --unify-grm lines...
# TODO: 2019-10-04 The script gets killed while trying to read in the .bim file from PLINK

##########
# Generate GRMs for individual chromosomes.
# These GRMs are directly useful as cis GRMs.
path_cis="$path_relation/cis"
rm -r $path_cis
mkdir $path_cis
# Iterate across autosomal chromosomes.
for cis in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    # For cis heritability, GRM must be specific to each gene's chromosome.
    # Filter by persons and minimal allelic frequence (MAF).
    # GCTA's format requirement for list of persons is text with tab delimiters.
    $path_gcta --bfile $path_genotype_ped --keep $path_persons --chr $cis --maf 0.01 --make-grm --out $path_cis/$cis
done

##########
# Define combinations of chromosomes for unification of trans GRMs.
path_combinations="$path_relation/combinations"
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

# TODO: in previous trials, this next part wrote .log files to the correct directory but dumped all other files with strange names
# in a different directory.
# TODO: is there an error in --unify-grm?

# TODO: 2019-10-03
# TODO: GCTA reads the --mgrm lists properly
# TODO: GCTA can't seem to actually read the files that the list specifies

# TODO: the cis GRMs aren't writing properly to the .../cis/ directory... only the .log file

path_trans="$path_relation/trans"
rm -r $path_trans
mkdir $path_trans
# Iterate across cis autosomal chromosomes.
for cis in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    # Define path to list of chromosomes for trans combination.
    path_combination="$path_combinations/$cis.txt"
    # Define path to union GRM.
    $path_gcta --mgrm $path_combination --unify-grm --out $path_trans/$cis
done
