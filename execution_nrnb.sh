#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_bin="/cellar/users/tcwaller/anaconda3/bin"
path_plink="/cellar/users/mpagadal/Programs"
path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"
subpath_dock="dock"
subpath_persons="Data/dock/selection/persons.txt"
subpath_genotype="gtex_genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"
path_dock="$path_user_nrnb/$subpath_dock"
path_persons="$path_user_cellar/$subpath_persons"
path_genotype="$path_user_nrnb/$subpath_genotype"

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

# Only convert data format.
#$path_bin/plink2 --vcf $path_genotype --out $path_dock/gtex-8_genotype

# Apply filters and convert data format.
$path_plink/plink2 --vcf $path_genotype --keep $path_persons --maf 0.01 --make-pgen --out $path_dock/gtex-8_genotype
