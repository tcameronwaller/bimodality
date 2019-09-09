#!/bin/bash

#chmod u+x script.sh

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"
subpath_dock="gtex_genotype/dock"
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

mkdir $path_dock

plink2 --vcf $path_genotype --keep $path_persons --maf 0.01 --make-pgen --out $path_dock/fileset
