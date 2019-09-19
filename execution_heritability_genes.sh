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
path_relation="$path_dock/relation/autosome"

path_heritability="$path_user_cellar/Data/heritability"
path_distribution="data_gene_families_persons_signals.tsv"

path_source_ENSG00000134184="$path_heritability/ENSG00000134184/$path_distribution"
path_source_ENSG00000164308="$path_heritability/ENSG00000164308/$path_distribution"
path_source_ENSG00000183793="$path_heritability/ENSG00000183793/$path_distribution"
path_source_ENSG00000185290="$path_heritability/ENSG00000185290/$path_distribution"
path_source_ENSG00000196436="$path_heritability/ENSG00000196436/$path_distribution"
path_source_ENSG00000197728="$path_heritability/ENSG00000197728/$path_distribution"
path_source_ENSG00000205571="$path_heritability/ENSG00000205571/$path_distribution"
path_source_ENSG00000231925="$path_heritability/ENSG00000231925/$path_distribution"
path_source_ENSG00000274512="$path_heritability/ENSG00000274512/$path_distribution"
path_source_ENSG00000280071="$path_heritability/ENSG00000280071/$path_distribution"
path_source_ENSG00000280670="$path_heritability/ENSG00000280670/$path_distribution"

path_product="$path_dock/product"




# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "The script is for heritability analysis."
echo "----------"
echo "--------------------------------------------------"

# Echo each command to console.
set -x

rm -r $path_product
mkdir $path_product

# Analysis
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000134184 --reml --out $path_product/ENSG00000134184
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000164308 --reml --out $path_product/ENSG00000164308
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000183793 --reml --out $path_product/ENSG00000183793
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000185290 --reml --out $path_product/ENSG00000185290
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000196436 --reml --out $path_product/ENSG00000196436
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000197728 --reml --out $path_product/ENSG00000197728
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000205571 --reml --out $path_product/ENSG00000205571
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000231925 --reml --out $path_product/ENSG00000231925
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000274512 --reml --out $path_product/ENSG00000274512
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000280071 --reml --out $path_product/ENSG00000280071
$path_gcta --grm $path_relation --pheno $path_source_ENSG00000280670 --reml --out $path_product/ENSG00000280670
