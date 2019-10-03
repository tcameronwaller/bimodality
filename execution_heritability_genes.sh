#!/bin/bash

#SBATCH --job-name=heritability
#SBATCH --output=/cellar/users/tcwaller/Data/dock/heritability_out.txt
#SBATCH --error=/cellar/users/tcwaller/Data/dock/heritability_error.txt
#SBATCH --mem=5G
#SBATCH --array=0-15445%10
#SBATCH --time=2-00:00:00 # days-hours:minutes:seconds

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_user_cellar="/cellar/users/tcwaller"
path_python="$path_user_cellar/anaconda3/bin/python"
path_gcta="$path_user_cellar/gcta_1.92.3beta3/gcta64"
path_plink="$path_user_cellar/plink2"
path_repository="$path_user_cellar/repository"
path_program="$path_repository/bimodality-master/bimodality"
path_dock="$path_user_cellar/Data/dock"
path_split="$path_dock/split"
path_distribution="$path_dock/distribution"
path_heritability="$path_dock/heritability"
path_relation="$path_heritability/relation"

# Define iteration variables.
readarray -t genes < "$path_split/genes.txt"
#count_genes=${#genes[@]}
#indices=$((count_genes-1))
#echo $indices
#echo $count_genes
#12824

#echo ${genes[@]}
#echo $genes[0]

gene=${genes[$SLURM_ARRAY_TASK_ID]}
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "gene: " $gene

# Execute program.
# Specify complete path to python installation.

# Access information about gene.
path_distribution_gene="$path_distribution/$gene"
path_chromosome="$path_distribution_gene/chromosome.txt"
chromosome=$(cat $path_chromosome)
path_relation_cis="$path_relation/cis/chromosome"
path_relation_trans="$path_relation/trans/chromosome"
path_gene_covariates=

path_phenotype="$path_distribution_gene/data_gene_families_persons_signals.tsv"

path_heritability_gene="$path_heritability/$gene"
rm -r $path_heritability_gene
mkdir $path_heritability_gene


path_gene_distribution="$path_heritability/$gene/$path_distribution"
path_gene_heritability="$path_dock/heritability/$gene"

# TODO: read the gene's chromosome...
# chromosome 6 --> 6_cis (chrom 6 only), 6_trans (all autosomes other than 6)



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

# TODO: read in array of genes from "split" directory...
# TODO: define references to gene's phenotype file...
# TODO: eventually include covariates


# Analysis
#$path_gcta --grm $path_relation/chromosome_6 --pheno $path_distribution --reml --out $path_result


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
