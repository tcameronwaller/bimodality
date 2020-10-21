#!/bin/bash

#SBATCH --job-name=permutation
#SBATCH --output=/cellar/users/tcwaller/Data/dock/out.txt
#SBATCH --error=/cellar/users/tcwaller/Data/dock/error.txt
#SBATCH --mem=10G
#SBATCH --array=0-15130%200
#SBATCH --time=5-00:00:00 # days-hours:minutes:seconds

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_user_cellar="/cellar/users/tcwaller"
path_python="$path_user_cellar/anaconda3/bin/python"
path_repository="$path_user_cellar/repository"
path_program="$path_repository/bimodality-master/bimodality"
path_dock="$path_user_cellar/Data/dock"
path_selection="$path_dock/selection/tight"
path_distribution="$path_dock/distribution"
path_permutation="$path_dock/permutation"

# Define iteration variables.
readarray -t genes < "$path_selection/genes_selection.txt"
#count_genes=${#genes[@]}
#indices=$((count_genes-1))
#echo $indices
#echo $count_genes
#12824
#echo ${genes[@]}
#echo $genes[0]

# Iterate on genes.
gene=${genes[$SLURM_ARRAY_TASK_ID]}
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "gene: " $gene

path_permutation_gene="$path_permutation/$gene"

hostname
date

# Determine whether a directory already exists for the current gene.
if [ ! -d $path_permutation_gene ]
then
    # Directory does not already exist for the gene.
    # Execute permutation procedure.
    $path_python $path_program/interface.py main --dock $path_dock --permutation --remote --gene $gene
fi
