#!/bin/bash

#SBATCH --job-name=permutation
#SBATCH --output=/cellar/users/tcwaller/Data/dock/out.txt
#SBATCH --error=/cellar/users/tcwaller/Data/dock/error.txt
#SBATCH --mem=10G
#SBATCH --array=100-200%1
#SBATCH --time=5-00:00:00 # days-hours:minutes:seconds

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH
path_project="/cellar/users/tcwaller"
subpath_repository="repository/bimodality-master/bimodality"
path_repository="$path_project/$subpath_repository"
subpath_program="repository/bimodality-master/bimodality"
path_program="$path_project/$subpath_program"
subpath_dock="Data/dock"
path_dock="$path_project/$subpath_dock"

# Define iteration variables.
readarray -t genes < "$path_dock/split/genes.txt"
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
