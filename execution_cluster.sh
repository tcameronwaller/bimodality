#!/bin/bash

#SBATCH --job-name=gene_distribution
#SBATCH -o=$path_dock/out.txt
#SBATCH -e=$path_dock/error.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5000
#SBATCH --array=0-12823%250

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
count_genes=${#genes[@]}
indices=$((count_genes-1))
echo $indices
echo $count_genes
#12824

#echo ${genes[@]}
#echo $genes[0]

gene=${genes[$SLURM_ARRAY_TASK_ID]}
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "gene: " $gene

# Execute program.
# Specify complete path to python installation.

hostname
date

$path_project/anaconda3/bin/python $path_program/interface.py main --dock $path_dock --distribution --remote --gene $gene

date
