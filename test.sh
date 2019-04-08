#!/bin/bash

# Organize paths.

export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

path_project="/media/tcameronwaller/primary/data/local/work/project/2019_bimodality/"
subpath_repository="repository/"
path_repository="$path_project$subpath_repository"
subpath_program="repository/bimodality/"
path_program="$path_project$subpath_program"
subpath_dock="dock/"
path_dock="$path_project$subpath_dock"

# Define iteration variables.
readarray genes < $path_dock/split/genes.txt
#echo ${genes[@]}
#echo $genes[0]
#echo $genes[1]
#echo $genes[2]

count_genes=${#genes[@]}

echo $count_genes

gene=${genes[$SGE_TASK_ID-1]}

# Execute program.
# Specify complete path to python installation.

hostname
date


date
