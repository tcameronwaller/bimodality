#!/bin/bash

# Specify arguments for qsub command.

# Specify shell.
#$ -S /bin/bash
# Specify name of job.
#$ -N shuffles
# Transfer system variables, such as PATH.
#$ -V
# Specify use of current working directory.
#$ -cwd
# Specify priority 0-15.
#$ -p -7
# Specify count of concurrent processes.
#$ -tc 250
# Specify long process.
#$ -l long
# Specify memory per process.
# Shuffle procedure with 1000 shuffles requires 1.2 Gigabytes maximum.
# Memory requirement might be greater with more shuffles.
#$ -l h_vmem=3G
# Specify maximal job duration, run time.
#$ -l h_rt=<48:00:00>
#$ -l s_rt=<48:00:00>
# Specify destinations for standard output and error.
#$ -o /cellar/users/tcwaller/out
#$ -e /cellar/users/tcwaller/error

# Specify index for iterative arguments for each processes.
# -t 1-$count_genes
#$ -t 1-18804

# Organize paths.

export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

path_project="/cellar/users/tcwaller"
subpath_repository="repository/bimodality-master/bimodality"
path_repository="$path_project/$subpath_repository"
subpath_program="repository/bimodality-master/bimodality"
path_program="$path_project/$subpath_program"
subpath_dock="Data/dock"
path_dock="$path_project/$subpath_dock"

#cd $path_project
#mkdir $path_dock

# Define iteration variables.
readarray genes < $path_dock/split/genes.txt
#echo ${genes[@]}
#echo $genes[0]
#echo $genes[1]
#echo $genes[2]

count_genes=${#genes[@]}

gene=${genes[$SGE_TASK_ID-1]}

# Execute program.
# Specify complete path to python installation.

hostname
date

$path_project/anaconda3/bin/python $path_program/interface.py main --dock $path_dock --pipe --remote --gene $gene

date
