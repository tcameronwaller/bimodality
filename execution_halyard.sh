#!/bin/bash

#chmod u+x script.sh

# Execution
# $ bash execution_halyard.sh

path_project="/media/tcameronwaller/primary/data/local/work/project/2019_bimodality/"
subpath_repository="repository/"
path_repository="$path_project$subpath_repository"
subpath_program="repository/bimodality/"
path_program="$path_project$subpath_program"
subpath_dock="dock/"
#path_dock="$path_project$subpath_dock"
path_dock="/home/tcameronwaller/dock/"

# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "Here are your working directories."
echo "project: $path_project"
echo "repository: $path_repository"
echo "dock: $path_dock"
echo "----------"
echo "--------------------------------------------------"

# Echo each command to console.
set -x

cd $path_project
mkdir $path_dock
cd $path_program

echo "Now set to call routine and procedures."


# - Preparation Routine

# -- essential procedures
#python3 interface.py main --dock $path_dock --access
#python3 interface.py main --dock $path_dock --assembly
#python3 interface.py main --dock $path_dock --selection

# -- nonessential, exploratory procedures
#python3 interface.py main --dock $path_dock --measurement
#python3 interface.py main --dock $path_dock --sample
#python3 interface.py main --dock $path_dock --tissue

# - Batch Routine

#python3 interface.py main --dock $path_dock --split

#python3 interface.py main --dock $path_dock --candidacy

python3 interface.py main --dock $path_dock --permutation --count 1000

#python3 interface.py main --dock $path_dock --distribution --local
#python3 interface.py main --dock $path_dock --distribution --remote --gene "ENSG00000186092"
#nohup python3 interface.py main --dock $path_dock --distribution --local > $path_dock/report.out &
#python3 interface.py main --dock $path_dock --organization
#python3 interface.py main --dock $path_dock --restriction
#python3 interface.py main --dock $path_dock --aggregation

#python3 interface.py main --dock $path_dock --collection

# Analysis Routine

#python3 interface.py main --dock $path_dock --combination
#python3 interface.py main --dock $path_dock --regression


#python3 interface.py main --dock $path_dock --analysis
#python3 interface.py main --dock $path_dock --function

#python3 interface.py main --dock $path_dock --metric
#python3 interface.py main --dock $path_dock --plot


#python3 interface.py main --dock $path_dock --test


#Rscript expression.R $path_dock

