#!/bin/bash

# cp -r /media/tcameronwaller/primary/data/local/work/project/2019_bimodality/dock_template /home/tcameronwaller/dock_template
# mv /home/tcameronwaller/dock_template /home/tcameronwaller/dock


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


##########
# - Preparation Routine

# -- essential procedures
#python3 interface.py main --dock $path_dock --access
#python3 interface.py main --dock $path_dock --assembly
#python3 interface.py main --dock $path_dock --selection

# -- nonessential, exploratory procedures
#python3 interface.py main --dock $path_dock --measurement
#python3 interface.py main --dock $path_dock --sample
#python3 interface.py main --dock $path_dock --tissue

##########
# - Batch Routine

#python3 interface.py main --dock $path_dock --split

#python3 interface.py main --dock $path_dock --distribution --local
#python3 interface.py main --dock $path_dock --distribution --remote --gene "ENSG00000186092"
#nohup python3 interface.py main --dock $path_dock --distribution --local > $path_dock/report.out &

#python3 interface.py main --dock $path_dock --shuffle --count 10000
#python3 interface.py main --dock $path_dock --permutation --local
#nohup python3 interface.py main --dock $path_dock --permutation --local > $path_dock/report.out &
#python3 interface.py main --dock $path_dock --permutation --remote --gene "ENSG00000186092"

python3 interface.py main --dock $path_dock --collection

# -- procedure for trial
#python3 interface.py main --dock $path_dock --metric

##########
# Analysis Routine

#python3 interface.py main --dock $path_dock --rank
#python3 interface.py main --dock $path_dock --category
#python3 interface.py main --dock $path_dock --plot

##########
# New stuff...
# TODO: I'll probably want a new, "integration" procedure
# TODO: bring together gene ranks from bimodality with genes from heritability?
#python3 interface.py main --dock $path_dock --integration

#python3 interface.py main --dock $path_dock --test

##########
# TODO: probably obsolete scrap now...
#python3 interface.py main --dock $path_dock --combination
#python3 interface.py main --dock $path_dock --analysis # <--- most of this stuff is now in "rank" procedure
#python3 interface.py main --dock $path_dock --function
#python3 interface.py main --dock $path_dock --expecto
#Rscript expression.R $path_dock
