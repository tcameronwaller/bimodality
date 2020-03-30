#!/bin/bash

# cp -r /media/tcameronwaller/primary/data/local/work/project/2019_bimodality/dock_template /home/tcameronwaller/dock_template
# mv /home/tcameronwaller/dock_template /home/tcameronwaller/dock


#chmod u+x script.sh

# Execution
# $ bash execution_halyard.sh

path_project="/media/tcameronwaller/primary/data/local/work/project/2020_bimodality/"
subpath_repository="repository/"
path_repository="$path_project$subpath_repository"
subpath_program="repository/bimodality/"
path_program="$path_project$subpath_program"
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

cd $path_program

echo "Now set to call routine and procedures."


##########
# - Preparation Routine

# -- Preparation of dock_template

# --- Preliminary run of selection procedure determines persons of interest.
# --- Calculate genotype principal components across only these persons.
# --- Execute on halyard
#scp /home/tcameronwaller/dock/selection/tight/heritability/simple/families_persons.tsv tcwaller@grenache.ucsd.edu:/cellar/users/tcwaller/Data/dock/access_private/families_persons.tsv
# --- execute on nrnb
#bash access_private.sh # run on 20 February 2020
#bash heritability_initial.sh # run on 20 February 2020
#cd /cellar/users/tcwaller/Data/dock/
#tar -czvf access_private.tar.gz access_private/
# --- execute on halyard
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/
#cp -r /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/dock_template_2020-02-06 dock_template_2020-02-20
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/dock_template_2020-02-20/
#rm access_private.tar.gz
#rm -r access_private
#scp tcwaller@grenache.ucsd.edu:/cellar/users/tcwaller/Data/dock/access_private.tar.gz access_private.tar.gz
#tar -xzvf access_private.tar.gz
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/
#cp -r dock_template_2020-02-20 /home/tcameronwaller/dock_template_2020-02-20
#cd /home/tcameronwaller/
#mv /home/tcameronwaller/dock_template_2020-02-20/ /home/tcameronwaller/dock/

# -- essential procedures
#python3 interface.py main --dock $path_dock --access # run on 21 February 2020
#python3 interface.py main --dock $path_dock --assembly # run on 24 March 2020
#python3 interface.py main --dock $path_dock --selection # run on 26 March 2020

##########
# - Batch Routine

#python3 interface.py main --dock $path_dock --split
#python3 interface.py main --dock $path_dock --distribution --local # run on 21 February 2020
#python3 interface.py main --dock $path_dock --distribution --remote --gene "ENSG00000186092"
#nohup python3 interface.py main --dock $path_dock --distribution --local > $path_dock/report.out &
#python3 interface.py main --dock $path_dock --candidacy

#python3 interface.py main --dock $path_dock --shuffle --count 10
# Current implementation requires approximately 3 Gigabytes of memory for
# 10,000 shuffles on each computational process or node.
# 100,000 shuffles require approximately 30 Gigabytes of memory on each node.
#python3 interface.py main --dock $path_dock --permutation --local
#nohup python3 interface.py main --dock $path_dock --permutation --local > $path_dock/report.out &
#python3 interface.py main --dock $path_dock --permutation --remote --gene "ENSG00000186092"
#python3 interface.py main --dock $path_dock --probability

##########
# Analysis Routine

#bash heritability_genes.sh # run on 4 March 2020
#python3 interface.py main --dock $path_dock --heritability

#python3 interface.py main --dock $path_dock --category # <-- (2019-11-24) this procedure is obsolete for now...
#python3 interface.py main --dock $path_dock --prediction

#python3 interface.py main --dock $path_dock --integration

#python3 interface.py main --dock $path_dock --structure

python3 interface.py main --dock $path_dock --plot


##########
# Nonessential procedures for exploration and testing... several obsolete

#python3 interface.py main --dock $path_dock --measurement
#python3 interface.py main --dock $path_dock --sample
#python3 interface.py main --dock $path_dock --tissue
#python3 interface.py main --dock $path_dock --metric
#python3 interface.py main --dock $path_dock --test

##########
# TODO: probably obsolete scrap now...
#python3 interface.py main --dock $path_dock --combination
#python3 interface.py main --dock $path_dock --analysis # <--- most of this stuff is now in "rank" procedure
#python3 interface.py main --dock $path_dock --function
#python3 interface.py main --dock $path_dock --expecto
#Rscript expression.R $path_dock
