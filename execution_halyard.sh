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

# --- Preliminary run of selection procedure determines persons of interest in
# --- each sub cohort (selection, respiration, ventilation).
# --- Calculate genotype principal components across only these persons.
# --- Execute on halyard
# ---- Execute access, assembly, and selection procedures
# ---- heritability_initial and heritability_genes procedures need information
# ---- from selection procedure
#tar -czvf dock_lite.tar.gz dock_lite/
#scp /home/tcameronwaller/dock_lite.tar.gz tcwaller@grenache.ucsd.edu:/cellar/users/tcwaller/Data/dock_lite.tar.gz

# --- execute on nrnb
#tar -xzvf /cellar/users/tcwaller/Data/dock_lite.tar.gz
#mv /cellar/users/tcwaller/Data/dock_lite/ /cellar/users/tcwaller/Data/dock/
#bash access_private.sh # run on 19 June 2020
#bash heritability_initial.sh # run on 19 June 2020
#cd /cellar/users/tcwaller/Data/dock/
#du -h selection/
#rm -r selection/
#tar -czvf access_private.tar.gz access_private/

# --- execute on halyard
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/
#cp -r /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/dock_template_2020-06-18 dock_template_2020-06-19
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/dock_template_2020-06-19/
#rm access_private.tar.gz
#rm -r access_private
#scp tcwaller@grenache.ucsd.edu:/cellar/users/tcwaller/Data/dock/access_private.tar.gz access_private.tar.gz
#tar -xzvf access_private.tar.gz
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/
#rm -r /home/tcameronwaller/dock/
#cp -r dock_template_2020-06-19 /home/tcameronwaller/dock_template_2020-06-19
#cd /home/tcameronwaller/
#mv /home/tcameronwaller/dock_template_2020-06-19/ /home/tcameronwaller/dock/

##########
# Organization routine

# --- execute on halyard
#cd /media/tcameronwaller/primary/data/local/work/project/2020_bimodality/repository/

# -- essential procedures
#python3 interface.py main --dock $path_dock --access # run on 19 June 2020
#python3 interface.py main --dock $path_dock --assembly # run on 1 August 2020
#python3 interface.py main --dock $path_dock --selection # run on 1 August 2020 for 15,742 genes (50% threshold)
#python3 interface.py main --dock $path_dock --collection # run on 19 August 2020 for 15,742 genes (50% threshold)

##########
# - Batch Routine

#python3 interface.py main --dock $path_dock --split # run on 15 July 2020 for 17,650 genes
#python3 interface.py main --dock $path_dock --distribution --local # run on 15 July 2020 for 17,650 genes (10% signal threshold); duration 4 hours
# run on 2 October 2020 for 15,742 genes in selection cohort only
#python3 interface.py main --dock $path_dock --distribution --remote --gene "ENSG00000186092"
#nohup python3 interface.py main --dock $path_dock --distribution --local > $path_dock/report.out &

# TODO: before candidacy procedure... determine which genes have significant distributions by permutations...

#python3 interface.py main --dock $path_dock --candidacy
#python3 interface.py main --dock $path_dock --stratification

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

#bash heritability_genes.sh # run on 2 October 2020
python3 interface.py main --dock $path_dock --heritability

#python3 interface.py main --dock $path_dock --category # <-- (2019-11-24) this procedure is obsolete for now...
#python3 interface.py main --dock $path_dock --prediction
#python3 interface.py main --dock $path_dock --function
#python3 interface.py main --dock $path_dock --integration

#python3 interface.py main --dock $path_dock --plot


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
