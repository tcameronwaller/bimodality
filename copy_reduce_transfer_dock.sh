#!/bin/bash

# cp -r /media/tcameronwaller/primary/data/local/work/project/2019_bimodality/dock_template /home/tcameronwaller/dock_template
# mv /home/tcameronwaller/dock_template /home/tcameronwaller/dock


#chmod u+x script.sh

# Execution
# $ bash execution_halyard.sh

path_home="/home/tcameronwaller"
path_dock="/home/tcameronwaller/dock"
path_dock_lite="/home/tcameronwaller/dock_lite"

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
# Copy and reduce dock

rm -r $path_dock_lite
mkdir $path_dock_lite
cp -r "$path_dock/selection" "$path_dock_lite/selection"
tar -czvf "$path_home/dock_lite.tar.gz" $path_dock_lite
scp "$path_home/dock_lite.tar.gz" tcwaller@grenache.ucsd.edu:/cellar/users/tcwaller/Data/dock_lite.tar.gz

##########
# Organize dock on server

#tar -xzvf dock_lite.tar.gz
#mv /cellar/users/tcwaller/Data/dock_lite/ /cellar/users/tcwaller/Data/dock/


##########
# Transfer from server to halyard
#tar -czvf "$path_home/dock_lite.tar.gz" $path_dock_lite
#scp tcwaller@grenache.ucsd.edu:/cellar/users/tcwaller/Data/dock_lite.tar.gz dock_lite.tar.gz
