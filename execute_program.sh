#!/bin/bash

#chmod u+x script.sh

path_project="/media/tcameronwaller/primary/data/local/work/project/2019_bimodality/"
subpath_repository="repository/"
path_repository="$path_project$subpath_repository"
subpath_program="repository/bimodality/"
path_program="$path_project$subpath_program"
subpath_dock="dock/"
path_dock="$path_project$subpath_dock"

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

#python3 interface.py main -d $path_dock -a
#python3 interface.py main -d $path_dock -b
#python3 interface.py main -d $path_dock -s
python3 interface.py main -d $path_dock -o
#python3 interface.py main -d $path_dock -n
#python3 interface.py main -d $path_dock -m
#python3 interface.py main -d $path_dock -t
