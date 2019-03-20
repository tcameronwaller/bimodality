#!/bin/bash

#chmod u+x script.sh

path_project="/media/tcameronwaller/primary/data/local/work/project/2019_bimodality/"
directory_repository="bimodality"
path_repository="$path_project$directory_repository/"
directory_dock="dock"
path_dock="$path_project$directory_dock/"

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
mkdir $directory_dock
cd $path_repository

echo "Now set to call routine and procedures."

#python3 interface.py main -d $path_dock -a
python3 interface.py main -d $path_dock -o

