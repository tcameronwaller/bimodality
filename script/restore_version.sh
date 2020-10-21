#!/bin/bash

#chmod u+x script.sh

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous version of the program..."

rm -r bimodality-master/
rm master.zip

# Access current version of the program.

echo "access current version of the program..."

wget https://github.com/tcameronwaller/bimodality/archive/master.zip
unzip master.zip
cd bimodality-master/

