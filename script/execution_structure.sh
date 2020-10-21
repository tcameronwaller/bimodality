#!/bin/bash

# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "This script accesses and organizes HiC-seq data for visualization."
echo "----------"
echo "--------------------------------------------------"

# Echo each command to console.
set -x

# Organize paths.
path_juicer="/home/tcameronwaller/Downloads/juicer_tools_1.14.08.jar"
path_download="/home/tcameronwaller/Downloads"
path_dock="/home/tcameronwaller/dock"
path_access="$path_dock/access"
path_structure="$path_dock/structure"

rm -r $path_structure
mkdir $path_structure

# Suppress echo each command to console.
set +x

##########
# Access source data.
# GSE95014
# GSM2494290
# GSM2494294
# GSM2494298
# Authors processed reads from HiC-seq in HiC-Pro.
# Data are a list of interactions between genomic fragments.

wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE95014&format=file&file=GSE95014_Hap1.validPairs.txt.gz
cp $path_download/GSE95014_Hap1.validPairs.txt.gz $path_access/GSE95014_Hap1.validPairs.txt.gz

##########
# Sort list of genomic interactions.
# 1. First chromosome must increase throughout file
# 2. First chromosome number must be less than second chromosome number

zcat $path_access/GSE95014_Hap1.validPairs.txt.gz | awk 'BEGIN{OFS="\t"} {if ($2 > $5){ print $1, $5, $6, $7, $2, $3, $4, $8}else { print $1, $2, $3, $4, $5, $6, $7, $8}}' | sort -k2,2d -k5,5d | gzip > $path_structure/data_pairs_sort_raw.txt.gz

##########
# Change format for compatibility with Juicer Tools Pre.

# Change designations of forward and reverse strands.
zcat $path_structure/data_pairs_sort_raw.txt.gz | sed -e 's/+/0/g' | sed -e 's/-/1/g' | gzip > $path_structure/data_pairs_sort_strand.txt.gz

# Arrange columns and insert dummy columns for restriction fragments.
zcat $path_structure/data_pairs_sort_strand.txt.gz | awk 'BEGIN{OFS="\t"} {print $4, $2, $3, 0, $7, $5, $6, 1}' | gzip > $path_structure/data_pairs_sort_format.txt.gz

##########
# Run Juicer Tools Pre function.
# Aggregate contacts by bins on basis of chromosomal regions.
# Compile high-efficiency binary matrix file.

java -Xms10g -Xmx10g -jar $path_juicer pre -d $path_structure/data_pairs_sort_format.txt.gz $path_structure/contacts.hic hg19
