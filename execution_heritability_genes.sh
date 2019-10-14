#!/bin/bash

# Suppress echo each command to console.
set +x

echo "--------------------------------------------------"
echo "----------"
echo "The script is for heritability analysis."
echo "----------"
echo "--------------------------------------------------"

# Echo each command to console.
set -x

# Organize paths.
path_gcta="/home/tcameronwaller/gcta_1.92.4beta/gcta64"

path_dock="/home/tcameronwaller/dock/"
path_gtex="$path_dock/gtex-8"
path_relation="$path_gtex/relation/autosome_common"
path_genes="$path_dock/split/genes.txt"
path_persons="$path_dock/selection/families_persons.tsv"
path_category="$path_dock/selection/persons_categories.tsv"
path_quantity="$path_dock/selection/persons_quantities.tsv"
path_distribution="$path_dock/distribution"
path_heritability="$path_dock/heritability"

rm -r $path_heritability
mkdir $path_heritability

# Suppress echo each command to console.
set +x

# Read genes.
readarray -t genes < $path_genes
# Report count of genes.
count_genes=${#genes[@]}
echo "count of genes: "
echo $count_genes
#echo ${genes[@]}
#echo $genes[0]

# Iterate on genes.
for gene in "${genes[@]}"
do
    echo "current gene: "
    echo $gene

    # Organize path.
    path_heritability_gene="$path_heritability/$gene"
    mkdir $path_heritability_gene
    path_simple="$path_heritability_gene/simple"
    path_complex="$path_heritability_gene/complex"
    mkdir $path_simple
    mkdir $path_complex

    # Access information about gene.
    path_distribution_gene="$path_distribution/$gene"
    path_phenotype="$path_distribution_gene/data_gene_families_persons_signals.tsv"

    # Execute heritability analysis.
    $path_gcta --grm $path_relation --keep $path_persons --pheno $path_phenotype --reml --out $path_simple/out --threads 5
    $path_gcta --grm $path_relation --keep $path_persons --pheno $path_phenotype --covar $path_category --qcovar $path_quantity --reml --out $path_complex/report --threads 5

done
