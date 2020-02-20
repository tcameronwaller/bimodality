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
path_gcta="/home/tcameronwaller/gcta_1.93.0beta/gcta64"

path_dock="/home/tcameronwaller/dock"
path_distribution="$path_dock/distribution/genes"
path_heritability="$path_dock/heritability/genes"

path_genes="$path_dock/selection/tight/genes_selection.txt"
path_persons="$path_dock/selection/tight/heritability/simple/families_persons.tsv"
path_simple_variables="$path_dock/selection/tight/heritability/simple/persons_variables.tsv"
path_complex_variables="$path_dock/selection/tight/heritability/complex/persons_variables.tsv"

path_relation_gcta="$path_dock/access_private/relation/gcta"
path_relation_gcta_bed_bim_fam="$path_relation_gcta/bed_bim_fam/autosome_common"
path_relation_gcta_pgen_pvar_psam="$path_relation_gcta/pgen_pvar_psam/autosome_common"
# Use either the Genetic Relationship Matrix from old or new formats.
path_relation=$path_relation_gcta_bed_bim_fam
#path_relation=$path_relation_gcta_pgen_pvar_psam

rm -r $path_heritability
mkdir -p $path_heritability

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
    $path_gcta --grm $path_relation --keep $path_persons --pheno $path_phenotype --qcovar $path_simple_variables --reml --out $path_simple/report --threads 5
    $path_gcta --grm $path_relation --keep $path_persons --pheno $path_phenotype --qcovar $path_complex_variables --reml --out $path_complex/report --threads 5
done
