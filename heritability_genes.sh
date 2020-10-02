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

# Define cohort of persons.
cohort="selection"

# Organize paths.
path_gcta="/home/tcameronwaller/gcta_1.93.2beta/gcta64"

path_dock="/home/tcameronwaller/dock"
path_distribution="$path_dock/distribution/$cohort/genes"
path_heritability="$path_dock/heritability/$cohort/genes"

#path_genes="$path_dock/selection/tight/samples_genes_signals/genes.txt" # <- needs to be distribution genes...
path_persons="$path_dock/selection/tight/persons_properties/$cohort/heritability/families_persons.tsv"
path_variables="$path_dock/selection/tight/persons_properties/$cohort/heritability/persons_variables.tsv"

path_relation_gcta="$path_dock/access_private/relation/$cohort/gcta"
path_relation_gcta_bed_bim_fam="$path_relation_gcta/bed_bim_fam/autosome_common"
#path_relation_gcta_pgen_pvar_psam="$path_relation_gcta/pgen_pvar_psam/autosome_common"
# Use either the Genetic Relationship Matrix from old or new formats.
path_relation=$path_relation_gcta_bed_bim_fam # <- this format works for GCTA
#path_relation=$path_relation_gcta_pgen_pvar_psam # <- this format does not work

rm -r $path_heritability
mkdir -p $path_heritability

# Suppress echo each command to console.
set +x

# Read genes.
#readarray -t genes < $path_genes
genes=($(ls $path_distribution))
# Report count of genes.
count_genes=${#genes[@]}
echo "count of genes: "
echo $count_genes
#echo ${genes[@]}
#echo $genes[0]

# Iterate on genes.
#for gene in "${genes[@]}"
for gene in "${genes[@]:0:100}"
do
    echo "current gene: "
    echo $gene

    # Organize path.
    path_heritability_gene="$path_heritability/$gene"
    mkdir $path_heritability_gene

    # Access information about gene.
    path_distribution_gene="$path_distribution/$gene"
    path_phenotype="$path_distribution_gene/data_gene_families_persons_signals.tsv"

    # Execute heritability analysis.
    $path_gcta --grm $path_relation --keep $path_persons --pheno $path_phenotype --qcovar $path_variables --reml --out $path_heritability_gene/report --threads 8
done
