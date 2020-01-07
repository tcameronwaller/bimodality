#!/bin/bash

#chmod u+x script.sh

# Copy data attributes for samples and persons from original GTEx download.

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

# Origin
path_pagadalla="/nrnb/users/mpagadal"
path_origin_ancestry="$path_pagadalla/gtex-ancestry/sklearn-gtex-ancestry"
path_origin_gtex_8="$path_pagadalla/gtex/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU"
path_origin_genotype="$path_origin_gtex_8/GenotypeFiles/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz"
path_origin_attributes="$path_origin_gtex_8/PhenotypeFiles"
path_origin_attribute_sample="$path_origin_attributes/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt"
path_origin_attribute_person="$path_origin_attributes/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"
path_origin_explanation_sample="$path_origin_attributes/phs000424.v8.pht002743.v8.GTEx_Sample_Attributes.data_dict.xml"
path_origin_explanation_person="$path_origin_attributes/phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.xml"

# Destination
path_user_cellar="/cellar/users/tcwaller"
path_user_nrnb="/nrnb/users/tcwaller"
path_destination_gtex_8="$path_user_nrnb/gtex-8"
path_destination_genotype="$path_destination_gtex_8/gtex-8_genotype.vcf.gz"
path_destination_dock_gtex_8="$path_user_cellar/Data/dock/gtex-8"
path_destination_attribute_sample="$path_destination_dock_gtex_8/sample_attribute.txt"
path_destination_attribute_person="$path_destination_dock_gtex_8/person_attribute.txt"
path_destination_explanation_sample="$path_destination_dock_gtex_8/sample_explanation.xml"
path_destination_explanation_person="$path_destination_dock_gtex_8/person_explanation.xml"

# Organize destination directories
rm -r $path_destination_gtex_8
mkdir $path_destination_gtex_8
rm -r $path_destination_dock_gtex_8
mkdir $path_destination_dock_gtex_8

# Copy files from origins to destinations.
cp -r $path_origin_genotype $path_destination_genotype
cp -r $path_origin_attribute_sample $path_destination_attribute_sample
cp -r $path_origin_attribute_person $path_destination_attribute_person
cp -r $path_origin_explanation_sample $path_destination_explanation_sample
cp -r $path_origin_explanation_person $path_destination_explanation_person
