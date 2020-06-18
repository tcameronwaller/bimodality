#!/bin/bash

#chmod u+x script.sh

# Copy data attributes for samples and persons from original GTEx download.

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

# Origin
path_gtex="/nrnb/data/controlled/2020_dbGaP_GTEx/75875/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU"
path_origin_ancestry="/nrnb/users/mpagadal/gtex/gtex-ancestry/sklearn-gtex-ancestry"
path_origin_attributes="$path_gtex/PhenotypeFiles"
path_origin_attribute_sample="$path_origin_attributes/phs000424.v8.pht002743.v8.p2.c1.GTEx_Sample_Attributes.GRU.txt.gz"
path_origin_attribute_person="$path_origin_attributes/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz"
path_origin_explanation_sample="$path_origin_attributes/phs000424.v8.pht002743.v8.GTEx_Sample_Attributes.data_dict.xml"
path_origin_explanation_person="$path_origin_attributes/phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.xml"

# Destination
path_user_cellar="/cellar/users/tcwaller"
path_access_private="$path_user_cellar/Data/dock/access_private"
path_destination_ancestry="$path_access_private/persons_ancestry.txt"
path_destination_attribute_sample="$path_access_private/sample_attribute.txt"
path_destination_attribute_person="$path_access_private/person_attribute.txt"
path_destination_explanation_sample="$path_access_private/sample_explanation.xml"
path_destination_explanation_person="$path_access_private/person_explanation.xml"

# Organize destination directories
rm -r $path_access_private
mkdir -p $path_access_private

# Copy files from origins to destinations.
cp -r $path_origin_ancestry $path_destination_ancestry
cp -r $path_origin_attribute_sample $path_destination_attribute_sample
cp -r $path_origin_attribute_person $path_destination_attribute_person
cp -r $path_origin_explanation_sample $path_destination_explanation_sample
cp -r $path_origin_explanation_person $path_destination_explanation_person
