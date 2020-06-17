#!/bin/bash

#chmod u+x script.sh

# Copy data attributes for samples and persons from original GTEx download.

# Organize paths.
export PATH=/cellar/users/tcwaller/anaconda3/bin:$PATH

# Origin
path_gtex="/nrnb/data/controlled/2020_dbGaP_GTEx/75875/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU"

# Remove old versions of genotype data
rm -r "$path_gtex/GenotypeFiles/phg000219.v3.GTEx_Illumina5M.genotype-original-submission.c1"
rm -r "$path_gtex/GenotypeFiles/phg000219.v3.GTEx_IlluminaExome.genotype-original-submission.c1"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_Illumina5M.genotype-calls-matrixfmt.c1"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_Illumina5M.genotype-calls-matrixfmt.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_Illumina5M.genotype-original-submission.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_IlluminaExome.genotype-calls-matrixfmt.c1"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_IlluminaExome.genotype-calls-matrixfmt.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_IlluminaExome.genotype-original-submission.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_Imputation.genotype-imputed-data.c1"
rm -r "$path_gtex/GenotypeFiles/phg000219.v4.GTEx_Pilot_Imputation.genotype-imputed-data.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Illumina25M.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Illumina25M.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Illumina25M.genotype-original-submission.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Illumina25M.genotype-original-submission.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Illumina25M.raw-data-idat.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Illumina25M.raw-data-idat.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WES_LoF_CNV.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WES_LoF_CNV.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WES_SNP_CNV.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WES_SNP_CNV.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000520.v2.GTEx_MidPoint_WGS_SNP_CNV.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WES.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WES.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WES.panel-of-normals.c1"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WES.panel-of-normals.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WGS_additional.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WGS_additional.genotype-calls-vcf.c1.GRU.tar"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1"
rm -r "$path_gtex/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar"
