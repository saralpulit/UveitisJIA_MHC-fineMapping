#!/bin/bash
#$ -cwd
#$ -l h_rt=36:00:00,h_vmem=100G
#$ -pe threaded 10

plink=/hpc/local/software/plink-1.90/plink

INPUT=$1
REFERENCE=$2

## Point to the imptuation pipeline, included in this repository
IMPUTE=/hpc/hers_en/sara/resources/mhc_imputation/ImputeHLA/SNP2HLA.csh

## Run imputation
$IMPUTE $INPUT $REFERENCE "$INPUT"_IMPUTED $plink 20000 1000
