#!/bin/sh
#$ -cwd
#$ -pe threaded 1
#$ -S /bin/bash

## Load plink (or point to the software in a directory)
module load plink

plink --noweb \
    --dosage ./jia.uveitis.hg18.merged.QC.mhc.IMPUTED.dosage.gz \ ## point to the dosages
    format=1 noheader \ ## state the format; this is the format output by the SNP2HLA pipeline
    --fam jia.uveitis.hg18.merged.QC.mhc.IMPUTED.fam \ ## load the fam file (list of sample identifiers)
    --pheno ./jia.uveitis.hg18.merged.QC.mhc.IMPUTED.uveitis.jia.pheno \ ## load the phenotype
    --out ./jia.uveitis.hg18.merged.QC.mhc.IMPUTED.uveitis.jia.phase_indicator \ ## name the output
    --map ./jia.uveitis.hg18.merged.QC.mhc.IMPUTED.map \ ## load a map file
    --covar ./jia.uveitis.hg18.merged.covariates.amino.acids.txt \ ## load the covariate file
    --covar-name PC1,PC2,PC3,PC4,PC5,SEX,PHASE ## list the covariates to use in the analysis
