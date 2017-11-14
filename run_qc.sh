#!/bin/sh
#$ -cwd
#$ -pe threaded 1
#$ -S /bin/bash

## Set Plink (https://www.cog-genomics.org/plink2) Stable Beta 4.2 (3 November 2017)
plink2=/hpc/local/software/plink-1.9/plink

#########################
##### SAMPLE QC
##########################

## Here's the prefix for the data
data=jia.uviitis.hg19.merge

## Calculate missingness across the samples; drop samples with missingness > 5%
$plink2 --bfile $data \
    --missing \
    --out $data

## Prepare a set of high-quality SNPs for the rest of sample QC
## Identify SNPs in regions that can be highly stratified by population
awk ' $1==2 && $4 > 129883530 && $4 < 140283530 { print $2 } ' $data.bim > $data.pca.exclude.snps ## SNPs in the lactase locus
awk ' $1==6 && $4 > 24092021 && $4 < 38892022 { print $2 } ' $data.bim >> $data.pca.exclude.snps ## SNPs in the MHC
awk ' $1==8 && $4 > 6612592 && $4 < 13455629 { print $2 } ' $data.bim >> $data.pca.exclude.snps ## Inversion on chr8
awk ' $1==17 && $4 > 40546474 && $4 < 44644684 { print $2 } ' $data.bim >> $data.pca.exclude.snps ## Inversion on chr17
## Identify non-autosomal SNPs
awk ' $1 > 22 { print $2 } ' $data.bim >> $data.pca.exclude.snps ## Non-autosomal SNPs
## Identify A/T and C/G SNPs (if frequency is near 50%, strandedness can be tricky here)
cat $data.bim | awk ' ($5=="A" && $6=="T") || ($5=="T" && $6=="A") || ($5=="C" && $6=="G") || ($5=="G" && $6=="C") { print $2 } ' >> $data.pca.exclude.snps 

## Drop the SNPs we just identified, and keep the remaining data for sample QC
## Also keep MAF > 5% and missingness < 0.1% SNPs, and LD prune (r2 = 0.2)
$plink2 --bfile $data \
    --maf 0.05 \
    --geno 0.001 \
    --exclude ./$data.pca.exclude.snps \
    --indep-pairwise 50 5 0.2 \
    --out $data.pca

## PREPARATION FOR PCA
## Extract overlap with HapMap3 and merge everything together
## hapmap.jia.overlap.snps contains a list of SNPs that overlap between the JIA data and HapMap 3
$plink2 --bfile $data \
    --extract ./hapmap.jia.overlap.snps \
    --make-bed \
    --out $data.pca

## Extract the same set of SNPs from the HapMap 3 data (lifted over to hg19)
$plink2 --bfile /hpc/cog_gonl/sara/hapmap3/hg19/hapmap3_r2.hg19.dbsnp138.fwd.founders \
    --extract hapmap.jia.overlap.snps \
    --make-bed \
    --out hapmap3.pca

## Merge together the study data and the HapMap 3 data
$plink2 --bfile hapmap3.pca \
    --bmerge $data.pca.bed $data.pca.bim $data.pca.fam \
    --make-bed \
    --out jia.hapmap3.merge.pca

## Merging may throw an error at SNPs that may need their strand flipped
## Note: if strand flipping doesn't resolve the error, just drop the SNPs with --exclude
$plink2 --bfile hapmap3.pca \
    --flip jia.hapmap3.merge.pca-merge.missnp \
    --make-bed \
    --out hapmap3.pca

## Recode the data to ped/map file format
$plink2 --bfile jia.hapmap3.merge.pca \
    --recode \
    --out jia.hapmap3.merge.pca

## Run principal component analysis
## See run_smartPCA.sh in the Afib-Stroke-Overlap repository
## Drop PCA outliers
./run_smartPCA.sh jia.hapmap3.merge.pca

## Check for sex mismatches between phenotype and genotype data; drop any mismatches
$plink2 --bfile $data \
    --check-sex \
    --out $data

## Calculate heterozygozity/inbreeding of samples
## Drop all samples > 3sd from the distribution
$plink2 --bfile $data \
    --extract $data.pca.prune.in \
    --het \
    --out $data

## Calculate relatedness; drop samples with pi-hat > 0.125 (cousin relationship or higher)
$plink2 --bfile $data \
    --extract $data.pca.prune.in \
    --genome \
    --out $data

## Assemble all sample QC failures in the following file: jia.uviitis.hg19.merge.sample.failures (format: FID IID)

#################
##### SNP QC
#################

## Calculate missingness across SNPs; drop SNPs with missingness > 5%
$plink2 --bfile $data \
    --remove ./jia.uviitis.hg19.merge.sample.failures \
    --missing \
    --out $data

## Calculate Hardy-Weinberg; drop SNPs with HWE p < 1e-6 in controls only
$plink2 --bfile $data \
    --remove ./jia.uviitis.hg19.merge.sample.failures \
    --hardy \
    --pheno ./jia.uviitis.hg19.merge.pheno.jia_vs_healthy.txt \
    --out $data

## Differential missingness between cases and controls; drop SNPs with p < 1e-6
$plink2 --bfile $data \
    --remove ./jia.uviitis.hg19.merge.sample.failures \
    --test-missing \
    --pheno ./jia.uviitis.hg19.merge.pheno.jia_vs_healthy.txt \
    --out $data

## Check frequency distribution
$plink2 --bfile $data \
    --remove ./jia.uviitis.hg19.merge.sample.failures \
    --freq \
    --out $data
