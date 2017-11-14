#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -pe threaded 5

## To run SHAPEIT, you need to point to the software, a genetic map, and (optionally) a reference dataset
## Note that we did not use a reference panel in this instance, because N > 100
shapeit=/hpc/local/software/shapeit_v2/shapeit.v2.r644.linux.x86_64
gen_map=/hpc/resources/1000Genomes/Phase3/IMPUTE2_format
kgp=/hpc/resources/1000Genomes/Phase3/IMPUTE2_format

## Make chromosome a command-line argument, so we can run in parallel
chr=$1

## Point to the data (in PLINK format) that we want to phase
input=/hpc/jia/imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr

## Print the data and time (to keep track of CPU time)
date +"%m/%d/%Y %H:%M:%S"

## Run without a reference panel (> 100 samples)

## Impute the autosome
if [ $chr -ne 23 ]; then

    $shapeit --input-bed $input.bed $input.bim $input.fam \
	    --input-map $gen_map/genetic_map_chr"$chr"_combined_b37.txt \
	    --output-max $input.phased.haps $input.phased.samples \
	    --thread 5 \
	    --output-log $input.shapeit2

fi

## Impute the X chromosome
if [ $chr -eq 23 ]; then

  ## SNPs are only in the nonPAR region

   $shapeit --input-bed $input.bed $input.bim $input.fam \
	    --input-map $gen_map/genetic_map_chrX_nonPAR_combined_b37.txt \
	    --output-max $input.phased.haps $input.phased.samples \
	    --chrX \
	    --thread 5 \
	    --output-log $input.shapeit2

fi

## Print the data and time
date +"%m/%d/%Y %H:%M:%S"
