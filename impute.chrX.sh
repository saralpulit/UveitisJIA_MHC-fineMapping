#!/bin/sh
#$ -cwd
#$ -pe threaded 2
#$ -S /bin/bash

## To impute, we have chunked the data up by chromosome and into 5Mb windows (so that e.g., start is "0" and stop is "5")
chr=$1
start=$2
stop=$3

## Point to IMPUTE2, the genetic map, and the reference panel we are going to use for imputation (1000 Genomes Phase 3)
impute2=/hpc/software/impute_v2.3.1_x86_64_static/impute2
gen_map=/hpc/resources/1000Genomes/Phase3/IMPUTE2_format
kgp=/hpc/resources/1000Genomes/Phase3/IMPUTE2_format

## here is the home directory we are working from
home=/hpc/sara/jia/imputation

## Print the data and time
date +"%m/%d/%Y %H:%M:%S"

## Call IMPUTE2
$impute2 \
    -m $gen_map/genetic_map_chrX_nonPAR_combined_b37.txt \ ## genetic map
    -h $kgp/1000GP_Phase3_chrX_NONPAR.hap.gz \ ## 1KG Phase 3 data
    -l $kgp/1000GP_Phase3_chrX_NONPAR.legend.gz \ ## 1KG Phase 3 data
    -known_haps_g $home/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.haps.gz \ ## the phased study data from SHAPEIT
    -chrX \ ## the specific chrX flag
    -int "$start"e6 "$stop"e6 \ ## start and stop of the imputation window
    -Ne 20000 \ ## effective sample size; see IMPUTE2 documentation for details
    -use_prephased_g \ ## let IMPUTE2 know you're using phased data
    -sample_g $home/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.samples \ ## study samples
    -k_hap 1000 \ ## number of haplotypes to search in imputation; set to roughly the number of European-ancestry haplotypes in 1KG Phase 3
    -buffer 250 \ ## buffer space around the imputation window (in kb)
    -r $home/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$start.$stop.log \ ## log file
    -o $home/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$start.$stop.imputed ## output

## Print the data and time
date +"%m/%d/%Y %H:%M:%S"
