#!/bin/sh
#$ -cwd
#$ -S /bin/bash
#$ -pe threaded 1
#$ -q short

## Input the chromosome, phenotype you want to analyze (see file passed to --pheno option in PLINK), window lower bound, window upper bound
## lower and upper correspond to the windows used in imputation
chr=$1
pheno=$2
lower=$3
upper=$4

## options for the "pheno" argument
## pheno: jia.healthy, uveitis.healthy, jia.uveitis

### The full Mb span of the chromosome
start=`awk ' $1=="'$chr'" { print $2 } ' /hpc/sara/jia/chromosome_windows.txt`
stop=`awk ' $1=="'$chr'" { print $3 } ' /hpc/sara/jia/chromosome_windows.txt`

## Find all the SNPs with imputation info score < 0.3; keep their identifiers so we can drop them in the GWAS
zcat ./imputation/chr$chr/info_by_snp/*.gz | awk ' $5 < 0.3 { print $2 } ' > ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.snps_info_below_0.3.drop

### Counter for each of the 5Mb intervals to be imputed
while [ $start -lt $stop ]; do

    lower=$start
    upper=$[$start+5]

    ## If the window of data exists, start the imputation
    if [ -f ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.gz ]; then

    echo $lower $upper
    ## Make a custom map file from the haps file (a bit of hack to get PLINK to run with a map file)
    zcat ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.gz | \
      awk ' { print "'$chr'", $2, 0, $3 } ' > \
        ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.map

    ## Now run the GWAS using Plink
    ## Be very careful here! You need to make sure that PLINK understands the data format, and that missing values (if there are any) are coded correctly.
    /hpc/local/plink2_dec14/plink \
	    --dosage ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.gz \ ## dosage data from imputation
	    noheader skip0=1 skip1=1 format=3 \ ## format of the data
	    --fam ./jia.uviitis.hg19.merge.QC.fam \ ## fam file with sample IDs
	    --map ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.map \ ## map file
	    --covar jia.uviitis.hg19.merge.QC.pca.pcs.txt \ ## file with covariates, containing columns "PC1", "PC2" and "SEX"
	    --covar-name PC1-PC2,SEX \ ## correct for PC1, PC2, and sex
	    --pheno ./jia.uviitis.hg19.merge.pheno.$pheno.txt \ ## pass PLINK the phenotype file
	    --out ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.$pheno ## write out the GWAS data
  
    ## Compress the GWAS data (these files are large)
    gzip ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.$pheno.assoc.dosage
    ## Remove the map file
    rm ./imputation/chr$chr/jia.uviitis.hg19.merge.QC.chr$chr.phased.$lower.$upper.imputed.map

    ## if there is no data covering this window of the genome, let us know!
    else
      echo "chr$chr $lower $upper NO DATA"
   fi

    ## update the start of the window and begin again
    start=$[$start+5]

done
