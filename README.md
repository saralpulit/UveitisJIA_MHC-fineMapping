# MHC fine-mapping in juvenile idiopathic arthritis (JIA)-associated uveitis

Scripts and data used in analyses for "An amino acid motif in HLA-DRβ1 distinguishes patients with uveitis in juvenile idiopathic arthritis."

This repository contains:

## Data

Summary-level datasets are available here:   
Phase 1 genome-wide association summary-level results: https://doi.org/10.5281/zenodo.1048977   
Phase 2 genome-wide association summary-level results: https://doi.org/10.5281/zenodo.1048979

#### chromosome_windows.txt
This text file can be useful for splitting up data into 5Mb windows or running imputation/GWAS in 5Mb chunks. The file contains:
 1. The chromosome number
 2. The starting position (in Mb) on that chromosome
 3. The ending position (in Mb) on that chromosome
 4. The total number of 5Mb chunks on that chromosome.

#### metal.mhc.txt
A parameters file to run a meta-analysis using METAL (https://genome.sph.umich.edu/wiki/METAL).

Usage:   
```/path/to/metal metal.mhc.txt```

## Scripts

#### run_qc.sh
A bash script containing a series of Plink commands used to run sample and variant QC

#### run_smartPCA.sh
Please see the Afib-Stroke-Overlap repository for this script. It uses EIGENSTRAT (https://www.hsph.harvard.edu/alkes-price/eigensoft-frequently-asked-questions/ to run principal component analysis (PCA) either using a reference set of data or just in your own samples.

Usage:   
```./run_smartPCA.sh data_prefix```

#### prephase.sh
A bash script to call SHAPEIT, which will phase your data before imputation.

Usage:   
```./prephase.sh chromosome```

#### impute.sh
A bash script to call IMPUTE2, which will impute your phased data (produced by SHAPEIT)

Usage:   
```./impute.sh chromosome window_start window_stop```

#### impute.chrX.sh
A bash script to call IMPUTE2, which will impute your phased data (produced by SHAPEIT) specifically for the X chromosome

Usage:   
```./impute.sh chromosome window_start window_stop```

#### run_gwas.sh
A bash script that will point Plink to the imputed data and run a genome-wide association study

Usage:   
```./run_gwas.sh chromosome phenotype window_start window_stop```







