# MHC fine-mapping in juvenile idiopathic arthritis (JIA)-associated uveitis

Scripts and data used in analyses for *"An amino acid motif in HLA-DRβ1 distinguishes patients with uveitis in juvenile idiopathic arthritis." (Preprint: https://www.biorxiv.org/content/early/2017/08/14/140954)*

This repository contains:

## Data

### Summary-level datasets are available here: 
 
 *Genome-wide data*
 1. Phase 1 genome-wide association summary-level results: https://doi.org/10.5281/zenodo.1048977   
 2. Phase 2 genome-wide association summary-level results: https://doi.org/10.5281/zenodo.1048979
 
  *MHC-specific data*
 1. MHC-wide association results comparing uveitis-JIA cases to non-uveitis JIA samples: https://doi.org/10.5281/zenodo.1049020   
 2. MHC-wide association results (uveitis-JIA cases vs. non-uveitis JIA samples), conditioning on aspartic acid (D) at position 11 in HLA-DRB1:https://doi.org/10.5281/zenodo.1049023
 3. MHC-wide association results (uveitis-JIA cases vs. non-uveitis JIA samples), conditioning on serine (S) at position 11 in HLA-DRB1: https://doi.org/10.5281/zenodo.1049025
 4. MHC-wide association results comparing all JIA cases versus disease-free controls: https://doi.org/10.5281/zenodo.1053793
 5. MHC-wide association results comparing uveitis-JIA cases versus disease-free controls: https://doi.org/10.5281/zenodo.1053797

   *MHC-specific data, split by sex*
 1. MHC-wide association results (uveitis-JIA cases vs. non-uveitis JIA samples) in *female samples only*: https://doi.org/10.5281/zenodo.1049027   
 2. MHC-wide association results (uveitis-JIA cases vs. non-uveitis JIA samples) in *male samples only*: https://doi.org/10.5281/zenodo.1049031
 
### Supplementary Tables

#### SuppTable2.HLAAssociations.csv
Supplementary Table 2 | Association results for HLA alleles with frequency >0.01 in JIA-uveitis cases versus JIA non-uveitis controls. For each amino acid, the frequency in cases and controls, the imputation-score (INFO; a score between 0 and 1, a higher value means a better imputation of that amino acid), odds ratio (OR), 95%-confidence interval (95%-CI) and p-values are given. 

#### SuppTable4.ProteinsforBindingAffinity.csv
Supplementary Table 4 | One-hundred and forty-seven proteins selected from iris tissue proteome data for peptide binding affinity prediction to common HLA-DRβ1 molecules. See Materials and Methods in the main article for methodological details. Mass spectrometric proteome data from human iris tissues (2,959 nonredundant proteins) was used as a representative source of proteins present in iris tissue. Protein accession numbers were extracted to filter in UniProt (Universal Protein Resource) for proteins with a subcellular location annotation limited to the nucleus. 147 proteins fulfilled these criteria and their full length amino acid sequences were fed into the neural network of the netMHCIIpan3.1 server.

### Other supporting files

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

#### 1. run_qc.sh    
A bash script containing a series of Plink commands used to run sample and variant QC

#### 2. run_smartPCA.sh
Please see the Afib-Stroke-Overlap repository for this script. It uses EIGENSTRAT (https://www.hsph.harvard.edu/alkes-price/eigensoft-frequently-asked-questions/ to run principal component analysis (PCA) either using a reference set of data or just in your own samples.

Usage:   
```./run_smartPCA.sh data_prefix```

#### 3. prephase.sh
A bash script to call SHAPEIT, which will phase your data before imputation.

Usage:   
```./prephase.sh chromosome```

#### 4. impute.sh
A bash script to call IMPUTE2, which will impute your phased data (produced by SHAPEIT)

Usage:   
```./impute.sh chromosome window_start window_stop```

#### 5. impute.chrX.sh
A bash script to call IMPUTE2, which will impute your phased data (produced by SHAPEIT) specifically for the X chromosome

Usage:   
```./impute.sh chromosome window_start window_stop```

#### 6. run_gwas.sh
A bash script that will point Plink to the imputed data and run a genome-wide association study

Usage:   
```./run_gwas.sh chromosome phenotype window_start window_stop```

#### 7. SNP2HLA.csh
A script for running the SNP2HLA imputation pipeline. Script authored by Xiaoming Jia

Usage:  
```./SNP2HLA.csh DATA (.bed/.bim/.fam) REFERENCE (.bgl.phased/.markers) OUTPUT plink {optional: java_max_memory[mb] marker_window_size}```

DATA: bim/bed/fam files of your data   
REFERENCE: your reference dataset (here, provided by the Type 1 Diabetes Genetics Consortium)   
OUTPUT: name of your output files    

#### 8. run_imputeHLA.sh
A bash script to run the HLA imputation.

Uage:
```./run_imputeHLA.sh input reference_data```

#### 9. mhc.assoc.sh
A script to run association testing (using PLINK) on the imputed MHC data

#### 10. uveitis.LRT.R
Code used for calculating likelihood ratio tests in the uveitis data (requires individual-level data)

#### 11. uveitis.sexual-dimorphism.R
Code used for testing the interaction of sex and the HLA-DRB1 serine-11 genotype (requires individual-level data) and for plotting association results in male vs. female samples (requires summary-level data; please see above for links).


