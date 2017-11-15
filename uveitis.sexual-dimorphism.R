####
#### Use an interaction test to see if there's a truly sexually-dimorphic signal at HLA-DRB1
####

## Read in the dosage data
## Header of this file (this is the same file used in uveitis.LRT.R, also in this repository)
## FID IID PC1 PC2 PC3 PC4 PC5 SEX PHASE PHENO AA_DRB1_11_32660115_SD.dosage ... [all additional columns are additional dosage data at MHC alleles]
data <- read.table("jia.uviitis.hg18.merged.covariates.amino.acids.txt",header=T)

## PHENO(uveitis case) = 2; PHENO(JIA non-uveitis case) = 1; PHENO(population-level control) = -9
data <- subset(data,data$PHENO!=-9) ## 192 uveitis cases and 330 non-uveitis JIA samples

## Standard association test for phenotype vs. HLA-DRB1 position 11 serine 
summary(glm(as.factor(PHENO) ~ PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(SEX) + PHASE + AA_DRB1_11_32660115_S.dosage, data=data, family="binomial"))

## Result:
## Coefficients:
##                              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                  -1.60317    0.30390  -5.275 1.32e-07 ***
## PC1                           2.52864    1.55325   1.628    0.104    
## PC2                          -1.82893    1.64854  -1.109    0.267    
## PC3                          -1.37088    1.77812  -0.771    0.441    
## PC4                          -0.81735    1.79260  -0.456    0.648    
## PC5                           0.38497    1.77819   0.216    0.829    
## as.factor(SEX)2              -0.08419    0.21036  -0.400    0.689    
## PHASE                        -0.16935    0.20441  -0.829    0.407    
## AA_DRB1_11_32660115_S.dosage  0.90518    0.14960   6.051 1.44e-09 ***

## Test of interaction for HLA-DRB1 serine at position 11 and sex
summary(glm(as.factor(PHENO) ~ PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(SEX) + PHASE + AA_DRB1_11_32660115_S.dosage + AA_DRB1_11_32660115_S.dosage*as.factor(SEX), data=data, family="binomial"))

## Result
## Coefficients:
##                                              Estimate Std. Error z value Pr(>|z|)   
## (Intercept)                                   -0.9487     0.3759  -2.524   0.0116 * 
## PC1                                            2.7764     1.5635   1.776   0.0758 . 
## PC2                                           -1.8223     1.6675  -1.093   0.2745   
## PC3                                           -1.3620     1.7984  -0.757   0.4488   
## PC4                                           -0.8445     1.8164  -0.465   0.6420   
## PC5                                            0.6126     1.7973   0.341   0.7332   
## as.factor(SEX)2                               -1.1654     0.4647  -2.508   0.0121 * 
## PHASE                                         -0.1416     0.2066  -0.685   0.4931   
## AA_DRB1_11_32660115_S.dosage                   0.3700     0.2479   1.492   0.1356   
## as.factor(SEX)2:AA_DRB1_11_32660115_S.dosage   0.8101     0.3128   2.590   0.0096 **

## Subset the data down to those samples that are poly RF-negative or oligo extended 
sample_subset <- read.table("polyRFneg_oligo.samples.list",header=T) ## list has columns FID, IID
data_subset <- merge(sample_subset,data,by="IID") ## N = 481 samples, 183 uveitis samples and 298 non-uveitis JIA samples

## Redo the original association test and test of interaction
## Association test with serine-11
summary(glm(as.factor(PHENO) ~ PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(SEX) + PHASE + AA_DRB1_11_32660115_S.dosage, data=data_subset, family="binomial"))

## Result
## Coefficients:
##                              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                   -1.1040     0.3317  -3.328 0.000875 ***
## PC1                            2.6253     1.6754   1.567 0.117132    
## PC2                           -1.6476     1.7785  -0.926 0.354227    
## PC3                           -1.6792     1.9250  -0.872 0.383025    
## PC4                           -1.7127     1.9174  -0.893 0.371734    
## PC5                            1.4209     1.9033   0.747 0.455329    
## as.factor(SEX)2               -0.3117     0.2233  -1.396 0.162711    
## PHASE                         -0.5185     0.2223  -2.332 0.019695 *  
## AA_DRB1_11_32660115_S.dosage   0.8856     0.1571   5.638 1.72e-08 ***

## Test of interaction
summary(glm(as.factor(PHENO) ~ PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(SEX) + PHASE + AA_DRB1_11_32660115_S.dosage + AA_DRB1_11_32660115_S.dosage*as.factor(SEX), data=data_subset, family="binomial"))

## Result
## Coefficients:
##                                              Estimate Std. Error z value Pr(>|z|)   
## (Intercept)                                   -0.4500     0.4051  -1.111   0.2667   
## PC1                                            2.8336     1.6913   1.675   0.0938 . 
## PC2                                           -1.7057     1.8048  -0.945   0.3446   
## PC3                                           -1.6164     1.9526  -0.828   0.4078   
## PC4                                           -1.9014     1.9396  -0.980   0.3269   
## PC5                                            1.5966     1.9207   0.831   0.4058   
## as.factor(SEX)2                               -1.4043     0.4879  -2.878   0.0040 **
## PHASE                                         -0.4972     0.2252  -2.208   0.0273 * 
## AA_DRB1_11_32660115_S.dosage                   0.3488     0.2581   1.352   0.1765   
## as.factor(SEX)2:AA_DRB1_11_32660115_S.dosage   0.8209     0.3270   2.510   0.0121 * 

## Is the p-value change specific to looking at these 183 uveitis samples and 298 non-uveitis JIA samples? 
## What if it were a random 183 uveitis cases and 298 non-uveitis JIA samples?

cases <- subset(data,data$PHENO==2)
controls <- subset(data,data$PHENO==1)

## Run a permutation; store the p-values of the interaction test from each permutation
pvals <- matrix(data=NA,ncol=1,nrow=1000)

## Run 1000 permutations
for(i in 1:1000) {
	
	## Randomly select 183 uveitis sampels
	case_subset <- cases[sample(1:192,183,replace=F),]
	control_subset <- controls[sample(1:330,298,replace=F),]
	
	permuted_data <- rbind(case_subset,control_subset)
	
	## Run the interaction test on the randomly-selected subset of cases and controls; store the interaction test p-value
	pvals[i,1] <- summary(glm(as.factor(PHENO) ~ PC1 + PC2 + PC3 + PC4 + PC5 + as.factor(SEX) + PHASE + AA_DRB1_11_32660115_S.dosage + AA_DRB1_11_32660115_S.dosage*as.factor(SEX), data=permuted_data, family="binomial"))$coef[[40]]
	
}

####
#### Plot: association testing in males vs. females in JIA-uveitis vs JIA non-uveitis
#### This can be done with the summary-level results denoted in the README file of this repository
####

## File containing summary-level results from female samples only
xx <- read.table("jia.uveitis.hg18.merged.QC.mhc.IMPUTED.phase_indicator.females.assoc.dosage.gz",header=T)

## File containing summary-level results from male samples only
xy <- read.table("jia.uveitis.hg18.merged.QC.mhc.IMPUTED.phase_indicator.males.assoc.dosage.gz",header=T)

## File containing summary-level results from female samples only, after conditioning on DRB1 position 233 residue T
xx.cond <- read.table("jia.uveitis.hg18.merged.QC.mhc.IMPUTED.phase_indicator.females.AA_DRB1_233_32656004_T.assoc.dosage.gz",header=T)

## Set up the plot so that it's one sex vs. the other
plot(0,col="white",xlim=c(29600000,33300000),axes=F,ylim=c(-10,10),xlab="Chromosome 6 position (Mb) (hg18)",ylab="-log10(p-value)"); box()
axis(side=2,at=c(seq(-9.5,-1.5,2),seq(1.5,9.5,2)),labels=c(rev(seq(1,9,2)),seq(1,9,2)),las=1)
axis(side=1,at=seq(29000000,33000000,1000000),labels=seq(29,33,1),las=1)

## Points from females
points(xx$BP,-log10(xx$P)+0.5,pch=23,bg="aliceblue",lwd=0.5)

## Points from males
points(xy$BP,log10(xy$P)-0.5,pch=23,bg="darkolivegreen2",lwd=0.5)

## Points from conditional testing in femlales
points(xx.cond$BP,-log10(xx.cond$P)+0.5,pch=23,bg="royalblue3",lwd=0.5)

## Add legends
legend("topright",legend="Females",cex=1.25,bty="n")
legend("topleft",pch=23,pt.bg=c("aliceblue","royalblue3"),legend=c("Initial association testing","Conditioned on DRB1 pos 233 (T)"),cex=0.75,bty="n")

legend("bottomright",legend="Males",cex=1.25,bty="n")
legend("bottomleft",legend="Initial association testing",pch=23,pt.bg="darkolivegreen2",cex=0.75,bty="n")
