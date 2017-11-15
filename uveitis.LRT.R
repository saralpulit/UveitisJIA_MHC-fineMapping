## Load in the data
## The header of this file looks like:
## FID IID PC1 PC2 PC3 PC4 PC5 SEX PHASE PHENO AA_DRB1_11_32660115_SD.dosage AA_DRB1_67_32659947_F.dosage AA_DRB1_11_32660115_S.dosage AA_DRB1_11_32660115_D.dosage [... additional columns of HLA-DRB1 residue dosages]

data <- read.table("jia.uviitis.hg18.me rged.covariates.amino.acids.txt",header=T)

## PHENO is coded as 2 (uveitis-JIA), 1 (non-uveitis JIA) or -9 (other)
data <- subset(data,data$PHENO!=-9)

#### Here is the base formula (phenotype and covariates)
covariates <- "PC1 + PC2 + PC3 + PC4 + PC5 + SEX + PHASE"
f0 <- formula(paste("as.factor(PHENO) ~ ", covariates, sep=""))

#### Test this compared to including either pos 10 (T), pos 11 (S), or pos 12 (Y)
#### These should all be the same (LD = 1 between all residue pairs)

result1 <- matrix(ncol=2,nrow=8,data=0)
colnames(result1) <- c("Residue","Goodness.P")

counter <- 1

## Point to all of the columns in the data that refer to an HLA-DRB1 amino acid residue, to test in the likelihood ratio test
for(aa in 11:18) {

	## Pick out the residue from the data
	residue <- colnames(data)[aa]
	
	## Build the second model, which is the covariates plus the residue
	f1 <- formula(paste("as.factor(PHENO) ~ ", covariates, " + ", residue, sep=""))
	print(f1)
	
	## Now run the regression for each of the two models; store it in m0 and m1, respectively		
	m0 <- glm(f0, data=data, family="binomial")
	m1 <- glm(f1, data=data, family="binomial")

	## Perform the goodness of fit (likelihood ratio) test
	goodness.test <- anova(m0,m1,test="Chisq")$P[2]
	print(goodness.test)
	
	## Store the result of the test
	result1[counter,1] <- residue; result1[counter,2] <- goodness.test
	
	counter <- counter+1
			
}

## Results for all the residues at positions 10-13
#    Residue                         Goodness.P            
#	[1,] "AA_DRB1_11_32660115_S.dosage"  "1.47596057315658e-10"
# 	[2,] "AA_DRB1_11_32660115_D.dosage"  "0.528814541985092"   
#	[3,] "AA_DRB1_13_32660109_SG.dosage" "1.46737885272758e-10"
#	[4,] "AA_DRB1_13_32660109_S.dosage"  "1.21318105498409e-05"
#	[5,] "AA_DRB1_13_32660109_G.dosage"  "0.00207722539784343" 
#	[6,] "AA_DRB1_10_32660118_Y.dosage"  "1.47596057315658e-10"
#	[7,] "AA_DRB1_12_32660112.dosage"    "1.47596057315658e-10"
#	[8,] "AA_DRB1_11_32660115_SD.dosage" "3.009625e-11"


## Now let's choose the model that also includes serine at position 11 as m0
f0 <- formula(paste("as.factor(PHENO) ~ ", covariates, " + AA_DRB1_11_32660115_S.dosage", sep=""))

#### Repeat the likelihood ratio testing to see if there's any better model
counter <- 1

result2 <- matrix(ncol=2,nrow=8,data=0)
colnames(result2) <- c("Residue","Goodness.P")

for(aa in 11:18) {

	residue <- colnames(data)[aa]
	f1 <- formula(paste("as.factor(PHENO) ~ ", covariates, " + AA_DRB1_11_32660115_S.dosage +", residue, sep=""))
	print(f1)
			
	m0 <- glm(f0, data=data, family="binomial")
	m1 <- glm(f1, data=data, family="binomial")

	goodness.test <- anova(m0,m1,test="Chisq")$P[2]
	print(goodness.test)
	
	result2[counter,1] <- residue
	result2[counter,2] <- goodness.test	
	
	counter <- counter+1
			
}

## Results for all residues at positions 10-13
#     Residue                         Goodness.P             
# [1,] "AA_DRB1_10_32660118_Y.dosage"  NA                     ## NA because of LD with serine-11
# [2,] "AA_DRB1_11_32660115_SD.dosage" "0.0757828440971554"   
# [3,] "AA_DRB1_11_32660115_S.dosage"  NA                     ## NA because serine-11 is already in the model
# [4,] "AA_DRB1_11_32660115_D.dosage"  "0.0764245163291888"   
# [5,] "AA_DRB1_12_32660112.dosage"    NA                     ## NA because of LD with serine-11
# [6,] "AA_DRB1_13_32660109_SG.dosage" "0.043064549768632"    
# [7,] "AA_DRB1_13_32660109_S.dosage"  "0.699680114306571"    
# [8,] "AA_DRB1_13_32660109_G.dosage"  "0.700453481628917"    

