#read the CMH result on empirical data
initial <- read.table('F0F60SNP.sync.cmh',header=FALSE, sep='\t',nrows=100)
classes <- sapply(initial,class)
CMHF0F60 <- read.table('F0F60SNP.sync.cmh',header=FALSE, sep='\t',colClasses=classes)

#read the CMH result on simulated data: autosomes
initial <- read.table('Dsim_Ne291_autosome_4539060independentSNPs_10reps_simulated.sync.cmh',header=FALSE, sep='\t',nrows=100)
classes <- sapply(initial,class)
CMHF0F60_null_auto <- read.table('Dsim_Ne291_autosome_4539060independentSNPs_10reps_simulated.sync.cmh',header=FALSE, sep='\t',colClasses=classes)

#read the CMH result on simulated data: X chromosome
initial <- read.table('Dsim_Ne262_X_557140independentSNPs_10reps_simulated.sync.cmh',header=FALSE, sep='\t',nrows=100)
classes <- sapply(initial,class)
CMHF0F60_null_X <- read.table('Dsim_Ne262_X_557140independentSNPs_10reps_simulated.sync.cmh',header=FALSE, sep='\t',colClasses=classes)

######################################################
# FDR corrected q-values of the CMH test for autosomes
######################################################

cmhNull <- CMHF0F60_null_auto$V24 #the 24th column of CMH result corresponds to -log10(p-value)
cmhF60 <- CMHF0F60$V24[CMHF0F60$V1 != 'X']

#this function computes FDR corrected q-values. emp_pval: -log10(pvalue) of empirical data,sim_pval: -log10(pvalue) of simulated data, thresh: FDR threshold, e.g. for 5% provide 0.05
compute_FDR <- function(emp_pval,sim_pval,thresh){
  ecdf_sim <- ecdf(sim_pval)
  ecdf_emp <- ecdf(emp_pval)
  fdr <- (1-ecdf_sim(emp_pval)) / (1-ecdf_emp(emp_pval))
  fdr[emp_pval > max(sim_pval)] <- 0
  fdr[sim_pval < min(emp_pval)] <- 1
  candidates <- emp_pval[fdr <= thresh]
  FDR <- min(candidates, na.rm = TRUE)
  return(FDR)
}

#compute 5% FDR
#SNPs with -log10(p-value) values above this threshold have allele frequency change that is more than expected under drift
FDR <- compute_FDR(cmhF60,cmhNull,0.05)
#########################################################
# FDR corrected q-values of the CMH test for X chromosome
#########################################################

cmhNull <- CMHF0F60_null_X$V24 #the 24th column of CMH result corresponds to -log10(p-value)
cmhF60 <- CMHF0F60$V24[CMHF0F60$V1 == 'X']

#this function computes FDR corrected q-values. emp_pval: -log10(pvalue) of empirical data,sim_pval: -log10(pvalue) of simulated data, thresh: FDR threshold, e.g. for 5% provide 0.05
compute_FDR <- function(emp_pval,sim_pval,thresh){
  ecdf_sim <- ecdf(sim_pval)
  ecdf_emp <- ecdf(emp_pval)
  fdr <- (1-ecdf_sim(emp_pval)) / (1-ecdf_emp(emp_pval))
  fdr[emp_pval > max(sim_pval)] <- 0
  fdr[sim_pval < min(emp_pval)] <- 1
  candidates <- emp_pval[fdr <= thresh]
  FDR <- min(candidates, na.rm = TRUE)
  return(FDR)
}

#compute 5% FDR
#SNPs with -log10(p-value) values above this threshold have allele frequency change that is more than expected under drift
FDR <- compute_FDR(cmhF60,cmhNull,0.05)
