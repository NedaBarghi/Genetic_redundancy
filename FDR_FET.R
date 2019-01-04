###########################################################################
#read the -log10(pvalue) of Fisher's exact test on simulated data autosomes
###########################################################################
#replicate1
Nullauto_rep1 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep1.sync.fet',header=FALSE, sep='\t')

Nullauto_rep1$V6 <- as.character(Nullauto_rep1$V6)
simAuto_rep1 <- sapply(strsplit(Nullauto_rep1$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep1 <- simAuto_rep1[!is.na(simAuto_rep1)]

#replicate2
Nullauto_rep2 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep2.sync.fet',header=FALSE, sep='\t')

Nullauto_rep2$V6 <- as.character(Nullauto_rep2$V6)
simAuto_rep2 <- sapply(strsplit(Nullauto_rep2$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep2 <- simAuto_rep2[!is.na(simAuto_rep2)]

#replicate3
Nullauto_rep3 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep3.sync.fet',header=FALSE, sep='\t')

Nullauto_rep3$V6 <- as.character(Nullauto_rep3$V6)
simAuto_rep3 <- sapply(strsplit(Nullauto_rep3$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep3 <- simAuto_rep3[!is.na(simAuto_rep3)]

#replicate4
Nullauto_rep4 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep4.sync.fet',header=FALSE, sep='\t')

Nullauto_rep4$V6 <- as.character(Nullauto_rep4$V6)
simAuto_rep4 <- sapply(strsplit(Nullauto_rep4$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep4 <- simAuto_rep4[!is.na(simAuto_rep4)]

#replicate5
Nullauto_rep5 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep5.sync.fet',header=FALSE, sep='\t')

Nullauto_rep5$V6 <- as.character(Nullauto_rep5$V6)
simAuto_rep5 <- sapply(strsplit(Nullauto_rep5$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep5 <- simAuto_rep5[!is.na(simAuto_rep5)]

#replicate6
Nullauto_rep6 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep6.sync.fet',header=FALSE, sep='\t')

Nullauto_rep6$V6 <- as.character(Nullauto_rep6$V6)
simAuto_rep6 <- sapply(strsplit(Nullauto_rep6$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep6 <- simAuto_rep6[!is.na(simAuto_rep6)]

#replicate7
Nullauto_rep7 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep7.sync.fet',header=FALSE, sep='\t')

Nullauto_rep7$V6 <- as.character(Nullauto_rep7$V6)
simAuto_rep7 <- sapply(strsplit(Nullauto_rep7$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep7 <- simAuto_rep7[!is.na(simAuto_rep7)]

#replicate8
Nullauto_rep8 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep8.sync.fet',header=FALSE, sep='\t')

Nullauto_rep8$V6 <- as.character(Nullauto_rep8$V6)
simAuto_rep8 <- sapply(strsplit(Nullauto_rep8$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep8 <- simAuto_rep8[!is.na(simAuto_rep8)]

#replicate9
Nullauto_rep9 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep9.sync.fet',header=FALSE, sep='\t')

Nullauto_rep9$V6 <- as.character(Nullauto_rep9$V6)
simAuto_rep9 <- sapply(strsplit(Nullauto_rep9$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep9 <- simAuto_rep9[!is.na(simAuto_rep9)]

#replicate10
Nullauto_rep10 <- read.table('Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep10.sync.fet',header=FALSE, sep='\t')

Nullauto_rep10$V6 <- as.character(Nullauto_rep10$V6)
simAuto_rep10 <- sapply(strsplit(Nullauto_rep10$V6,'='), function(x) as.numeric(x[2]))
simAuto_rep10 <- simAuto_rep10[!is.na(simAuto_rep10)]

##############################################################################
#read the -log10(pvalue) of Fisher's exact test on simulated data X chromosome
##############################################################################
#replicate1
NullX_rep1 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep1.sync.fet',header=FALSE, sep='\t')

NullX_rep1$V6 <- as.character(NullX_rep1$V6)
simX_rep1 <- sapply(strsplit(NullX_rep1$V6,'='), function(x) as.numeric(x[2]))
simX_rep1 <- simX_rep1[!is.na(simX_rep1)]

#replicate2
NullX_rep2 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep2.sync.fet',header=FALSE, sep='\t')

NullX_rep2$V6 <- as.character(NullX_rep2$V6)
simX_rep2 <- sapply(strsplit(NullX_rep2$V6,'='), function(x) as.numeric(x[2]))
simX_rep2 <- simX_rep2[!is.na(simX_rep2)]

#replicate3
NullX_rep3 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep3.sync.fet',header=FALSE, sep='\t')

NullX_rep3$V6 <- as.character(NullX_rep3$V6)
simX_rep3 <- sapply(strsplit(NullX_rep3$V6,'='), function(x) as.numeric(x[2]))
simX_rep3 <- simX_rep3[!is.na(simX_rep3)]

#replicate4
NullX_rep4 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep4.sync.fet',header=FALSE, sep='\t')

NullX_rep4$V6 <- as.character(NullX_rep4$V6)
simX_rep4 <- sapply(strsplit(NullX_rep4$V6,'='), function(x) as.numeric(x[2]))
simX_rep4 <- simX_rep4[!is.na(simX_rep4)]

#replicate5
NullX_rep5 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep5.sync.fet',header=FALSE, sep='\t')

NullX_rep5$V6 <- as.character(NullX_rep5$V6)
simX_rep5 <- sapply(strsplit(NullX_rep5$V6,'='), function(x) as.numeric(x[2]))
simX_rep5 <- simX_rep5[!is.na(simX_rep5)]

#replicate6
NullX_rep6 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep6.sync.fet',header=FALSE, sep='\t')

NullX_rep6$V6 <- as.character(NullX_rep6$V6)
simX_rep6 <- sapply(strsplit(NullX_rep6$V6,'='), function(x) as.numeric(x[2]))
simX_rep6 <- simX_rep6[!is.na(simX_rep6)]

#replicate7
NullX_rep7 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep7.sync.fet',header=FALSE, sep='\t')

NullX_rep7$V6 <- as.character(NullX_rep7$V6)
simX_rep7 <- sapply(strsplit(NullX_rep7$V6,'='), function(x) as.numeric(x[2]))
simX_rep7 <- simX_rep7[!is.na(simX_rep7)]

#replicate8
NullX_rep8 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep8.sync.fet',header=FALSE, sep='\t')

NullX_rep8$V6 <- as.character(NullX_rep8$V6)
simX_rep8 <- sapply(strsplit(NullX_rep8$V6,'='), function(x) as.numeric(x[2]))
simX_rep8 <- simX_rep8[!is.na(simX_rep8)]

#replicate9
NullX_rep9 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep9.sync.fet',header=FALSE, sep='\t')

NullX_rep9$V6 <- as.character(NullX_rep9$V6)
simX_rep9 <- sapply(strsplit(NullX_rep9$V6,'='), function(x) as.numeric(x[2]))
simX_rep9 <- simX_rep9[!is.na(simX_rep9)]

#replicate10
NullX_rep10 <- read.table('Dsim_Ne262_X_557140independentSNPs_simulated_Rep10.sync.fet',header=FALSE, sep='\t')

NullX_rep10$V6 <- as.character(NullX_rep10$V6)
simX_rep10 <- sapply(strsplit(NullX_rep10$V6,'='), function(x) as.numeric(x[2]))
simX_rep10 <- simX_rep10[!is.na(simX_rep10)]

#############################################################################################
#read the -log10(pvalue) of Fisher's exact test on empirical data autosomes and X chromosome
#############################################################################################
#replicate1
rep1 <- read.table('F0F60_Rep1.sync.fet',header=FALSE, sep='\t')

auto_rep1 <- as.character(rep1$V6[rep1$V1 != 'X'])
Auto_rep1 <- sapply(strsplit(auto_rep1,'='), function(x) as.numeric(x[2]))
Auto_rep1 <- Auto_rep1[!is.na(Auto_rep1)]

x_rep1 <- as.character(rep1$V6[rep1$V1 == 'X'])
X_rep1 <- sapply(strsplit(x_rep1,'='), function(x) as.numeric(x[2]))
X_rep1 <- X_rep1[!is.na(X_rep1)]

#replicate2
rep2 <- read.table('F0F60_Rep2.sync.fet',header=FALSE, sep='\t')

auto_rep2 <- as.character(rep2$V6[rep2$V1 != 'X'])
Auto_rep2 <- sapply(strsplit(auto_rep2,'='), function(x) as.numeric(x[2]))
Auto_rep2 <- Auto_rep2[!is.na(Auto_rep2)]

x_rep2 <- as.character(rep2$V6[rep2$V1 == 'X'])
X_rep2 <- sapply(strsplit(x_rep2,'='), function(x) as.numeric(x[2]))
X_rep2 <- X_rep2[!is.na(X_rep2)]

#replicate3
rep3 <- read.table('F0F60_Rep3.sync.fet',header=FALSE, sep='\t')

auto_rep3 <- as.character(rep3$V6[rep3$V1 != 'X'])
Auto_rep3 <- sapply(strsplit(auto_rep3,'='), function(x) as.numeric(x[2]))
Auto_rep3 <- Auto_rep3[!is.na(Auto_rep3)]

x_rep3 <- as.character(rep3$V6[rep3$V1 == 'X'])
X_rep3 <- sapply(strsplit(x_rep3,'='), function(x) as.numeric(x[2]))
X_rep3 <- X_rep3[!is.na(X_rep3)]

#replicate4
rep4 <- read.table('F0F60_Rep4.sync.fet',header=FALSE, sep='\t')

auto_rep4 <- as.character(rep4$V6[rep4$V1 != 'X'])
Auto_rep4 <- sapply(strsplit(auto_rep4,'='), function(x) as.numeric(x[2]))
Auto_rep4 <- Auto_rep4[!is.na(Auto_rep4)]

x_rep4 <- as.character(rep4$V6[rep4$V1 == 'X'])
X_rep4 <- sapply(strsplit(x_rep4,'='), function(x) as.numeric(x[2]))
X_rep4 <- X_rep4[!is.na(X_rep4)]

#replicate5
rep5 <- read.table('F0F60_Rep5.sync.fet',header=FALSE, sep='\t')

auto_rep5 <- as.character(rep5$V6[rep5$V1 != 'X'])
Auto_rep5 <- sapply(strsplit(auto_rep5,'='), function(x) as.numeric(x[2]))
Auto_rep5 <- Auto_rep5[!is.na(Auto_rep5)]

x_rep5 <- as.character(rep5$V6[rep5$V1 == 'X'])
X_rep5 <- sapply(strsplit(x_rep5,'='), function(x) as.numeric(x[2]))
X_rep5 <- X_rep5[!is.na(X_rep5)]

#replicate6
rep6 <- read.table('F0F60_Rep6.sync.fet',header=FALSE, sep='\t')

auto_rep6 <- as.character(rep6$V6[rep6$V1 != 'X'])
Auto_rep6 <- sapply(strsplit(auto_rep6,'='), function(x) as.numeric(x[2]))
Auto_rep6 <- Auto_rep6[!is.na(Auto_rep6)]

x_rep6 <- as.character(rep6$V6[rep6$V1 == 'X'])
X_rep6 <- sapply(strsplit(x_rep6,'='), function(x) as.numeric(x[2]))
X_rep6 <- X_rep6[!is.na(X_rep6)]

#replicate7
rep7 <- read.table('F0F60_Rep7.sync.fet',header=FALSE, sep='\t')

auto_rep7 <- as.character(rep7$V6[rep7$V1 != 'X'])
Auto_rep7 <- sapply(strsplit(auto_rep7,'='), function(x) as.numeric(x[2]))
Auto_rep7 <- Auto_rep7[!is.na(Auto_rep7)]

x_rep7 <- as.character(rep7$V6[rep7$V1 == 'X'])
X_rep7 <- sapply(strsplit(x_rep7,'='), function(x) as.numeric(x[2]))
X_rep7 <- X_rep7[!is.na(X_rep7)]

#replicate8
rep8 <- read.table('F0F60_Rep8.sync.fet',header=FALSE, sep='\t')

auto_rep8 <- as.character(rep8$V6[rep8$V1 != 'X'])
Auto_rep8 <- sapply(strsplit(auto_rep8,'='), function(x) as.numeric(x[2]))
Auto_rep8 <- Auto_rep8[!is.na(Auto_rep8)]

x_rep8 <- as.character(rep8$V6[rep8$V1 == 'X'])
X_rep8 <- sapply(strsplit(x_rep8,'='), function(x) as.numeric(x[2]))
X_rep8 <- X_rep8[!is.na(X_rep8)]

#replicate9
rep9 <- read.table('F0F60_Rep9.sync.fet',header=FALSE, sep='\t')

auto_rep9 <- as.character(rep9$V6[rep9$V1 != 'X'])
Auto_rep9 <- sapply(strsplit(auto_rep9,'='), function(x) as.numeric(x[2]))
Auto_rep9 <- Auto_rep9[!is.na(Auto_rep9)]

x_rep9 <- as.character(rep9$V6[rep9$V1 == 'X'])
X_rep9 <- sapply(strsplit(x_rep9,'='), function(x) as.numeric(x[2]))
X_rep9 <- X_rep9[!is.na(X_rep9)]

#replicate10
rep10 <- read.table('F0F60_Rep10.sync.fet',header=FALSE, sep='\t')

auto_rep10 <- as.character(rep10$V6[rep10$V1 != 'X'])
Auto_rep10 <- sapply(strsplit(auto_rep10,'='), function(x) as.numeric(x[2]))
Auto_rep10 <- Auto_rep10[!is.na(Auto_rep10)]

x_rep10 <- as.character(rep10$V6[rep10$V1 == 'X'])
X_rep10 <- sapply(strsplit(x_rep10,'='), function(x) as.numeric(x[2]))
X_rep10 <- X_rep10[!is.na(X_rep10)]

###############################
#compute FDR corrected q-values
###############################

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

#compute 5% FDR for autosomes and chromosome X
#SNPs with -log10(p-value) values above this threshold have allele frequency change that is more than expected under drift

FDR_auto_rep1 <- compute_FDR(Auto_rep1,simAuto_rep1,0.05)
FDR_x_rep1 <- compute_FDR(X_rep1,simX_rep1,0.05)

FDR_auto_rep2 <- compute_FDR(Auto_rep2,simAuto_rep2,0.05)
FDR_x_rep2 <- compute_FDR(X_rep2,simX_rep2,0.05)

FDR_auto_rep3 <- compute_FDR(Auto_rep3,simAuto_rep3,0.05)
FDR_x_rep3 <- compute_FDR(X_rep3,simX_rep3,0.05)

FDR_auto_rep4 <- compute_FDR(Auto_rep4,simAuto_rep4,0.05)
FDR_x_rep4 <- compute_FDR(X_rep4,simX_rep4,0.05)

FDR_auto_rep5 <- compute_FDR(Auto_rep5,simAuto_rep5,0.05)
FDR_x_rep5 <- compute_FDR(X_rep5,simX_rep5,0.05)

FDR_auto_rep6 <- compute_FDR(Auto_rep6,simAuto_rep6,0.05)
FDR_x_rep6 <- compute_FDR(X_rep6,simX_rep6,0.05)

FDR_auto_rep7 <- compute_FDR(Auto_rep7,simAuto_rep7,0.05)
FDR_x_rep7 <- compute_FDR(X_rep7,simX_rep7,0.05)

FDR_auto_rep8 <- compute_FDR(Auto_rep8,simAuto_rep8,0.05)
FDR_x_rep8 <- compute_FDR(X_rep8,simX_rep8,0.05)

FDR_auto_rep9 <- compute_FDR(Auto_rep9,simAuto_rep9,0.05)
FDR_x_rep9 <- compute_FDR(X_rep9,simX_rep9,0.05)

FDR_auto_rep60 <- compute_FDR(Auto_rep60,simAuto_rep60,0.05)
FDR_x_rep60 <- compute_FDR(X_rep60,simX_rep60,0.05)
