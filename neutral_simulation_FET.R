#Read allele frequency of SNPs
initial <- read.table('F0F60SNP_AFrising_cov.txt',header=TRUE, sep='\t',nrows=100)
classes <- sapply(initial,class)
SNP_freqCov <- read.table('F0F60SNP_AFrising_cov.txt',header=TRUE, sep='\t',colClasses=classes)

# simulate drift trajectories
fast.trajectory <- function(p0, Ne, t, s=0, h=0.5) {
  traj <- matrix(NA, ncol=length(t), nrow=max(length(p0), length(s), length(h)), dimnames=list(c(), paste0("F", t)))
  if(0 %in% t)
    traj[,"F0"] <- p0
  
  g <- 1
  p <- p0
  q <- 1-p0
  wAA <- 1+s
  wAa <- 1+h*s
  waa <- 1
  # simulate allele frequencies across time
  while(g <= max(t)) {
    # apply selection and random drift
    p <- (wAA*p^2 + wAa*p*q) / (wAA*p^2 + wAa*2*p*q + waa*q^2)
    if(!is.na(Ne))
      p <- rbinom(length(p), Ne, p)/Ne
    q <- 1-p
    
    # if necessary then save current allele frequency to results
    if(g %in% t)
      traj[,paste0("F", g)] <- p
    
    g <- g+1
  }
  
  if(nrow(traj) == 1)
    return(as.vector(traj))
  
  return(traj)
}

#############################
#simulate drift for autosomes
#############################
gen <- 60
allreps_Ne <- c(381,336,247,255,246,307,397,244,310,287) #these estimated Ne are all_Ne_auto in Ne_estimation.R

system.time(simTraj_1 <- fast.trajectory(p0=SNP_freqCov$F0R1[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[1]*2, t=c(0,gen)))
system.time(simTraj_2 <- fast.trajectory(p0=SNP_freqCov$F0R2[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[2]*2, t=c(0,gen)))
system.time(simTraj_3 <- fast.trajectory(p0=SNP_freqCov$F0R3[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[3]*2, t=c(0,gen)))
system.time(simTraj_4 <- fast.trajectory(p0=SNP_freqCov$F0R4[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[4]*2, t=c(0,gen)))
system.time(simTraj_5 <- fast.trajectory(p0=SNP_freqCov$F0R5[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[5]*2, t=c(0,gen)))
system.time(simTraj_6 <- fast.trajectory(p0=SNP_freqCov$F0R6[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[6]*2, t=c(0,gen)))
system.time(simTraj_7 <- fast.trajectory(p0=SNP_freqCov$F0R7[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[7]*2, t=c(0,gen)))
system.time(simTraj_8 <- fast.trajectory(p0=SNP_freqCov$F0R8[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[8]*2, t=c(0,gen)))
system.time(simTraj_9 <- fast.trajectory(p0=SNP_freqCov$F0R9[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[9]*2, t=c(0,gen)))
system.time(simTraj_10 <- fast.trajectory(p0=SNP_freqCov$F0R10[SNP_freqCov$chr != 'X'], Ne=allreps_Ne[10]*2, t=c(0,gen)))

# sample alleles at certain sequence coverage for all 10 replicates
Bcov <- c(193,172,255,178,300,186,184,156,274,169) #average genomic coverage of 10 hot-evolved replicates at F0
Ecov <- c(128,88,74,61,84,75,103,99,115,114) #average genomic coverage of 10 hot-evolved replicates at F60

covBase_1 <- rpois(nrow(simTraj_1), Bcov[1])
cntBase_1 <- rbinom(nrow(simTraj_1), size=covBase_1, prob=simTraj_1[,1])
covEvol_1 <- rpois(nrow(simTraj_1), Ecov[1])
cntEvol_1 <- rbinom(nrow(simTraj_1), size=covEvol_1, prob=simTraj_1[,2])

covBase_2 <- rpois(nrow(simTraj_2), Bcov[2])
cntBase_2 <- rbinom(nrow(simTraj_2), size=covBase_2, prob=simTraj_2[,1])
covEvol_2 <- rpois(nrow(simTraj_2), Ecov[2])
cntEvol_2 <- rbinom(nrow(simTraj_2), size=covEvol_2, prob=simTraj_2[,2])

covBase_3 <- rpois(nrow(simTraj_3), Bcov[3])
cntBase_3 <- rbinom(nrow(simTraj_3), size=covBase_3, prob=simTraj_3[,1])
covEvol_3 <- rpois(nrow(simTraj_3), Ecov[3])
cntEvol_3 <- rbinom(nrow(simTraj_3), size=covEvol_3, prob=simTraj_3[,2])

covBase_4 <- rpois(nrow(simTraj_4), Bcov[4])
cntBase_4 <- rbinom(nrow(simTraj_4), size=covBase_4, prob=simTraj_4[,1])
covEvol_4 <- rpois(nrow(simTraj_4), Ecov[4])
cntEvol_4 <- rbinom(nrow(simTraj_4), size=covEvol_4, prob=simTraj_4[,2])

covBase_5 <- rpois(nrow(simTraj_5), Bcov[5])
cntBase_5 <- rbinom(nrow(simTraj_5), size=covBase_5, prob=simTraj_5[,1])
covEvol_5 <- rpois(nrow(simTraj_5), Ecov[5])
cntEvol_5 <- rbinom(nrow(simTraj_5), size=covEvol_5, prob=simTraj_5[,2])

covBase_6 <- rpois(nrow(simTraj_6), Bcov[6])
cntBase_6 <- rbinom(nrow(simTraj_6), size=covBase_6, prob=simTraj_6[,1])
covEvol_6 <- rpois(nrow(simTraj_6), Ecov[6])
cntEvol_6 <- rbinom(nrow(simTraj_6), size=covEvol_6, prob=simTraj_6[,2])

covBase_7 <- rpois(nrow(simTraj_7), Bcov[7])
cntBase_7 <- rbinom(nrow(simTraj_7), size=covBase_7, prob=simTraj_7[,1])
covEvol_7 <- rpois(nrow(simTraj_7), Ecov[7])
cntEvol_7 <- rbinom(nrow(simTraj_7), size=covEvol_7, prob=simTraj_7[,2])

covBase_8 <- rpois(nrow(simTraj_8), Bcov[8])
cntBase_8 <- rbinom(nrow(simTraj_8), size=covBase_8, prob=simTraj_8[,1])
covEvol_8 <- rpois(nrow(simTraj_8), Ecov[8])
cntEvol_8 <- rbinom(nrow(simTraj_8), size=covEvol_8, prob=simTraj_8[,2])

covBase_9 <- rpois(nrow(simTraj_9), Bcov[9])
cntBase_9 <- rbinom(nrow(simTraj_9), size=covBase_9, prob=simTraj_9[,1])
covEvol_9 <- rpois(nrow(simTraj_9), Ecov[9])
cntEvol_9 <- rbinom(nrow(simTraj_9), size=covEvol_9, prob=simTraj_9[,2])

covBase_10 <- rpois(nrow(simTraj_10), Bcov[10])
cntBase_10 <- rbinom(nrow(simTraj_10), size=covBase_10, prob=simTraj_10[,1])
covEvol_10 <- rpois(nrow(simTraj_10), Ecov[10])
cntEvol_10 <- rbinom(nrow(simTraj_10), size=covEvol_10, prob=simTraj_10[,2])

#concatenate counts and coverages for all 10 replicates
syncTable <- c()
ListCnt <- list(cntBase_1,cntEvol_1, cntBase_2, cntEvol_2, cntBase_3, cntEvol_3,cntBase_4, cntEvol_4,cntBase_5, cntEvol_5,cntBase_6, cntEvol_6,cntBase_7, cntEvol_7,cntBase_8, cntEvol_8,cntBase_9, cntEvol_9,cntBase_10, cntEvol_10)
ListCov <- list(covBase_1,covEvol_1, covBase_2, covEvol_2, covBase_3, covEvol_3,covBase_4, covEvol_4,covBase_5, covEvol_5,covBase_6, covEvol_6,covBase_7, covEvol_7,covBase_8, covEvol_8,covBase_9, covEvol_9,covBase_10, covEvol_10)

#this function converts the counts of simulated drift to sync format. For simplicity all polymorphic sites are A/C combinations.
f <- function(x){paste(ListCnt[[i]][x], 0,ListCov[[i]][x]-ListCnt[[i]][x] , 0, 0, 0, sep=":")}

snpNum <- length(SNP_freqCov$pos[SNP_freqCov$chr != 'X'])

for (i in 1:length(ListCnt)) {
  syncTable <- cbind(syncTable, sapply(1:snpNum,f))
}
base <- c(rep('A',time=snpNum))

rep1_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,1:2])
rep2_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,3:4])
rep3_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,5:6])
rep4_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,7:8])
rep5_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,9:10])
rep6_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,11:12])
rep7_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,13:14])
rep8_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,15:16])
rep9_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,17:18])
rep10_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr != 'X']),SNP_freqCov$pos[SNP_freqCov$chr != 'X'], base,syncTable[,19:20])

#write the results in a file. Note that the order of samples in each file is F0R*, F60R*
write.table(rep1_syncTable, file='/path_to_save/Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep1.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep2_syncTable, file='/path_to_save/Dsim_Ne336_autosome_4539060independentSNPs_simulated_rep2.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep3_syncTable, file='/path_to_save/Dsim_Ne247_autosome_4539060independentSNPs_simulated_rep3.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep4_syncTable, file='/path_to_save/Dsim_Ne255_autosome_4539060independentSNPs_simulated_rep4.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep5_syncTable, file='/path_to_save/Dsim_Ne246_autosome_4539060independentSNPs_simulated_rep5.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep6_syncTable, file='/path_to_save/Dsim_Ne307_autosome_4539060independentSNPs_simulated_rep6.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep7_syncTable, file='/path_to_save/Dsim_Ne297_autosome_4539060independentSNPs_simulated_rep7.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep8_syncTable, file='/path_to_save/Dsim_Ne244_autosome_4539060independentSNPs_simulated_rep8.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep9_syncTable, file='/path_to_save/Dsim_Ne310_autosome_4539060independentSNPs_simulated_rep9.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep10_syncTable, file='/path_to_save/Dsim_Ne287_autosome_4539060independentSNPs_simulated_rep10.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)

#################################
#simulate drift for chromosome X
#################################
est_Ne <- 262 #this value is F0F60_Np.planI_X in Ne_estimation.R
gen <- 60

# simulate drift trajectories
fast.trajectory <- function(p0, Ne, t, s=0, h=0.5) {
  traj <- matrix(NA, ncol=length(t), nrow=max(length(p0), length(s), length(h)), dimnames=list(c(), paste0("F", t)))
  if(0 %in% t)
    traj[,"F0"] <- p0
  
  g <- 1
  p <- p0
  q <- 1-p0
  wAA <- 1+s
  wAa <- 1+h*s
  waa <- 1
  # simulate allele frequencies across time
  while(g <= max(t)) {
    # apply selection and random drift
    p <- (wAA*p^2 + wAa*p*q) / (wAA*p^2 + wAa*2*p*q + waa*q^2)
    if(!is.na(Ne))
      p <- rbinom(length(p), Ne, p)/Ne
    q <- 1-p
    
    # if necessary then save current allele frequency to results
    if(g %in% t)
      traj[,paste0("F", g)] <- p
    
    g <- g+1
  }
  
  if(nrow(traj) == 1)
    return(as.vector(traj))
  
  return(traj)
}

system.time(simTraj_1 <- fast.trajectory(p0=SNP_freqCov$F0R1[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_2 <- fast.trajectory(p0=SNP_freqCov$F0R2[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_3 <- fast.trajectory(p0=SNP_freqCov$F0R3[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_4 <- fast.trajectory(p0=SNP_freqCov$F0R4[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_5 <- fast.trajectory(p0=SNP_freqCov$F0R5[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_6 <- fast.trajectory(p0=SNP_freqCov$F0R6[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_7 <- fast.trajectory(p0=SNP_freqCov$F0R7[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_8 <- fast.trajectory(p0=SNP_freqCov$F0R8[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_9 <- fast.trajectory(p0=SNP_freqCov$F0R9[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))
system.time(simTraj_10 <- fast.trajectory(p0=SNP_freqCov$F0R10[SNP_freqCov$chr == 'X'], Ne=est_Ne*2*2, t=c(0,gen)))

# sample alleles at certain sequence coverage for all 10 replicates
Bcov <- c(193,172,255,178,300,186,184,156,274,169) #average genomic coverage of 10 hot-evolved replicates at F0
Ecov <- c(128,88,74,61,84,75,103,99,115,114) #average genomic coverage of 10 hot-evolved replicates at F60

covBase_1 <- rpois(nrow(simTraj_1), Bcov[1])
cntBase_1 <- rbinom(nrow(simTraj_1), size=covBase_1, prob=simTraj_1[,1])
covEvol_1 <- rpois(nrow(simTraj_1), Ecov[1])
cntEvol_1 <- rbinom(nrow(simTraj_1), size=covEvol_1, prob=simTraj_1[,2])

covBase_2 <- rpois(nrow(simTraj_2), Bcov[2])
cntBase_2 <- rbinom(nrow(simTraj_2), size=covBase_2, prob=simTraj_2[,1])
covEvol_2 <- rpois(nrow(simTraj_2), Ecov[2])
cntEvol_2 <- rbinom(nrow(simTraj_2), size=covEvol_2, prob=simTraj_2[,2])

covBase_3 <- rpois(nrow(simTraj_3), Bcov[3])
cntBase_3 <- rbinom(nrow(simTraj_3), size=covBase_3, prob=simTraj_3[,1])
covEvol_3 <- rpois(nrow(simTraj_3), Ecov[3])
cntEvol_3 <- rbinom(nrow(simTraj_3), size=covEvol_3, prob=simTraj_3[,2])

covBase_4 <- rpois(nrow(simTraj_4), Bcov[4])
cntBase_4 <- rbinom(nrow(simTraj_4), size=covBase_4, prob=simTraj_4[,1])
covEvol_4 <- rpois(nrow(simTraj_4), Ecov[4])
cntEvol_4 <- rbinom(nrow(simTraj_4), size=covEvol_4, prob=simTraj_4[,2])

covBase_5 <- rpois(nrow(simTraj_5), Bcov[5])
cntBase_5 <- rbinom(nrow(simTraj_5), size=covBase_5, prob=simTraj_5[,1])
covEvol_5 <- rpois(nrow(simTraj_5), Ecov[5])
cntEvol_5 <- rbinom(nrow(simTraj_5), size=covEvol_5, prob=simTraj_5[,2])

covBase_6 <- rpois(nrow(simTraj_6), Bcov[6])
cntBase_6 <- rbinom(nrow(simTraj_6), size=covBase_6, prob=simTraj_6[,1])
covEvol_6 <- rpois(nrow(simTraj_6), Ecov[6])
cntEvol_6 <- rbinom(nrow(simTraj_6), size=covEvol_6, prob=simTraj_6[,2])

covBase_7 <- rpois(nrow(simTraj_7), Bcov[7])
cntBase_7 <- rbinom(nrow(simTraj_7), size=covBase_7, prob=simTraj_7[,1])
covEvol_7 <- rpois(nrow(simTraj_7), Ecov[7])
cntEvol_7 <- rbinom(nrow(simTraj_7), size=covEvol_7, prob=simTraj_7[,2])

covBase_8 <- rpois(nrow(simTraj_8), Bcov[8])
cntBase_8 <- rbinom(nrow(simTraj_8), size=covBase_8, prob=simTraj_8[,1])
covEvol_8 <- rpois(nrow(simTraj_8), Ecov[8])
cntEvol_8 <- rbinom(nrow(simTraj_8), size=covEvol_8, prob=simTraj_8[,2])

covBase_9 <- rpois(nrow(simTraj_9), Bcov[9])
cntBase_9 <- rbinom(nrow(simTraj_9), size=covBase_9, prob=simTraj_9[,1])
covEvol_9 <- rpois(nrow(simTraj_9), Ecov[9])
cntEvol_9 <- rbinom(nrow(simTraj_9), size=covEvol_9, prob=simTraj_9[,2])

covBase_10 <- rpois(nrow(simTraj_10), Bcov[10])
cntBase_10 <- rbinom(nrow(simTraj_10), size=covBase_10, prob=simTraj_10[,1])
covEvol_10 <- rpois(nrow(simTraj_10), Ecov[10])
cntEvol_10 <- rbinom(nrow(simTraj_10), size=covEvol_10, prob=simTraj_10[,2])

#concatenate counts and coverages for all 10 replicates
syncTable <- c()
ListCnt <- list(cntBase_1,cntEvol_1, cntBase_2, cntEvol_2, cntBase_3, cntEvol_3,cntBase_4, cntEvol_4,cntBase_5, cntEvol_5,cntBase_6, cntEvol_6,cntBase_7, cntEvol_7,cntBase_8, cntEvol_8,cntBase_9, cntEvol_9,cntBase_10, cntEvol_10)
ListCov <- list(covBase_1,covEvol_1, covBase_2, covEvol_2, covBase_3, covEvol_3,covBase_4, covEvol_4,covBase_5, covEvol_5,covBase_6, covEvol_6,covBase_7, covEvol_7,covBase_8, covEvol_8,covBase_9, covEvol_9,covBase_10, covEvol_10)

#this function converts the counts of simulated drift to sync format. For simplicity all polymorphic sites are A/C combinations.
f <- function(x){paste(ListCnt[[i]][x], 0,ListCov[[i]][x]-ListCnt[[i]][x] , 0, 0, 0, sep=":")}

snpNum <- length(SNP_freqCov$pos[SNP_freqCov$chr == 'X'])

for (i in 1:length(ListCnt)) {
  syncTable <- cbind(syncTable, sapply(1:snpNum,f))
}
base <- c(rep('A',time=snpNum))

rep1_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,1:2])
rep2_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,3:4])
rep3_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,5:6])
rep4_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,7:8])
rep5_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,9:10])
rep6_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,11:12])
rep7_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,13:14])
rep8_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,15:16])
rep9_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,17:18])
rep10_syncTable <- cbind(as.character(SNP_freqCov$chr[SNP_freqCov$chr == 'X']),SNP_freqCov$pos[SNP_freqCov$chr == 'X'], base,syncTable[,19:20])

#write the results in a file. Note that the order of samples in each file is F0R*, F60R*
write.table(rep1_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep1.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep2_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep2.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep3_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep3.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep4_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep4.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep5_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep5.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep6_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep6.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep7_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep7.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep8_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep8.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep9_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep9.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)
write.table(rep10_syncTable, file='/path_to_save/Dsim_Ne262_X_557140independentSNPs_simulated_rep10.sync', sep = '\t', col.names = FALSE, row.names = FALSE, quote=FALSE)

