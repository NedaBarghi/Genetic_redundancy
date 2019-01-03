library(Nest)

#Read allele frequency and coverage of SNPs
initial <- read.table('F0F60SNP_AFrising_cov.txt',header=TRUE, sep='\t',nrows=100)
classes <- sapply(initial,class)
SNP_freqCov <- read.table('F0F60SNP_AFrising_cov.txt',header=TRUE, sep='\t',colClasses=classes)

pop_census <- 1000

#compute Ne for each replicate separately in windows of 1000 SNPs
F0F60Rep1 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R1,pt = SNP_freqCov$F60R1, cov0=SNP_freqCov$F0R1_cov,covt=SNP_freqCov$F60R1_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep2 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R2,pt = SNP_freqCov$F60R2, cov0=SNP_freqCov$F0R2_cov,covt=SNP_freqCov$F60R2_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep3 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R3,pt = SNP_freqCov$F60R3, cov0=SNP_freqCov$F0R3_cov,covt=SNP_freqCov$F60R3_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep4 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R4,pt = SNP_freqCov$F60R4, cov0=SNP_freqCov$F0R4_cov,covt=SNP_freqCov$F60R4_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep5 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R5,pt = SNP_freqCov$F60R5, cov0=SNP_freqCov$F0R5_cov,covt=SNP_freqCov$F60R5_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep6 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R6,pt = SNP_freqCov$F60R6, cov0=SNP_freqCov$F0R6_cov,covt=SNP_freqCov$F60R6_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep7 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R7,pt = SNP_freqCov$F60R7, cov0=SNP_freqCov$F0R7_cov,covt=SNP_freqCov$F60R7_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep8 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R8,pt = SNP_freqCov$F60R8, cov0=SNP_freqCov$F0R8_cov,covt=SNP_freqCov$F60R8_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep9 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R9,pt = SNP_freqCov$F60R9, cov0=SNP_freqCov$F0R9_cov,covt=SNP_freqCov$F60R9_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)
F0F60Rep10 <- estimateWndNe(chr=as.character(SNP_freqCov$chr), pos=SNP_freqCov$pos, wndSize=1000,p0 = SNP_freqCov$F0R10,pt = SNP_freqCov$F60R10, cov0=SNP_freqCov$F0R10_cov,covt=SNP_freqCov$F60R10_cov,t=60,unit = 'SNP',ploidy=2,truncAF=NA,method=c("P.planI"), poolSize=rep(pop_census, times=2),Ncensus=pop_census)

#median of autosomes for all 10 replicates
all_Ne_auto <- c(median(F0F60Rep1$Np.planI[F0F60Rep1$chr != 'X']),median(F0F60Rep2$Np.planI[F0F60Rep2$chr != 'X']),median(F0F60Rep3$Np.planI[F0F60Rep3$chr != 'X']),median(F0F60Rep4$Np.planI[F0F60Rep4$chr != 'X']),median(F0F60Rep5$Np.planI[F0F60Rep5$chr != 'X']),
            median(F0F60Rep6$Np.planI[F0F60Rep6$chr != 'X']),median(F0F60Rep7$Np.planI[F0F60Rep7$chr != 'X']),median(F0F60Rep8$Np.planI[F0F60Rep8$chr != 'X']),median(F0F60Rep9$Np.planI[F0F60Rep9$chr != 'X']),median(F0F60Rep10$Np.planI[F0F60Rep10$chr != 'X']))
all_Ne_auto <- round(all_Ne_auto)

#mean of median of autosomes for all 10 replicates
F0F60_Np.planI_woX <- mean(c(median(F0F60Rep1$Np.planI[F0F60Rep1$chr != 'X']),median(F0F60Rep2$Np.planI[F0F60Rep2$chr != 'X']),median(F0F60Rep3$Np.planI[F0F60Rep3$chr != 'X']),median(F0F60Rep4$Np.planI[F0F60Rep4$chr != 'X']),median(F0F60Rep5$Np.planI[F0F60Rep5$chr != 'X']),
                             median(F0F60Rep6$Np.planI[F0F60Rep6$chr != 'X']),median(F0F60Rep7$Np.planI[F0F60Rep7$chr != 'X']),median(F0F60Rep8$Np.planI[F0F60Rep8$chr != 'X']),median(F0F60Rep9$Np.planI[F0F60Rep9$chr != 'X']),median(F0F60Rep10$Np.planI[F0F60Rep10$chr != 'X'])))

#median of chromosome X for all 10 replicates
all_Ne_x <- c(median(F0F60Rep1$Np.planI[F0F60Rep1$chr == 'X']),median(F0F60Rep2$Np.planI[F0F60Rep2$chr == 'X']),median(F0F60Rep3$Np.planI[F0F60Rep3$chr == 'X']),median(F0F60Rep4$Np.planI[F0F60Rep4$chr == 'X']),median(F0F60Rep5$Np.planI[F0F60Rep5$chr == 'X']),
             median(F0F60Rep6$Np.planI[F0F60Rep6$chr == 'X']),median(F0F60Rep7$Np.planI[F0F60Rep7$chr == 'X']),median(F0F60Rep8$Np.planI[F0F60Rep8$chr == 'X']),median(F0F60Rep9$Np.planI[F0F60Rep9$chr == 'X']),median(F0F60Rep10$Np.planI[F0F60Rep10$chr =='X']))
all_Ne_x <- round(all_Ne_x)

#mean of median of chromosome X for all 10 replicates
F0F60_Np.planI_X <- mean(c(median(F0F60Rep1$Np.planI[F0F60Rep1$chr == 'X']),median(F0F60Rep2$Np.planI[F0F60Rep2$chr == 'X']),median(F0F60Rep3$Np.planI[F0F60Rep3$chr == 'X']),median(F0F60Rep4$Np.planI[F0F60Rep4$chr == 'X']),median(F0F60Rep5$Np.planI[F0F60Rep5$chr == 'X']),
                           median(F0F60Rep6$Np.planI[F0F60Rep6$chr == 'X']),median(F0F60Rep7$Np.planI[F0F60Rep7$chr == 'X']),median(F0F60Rep8$Np.planI[F0F60Rep8$chr == 'X']),median(F0F60Rep9$Np.planI[F0F60Rep9$chr == 'X']),median(F0F60Rep10$Np.planI[F0F60Rep10$chr =='X'])))

