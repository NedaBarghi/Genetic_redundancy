###############################
#FDR correction for CMH test#
##############################

#1. Run this script in Rstudio.
#For this step you need F0F60SNP_AFrising_cov.txt file from estimate_Ne.sh to get the starting frequency of SNPs and also estimated Ne to be used for drift simulations. This script simulates drift using the starting frequency of alleles and estimated Ne for autosomes and chromosome X separately.

neutral_simulation_CMH.R

#2. Run CMH test on simulated data for autosomes and chromosome X

perl /popoolation2/cmh-test.pl --input Dsim_Ne291_autosome_4539060independentSNPs_10reps_simulated.sync --output Dsim_Ne291_autosome_4539060independentSNPs_10reps_simulated.sync.cmh --min-count 10 --min-coverage 30 --max-coverage 423 --population 1-2,3-4,5-6,7-8,9-10,11-12,13-14,15-16,17-18,19-20 --remove-temp

perl /popoolation2/cmh-test.pl --input Dsim_Ne262_X_557140independentSNPs_10reps_simulated.sync --output Dsim_Ne262_X_557140independentSNPs_10reps_simulated.sync.cmh --min-count 10 --min-coverage 30 --max-coverage 423 --population 1-2,3-4,5-6,7-8,9-10,11-12,13-14,15-16,17-18,19-20 --remove-temp

3. Compute false discovery rate corrected q-values of the CMH test. At this step you need the result of CMH test on empirical data from SNP_Calling_CMH_FET.sh step 6b. Run this script in Rstudio.

FDR_CMH.R

########################################
#FDR correction for Fisher's exact test#
#########################################

1. Run this script in Rstudio.
#For this step you need F0F60SNP_AFrising_cov.txt file from estimate_Ne.sh to get the starting frequency of SNPs and also estimated Ne to be used for drift simulations. This script simulates drift using the starting frequency of alleles and estimated Ne for autosomes and chromosome X separately.

neutral_simulation_FET.R

#2. Run Fisher's exact test on simulated data for autosomes and chromosome X

#Autosomes
perl /popoolation2/fisher-test.pl --input Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep1.sync --output Dsim_Ne381_autosome_4539060independentSNPs_simulated_rep1.sync.fet --min-count 5 --min-coverage 16 --max-coverage 270

perl /popoolation2/fisher-test.pl --input Dsim_Ne336_autosome_4539060independentSNPs_simulated_rep2.sync --output Dsim_Ne336_autosome_4539060independentSNPs_simulated_rep2.sync.fet --min-count 5 --min-coverage 13 --max-coverage 243

perl /popoolation2/fisher-test.pl --input Dsim_Ne247_autosome_4539060independentSNPs_simulated_rep3.sync --output Dsim_Ne247_autosome_4539060independentSNPs_simulated_rep3.sync.fet --min-count 5 --min-coverage 16 --max-coverage 359

perl /popoolation2/fisher-test.pl --input Dsim_Ne255_autosome_4539060independentSNPs_simulated_rep4.sync --output Dsim_Ne255_autosome_4539060independentSNPs_simulated_rep4.sync.fet --min-count 5 --min-coverage 12 --max-coverage 253

perl /popoolation2/fisher-test.pl --input Dsim_Ne246_autosome_4539060independentSNPs_simulated_rep5.sync --output Dsim_Ne246_autosome_4539060independentSNPs_simulated_rep5.sync.fet --min-count 5 --min-coverage 19 --max-coverage 423

perl /popoolation2/fisher-test.pl --input Dsim_Ne307_autosome_4539060independentSNPs_simulated_rep6.sync --output Dsim_Ne307_autosome_4539060independentSNPs_simulated_rep6.sync.fet --min-count 5 --min-coverage 13 --max-coverage 271

perl /popoolation2/fisher-test.pl --input Dsim_Ne297_autosome_4539060independentSNPs_simulated_rep7.sync --output Dsim_Ne297_autosome_4539060independentSNPs_simulated_rep7.sync.fet --min-count 5 --min-coverage 14 --max-coverage 260

perl /popoolation2/fisher-test.pl --input Dsim_Ne244_autosome_4539060independentSNPs_simulated_rep8.sync --output Dsim_Ne244_autosome_4539060independentSNPs_simulated_rep8.sync.fet --min-count 5 --min-coverage 13 --max-coverage 226

perl /popoolation2/fisher-test.pl --input Dsim_Ne310_autosome_4539060independentSNPs_simulated_rep9.sync --output Dsim_Ne310_autosome_4539060independentSNPs_simulated_rep9.sync.fet --min-count 5 --min-coverage 19 --max-coverage 383

perl /popoolation2/fisher-test.pl --input Dsim_Ne287_autosome_4539060independentSNPs_simulated_rep10.sync --output Dsim_Ne287_autosome_4539060independentSNPs_simulated_rep10.sync.fet --min-count 5 --min-coverage 14 --max-coverage 240

#chromosome X
perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep1.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep1.sync.fet --min-count 5 --min-coverage 16 --max-coverage 270

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep2.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep2.sync.fet --min-count 5 --min-coverage 13 --max-coverage 243

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep3.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep3.sync.fet --min-count 5 --min-coverage 16 --max-coverage 359

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep4.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep4.sync.fet --min-count 5 --min-coverage 12 --max-coverage 253

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep5.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep5.sync.fet --min-count 5 --min-coverage 19 --max-coverage 423

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep6.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep6.sync.fet --min-count 5 --min-coverage 13 --max-coverage 271

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep7.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep7.sync.fet --min-count 5 --min-coverage 14 --max-coverage 260

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep8.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep8.sync.fet --min-count 5 --min-coverage 13 --max-coverage 226

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep9.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep9.sync.fet --min-count 5 --min-coverage 19 --max-coverage 383

perl /popoolation2/fisher-test.pl --input Dsim_Ne262_X_557140independentSNPs_simulated_Rep10.sync --output Dsim_Ne262_X_557140independentSNPs_simulated_Rep10.sync.fet --min-count 5 --min-coverage 14 --max-coverage 240

#3. Compute false discovery rate corrected q-values of the Fisher exact tests. At this step you need the result of Fisher exact tests on empirical data from SNP_Calling_CMH_FET.sh step 7b. Run this script in Rstudio.

FDR_FET.R
