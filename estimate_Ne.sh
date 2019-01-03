###############
#Ne estimation#
###############

#1. Compute allele frequencies and coverage of SNPs
#1a. Separate F0 and F60 samples from F0-F60SNP_CMH_FET_blockID.sync file (Dryad Digital Repository: https://doi.org/10.5061/dryad.rr137kn) 
cut -f 1-13,64-73 F0-F60SNP_CMH_FET_blockID.sync > F0F60SNP.sync

python compute_AF_cov_Ne.py --input F0F60SNP.sync --output F0F60SNP_AFrising_cov.txt

#2. Estimate Ne, run this script in Rstudio.
Ne_estimation.R