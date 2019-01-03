#############
#SNP calling#
#############

#1. The raw reads for each samples of founder population (F0) can be downloaded using the SRA run accession numbers provided in S7A Table in the manuscript. Trimming, mapping and filtering of reads was performed as described in Barghi et al. (Barghi N, Tobler R, Nolte V, Schlötterer C. Drosophila simulans : A Species with Improved Resolution in Evolve and Resequence Studies. G3 (Bethesda). 2017;7:2337–43). 

#2. A mpileup file was made for bam files of all 10 libraries of F0 and converted to sync file with minimum quality of 40 (bam file nomenclature is the same as in S7A Table). Samtools version version 1.2 was used. popoolation2 software is from Kofler R, Pandey RV, Schlötterer C. PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics. 2011;27(24):3435–6.

samtools mpileup --max-depth 8000 -6 -B -Q 0 -f reference_genome.fa Dsim_Fl_Base_1.bam Dsim_Fl_Base_2.bam Dsim_Fl_Base_3.bam Dsim_Fl_Base_4.bam Dsim_Fl_Base_5.bam Dsim_Fl_Base_6.bam Dsim_Fl_Base_7.bam Dsim_Fl_Base_8.bam Dsim_Fl_Base_9.bam Dsim_Fl_Base_10.bam | java -Xmx20g -jar /popoolation2_1201/mpileup2sync.jar --input /dev/stdin --output F0_Q40.sync --fastq-type sanger --min-qual 40 --threads 12

#3. SNPs can be called using call_SNP.py

python call_SNP.py --input F0_Q40.sync --output F0_Q40_SNPs.sync

#4. Identify indels and mask indels, TEs and Y-translocated genes
#4a. The raw reads for each samples of founder and evolved populations can be downloaded using the SRA run accession numbers provided in S7A Table in the manuscript. Trimming, mapping and filtering of reads was performed as described in Barghi et al. (Barghi N, Tobler R, Nolte V, Schlötterer C. Drosophila simulans : A Species with Improved Resolution in Evolve and Resequence Studies. G3 (Bethesda). 2017;7:2337–43). Make a mpileup file of all libraries (10 replicates and 7 timepoints):

samtools mpileup --max-depth 8000 -6 -B -Q 0 -f reference_genome.fa Dsim_Fl_Base_1.bam Dsim_Fl_Base_2.bam Dsim_Fl_Base_3.bam Dsim_Fl_Base_4.bam Dsim_Fl_Base_5.bam Dsim_Fl_Base_6.bam Dsim_Fl_Base_7.bam Dsim_Fl_Base_8.bam Dsim_Fl_Base_9.bam Dsim_Fl_Base_10.bam Dsim_Fl_Hot_F10_1.bam Dsim_Fl_Hot_F10_2.bam Dsim_Fl_Hot_F10_3.bam Dsim_Fl_Hot_F10_4.bam Dsim_Fl_Hot_F10_5.bam Dsim_Fl_Hot_F10_6.bam Dsim_Fl_Hot_F10_7.bam Dsim_Fl_Hot_F10_8.bam Dsim_Fl_Hot_F10_9.bam Dsim_Fl_Hot_F10_10.bam Dsim_Fl_Hot_F20_1.bam Dsim_Fl_Hot_F20_2.bam Dsim_Fl_Hot_F20_3.bam Dsim_Fl_Hot_F20_4.bam Dsim_Fl_Hot_F20_5.bam Dsim_Fl_Hot_F20_6.bam Dsim_Fl_Hot_F20_7.bam Dsim_Fl_Hot_F20_8.bam Dsim_Fl_Hot_F20_9.bam Dsim_Fl_Hot_F20_10.bam Dsim_Fl_Hot_F30_1.bam Dsim_Fl_Hot_F30_2.bam Dsim_Fl_Hot_F30_3.bam Dsim_Fl_Hot_F30_4.bam Dsim_Fl_Hot_F30_5.bam Dsim_Fl_Hot_F30_6.bam Dsim_Fl_Hot_F30_7.bam Dsim_Fl_Hot_F30_8.bam Dsim_Fl_Hot_F30_9.bam Dsim_Fl_Hot_F30_10.bam Dsim_Fl_Hot_F40_1.bam Dsim_Fl_Hot_F40_2.bam Dsim_Fl_Hot_F40_3.bam Dsim_Fl_Hot_F40_4.bam Dsim_Fl_Hot_F40_5.bam Dsim_Fl_Hot_F40_6.bam Dsim_Fl_Hot_F40_7.bam Dsim_Fl_Hot_F40_8.bam Dsim_Fl_Hot_F40_9.bam Dsim_Fl_Hot_F40_10.bam Dsim_Fl_Hot_F50_1.bam Dsim_Fl_Hot_F50_2.bam Dsim_Fl_Hot_F50_3.bam Dsim_Fl_Hot_F50_4.bam Dsim_Fl_Hot_F50_5.bam Dsim_Fl_Hot_F50_6.bam Dsim_Fl_Hot_F50_7.bam Dsim_Fl_Hot_F50_8.bam Dsim_Fl_Hot_F50_9.bam Dsim_Fl_Hot_F50_10.bam Dsim_Fl_Hot_F60_1.bam Dsim_Fl_Hot_F60_2.bam Dsim_Fl_Hot_F60_3.bam Dsim_Fl_Hot_F60_4.bam Dsim_Fl_Hot_F60_5.bam Dsim_Fl_Hot_F60_6.bam Dsim_Fl_Hot_F60_7.bam Dsim_Fl_Hot_F60_8.bam Dsim_Fl_Hot_F60_9.bam Dsim_Fl_Hot_F60_10.bam > F0-F60.mpileup

#4b. Identify indels
perl /popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --input F0-F60.mpileup --output F0-F60_indelregions.gtf --indel-window 5 --min-count 167


#4c. combine indels, TEs and Y-translocated genes into a file. TEfinal.gff was annotated using the pipeline described in Kofler et al. (Kofler R, Nolte V, Schlötterer C. Tempo and Mode of Transposable Element Activity in Drosophila. Plos Genetics. 2015; 117: e1005406) and Ytransloc_regions_200bp.gff contains 200-bp flanking the SNPs specific to autosomal genes translocated to the Y chromosome in Tobler et al. (Tobler R, Nolte V, Schlötterer C. High rate of translocation-based gene birth on the Drosophila Y chromosome. Proc Natl Acad Sci USA. 2017; 114:44: 11721-11726).

cat F0-F60_indelregions.gtf TEfinal.gff Ytransloc_regions_200bp.gff > F0-F60_indel_TE_Ytransloc.gtf

#4d. mask indels, TEs and Y-translocated genes
perl /popoolation2_1201/indel_filtering/filter-sync-by-gtf.pl --gtf F0-F60_indel_TE_Ytransloc.gtf --input F0_Q40_SNPs.sync --output F0_Q40_SNPs_indel_TE_masked.sync

#5. Extract the SNPs for all samples
#5a. A mpileup file was made for bam files of all libraries (10 replicates and 7 timepoints) and converted to sync file as follows (bam file nomenclature is the same as in S7A Table):

samtools mpileup --max-depth 8000 -6 -B -Q 0 -f reference_genome.fa Dsim_Fl_Base_1.bam Dsim_Fl_Base_2.bam Dsim_Fl_Base_3.bam Dsim_Fl_Base_4.bam Dsim_Fl_Base_5.bam Dsim_Fl_Base_6.bam Dsim_Fl_Base_7.bam Dsim_Fl_Base_8.bam Dsim_Fl_Base_9.bam Dsim_Fl_Base_10.bam Dsim_Fl_Hot_F10_1.bam Dsim_Fl_Hot_F10_2.bam Dsim_Fl_Hot_F10_3.bam Dsim_Fl_Hot_F10_4.bam Dsim_Fl_Hot_F10_5.bam Dsim_Fl_Hot_F10_6.bam Dsim_Fl_Hot_F10_7.bam Dsim_Fl_Hot_F10_8.bam Dsim_Fl_Hot_F10_9.bam Dsim_Fl_Hot_F10_10.bam Dsim_Fl_Hot_F20_1.bam Dsim_Fl_Hot_F20_2.bam Dsim_Fl_Hot_F20_3.bam Dsim_Fl_Hot_F20_4.bam Dsim_Fl_Hot_F20_5.bam Dsim_Fl_Hot_F20_6.bam Dsim_Fl_Hot_F20_7.bam Dsim_Fl_Hot_F20_8.bam Dsim_Fl_Hot_F20_9.bam Dsim_Fl_Hot_F20_10.bam Dsim_Fl_Hot_F30_1.bam Dsim_Fl_Hot_F30_2.bam Dsim_Fl_Hot_F30_3.bam Dsim_Fl_Hot_F30_4.bam Dsim_Fl_Hot_F30_5.bam Dsim_Fl_Hot_F30_6.bam Dsim_Fl_Hot_F30_7.bam Dsim_Fl_Hot_F30_8.bam Dsim_Fl_Hot_F30_9.bam Dsim_Fl_Hot_F30_10.bam Dsim_Fl_Hot_F40_1.bam Dsim_Fl_Hot_F40_2.bam Dsim_Fl_Hot_F40_3.bam Dsim_Fl_Hot_F40_4.bam Dsim_Fl_Hot_F40_5.bam Dsim_Fl_Hot_F40_6.bam Dsim_Fl_Hot_F40_7.bam Dsim_Fl_Hot_F40_8.bam Dsim_Fl_Hot_F40_9.bam Dsim_Fl_Hot_F40_10.bam Dsim_Fl_Hot_F50_1.bam Dsim_Fl_Hot_F50_2.bam Dsim_Fl_Hot_F50_3.bam Dsim_Fl_Hot_F50_4.bam Dsim_Fl_Hot_F50_5.bam Dsim_Fl_Hot_F50_6.bam Dsim_Fl_Hot_F50_7.bam Dsim_Fl_Hot_F50_8.bam Dsim_Fl_Hot_F50_9.bam Dsim_Fl_Hot_F50_10.bam Dsim_Fl_Hot_F60_1.bam Dsim_Fl_Hot_F60_2.bam Dsim_Fl_Hot_F60_3.bam Dsim_Fl_Hot_F60_4.bam Dsim_Fl_Hot_F60_5.bam Dsim_Fl_Hot_F60_6.bam Dsim_Fl_Hot_F60_7.bam Dsim_Fl_Hot_F60_8.bam Dsim_Fl_Hot_F60_9.bam Dsim_Fl_Hot_F60_10.bam | java -Xmx20g -jar  /popoolation2_1201/mpileup2sync.jar --input /dev/stdin --output F0-F60_Q20.sync --fastq-type sanger --min-qual 20 --threads 12

#5b. Extract the called SNPs (at stage 4d) for all samples

awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print f1[$1$2], $0}' F0_Q40_SNPs_indel_TE_masked.sync F0-F60_Q20.sync | cut -f 1-3,17-86 > F0-F60SNP.sync

#Note that F0-F60SNP.sync is still NOT the final set of SNPs. F0-F60SNP.sync will be filtered during CMH test and thus the final set of SNPs (see below)

#############
#CMH test#
#############

#6. Run Cochran-Mantel-Haenszel (CMH) test on F0 and F60 samples
#6a. Separate F0 and F60 samples
cut -f 1-13,64-73 F0-F60SNP.sync > F0F60SNP.sync

#6b. Run CMH test
perl /popoolation2/cmh-test.pl --input F0F60SNP.sync --output F0F60SNP.sync.cmh --min-count 10 --min-coverage 30 --max-coverage 423 --population 1-11,2-12,3-13,4-14,5-15,6-16,7-17,8-18,9-19,10-20 --remove-temp

#F0F60SNP.sync.cmh is already filtered for low and high coverage SNPs and is the final SNP set. 

#############
#FET test#
#############

#7. Run Fisher's exact test (FET) on all 10 hot-evolved replicates 
#7a. Separate F0 and F60 samples for each replicate in a separate file
cut -f 1-4,14  F0F60SNP.sync > F0F60_Rep1.sync

cut -f 1-3,5,15 F0F60SNP.sync > F0F60_Rep2.sync

cut -f 1-3,6,16 F0F60SNP.sync > F0F60_Rep3.sync

cut -f 1-3,7,17 F0F60SNP.sync > F0F60_Rep4.sync

cut -f 1-3,8,18 F0F60SNP.sync > F0F60_Rep5.sync

cut -f 1-3,9,19 F0F60SNP.sync > F0F60_Rep6.sync

cut -f 1-3,10,20 F0F60SNP.sync > F0F60_Rep7.sync

cut -f 1-3,11,21 F0F60SNP.sync > F0F60_Rep8.sync

cut -f 1-3,12,22 F0F60SNP.sync > F0F60_Rep9.sync

cut -f 1-3,13,23 F0F60SNP.sync > F0F60_Rep10.sync

#7b. Run FET
perl /popoolation2/fisher-test.pl --input F0F60_Rep1.sync --output F0F60_Rep1.sync.fet --min-count 5 --min-coverage 16 --max-coverage 270

perl /popoolation2/fisher-test.pl --input F0F60_Rep2.sync --output F0F60_Rep2.sync.fet --min-count 5 --min-coverage 13 --max-coverage 243

perl /popoolation2/fisher-test.pl --input F0F60_Rep3.sync --output F0F60_Rep3.sync.fet --min-count 5 --min-coverage 16 --max-coverage 359

perl /popoolation2/fisher-test.pl --input F0F60_Rep4.sync --output F0F60_Rep4.sync.fet --min-count 5 --min-coverage 12 --max-coverage 253

perl /popoolation2/fisher-test.pl --input F0F60_Rep5.sync --output F0F60_Rep5.sync.fet --min-count 5 --min-coverage 19 --max-coverage 423

perl /popoolation2/fisher-test.pl --input F0F60_Rep6.sync --output F0F60_Rep6.sync.fet --min-count 5 --min-coverage 13 --max-coverage 271

perl /popoolation2/fisher-test.pl --input F0F60_Rep7.sync --output F0F60_Rep7.sync.fet --min-count 5 --min-coverage 14 --max-coverage 260

perl /popoolation2/fisher-test.pl --input F0F60_Rep8.sync --output F0F60_Rep8.sync.fet --min-count 5 --min-coverage 13 --max-coverage 226

perl /popoolation2/fisher-test.pl --input F0F60_Rep9.sync --output F0F60_Rep9.sync.fet --min-count 5 --min-coverage 19 --max-coverage 383

perl /popoolation2/fisher-test.pl --input F0F60_Rep10.sync --output F0F60_Rep10.sync.fet --min-count 5 --min-coverage 14 --max-coverage 240

#The final SNP set, -log10-transformed p-values of CMH test and -log10-transformed p-values of FETs are provided in F0-F60SNP_CMH_FET_blockID.sync.zip file (Dryad Digital Repository: https://doi.org/10.5061/dryad.rr137kn) columns 1 to 84.

