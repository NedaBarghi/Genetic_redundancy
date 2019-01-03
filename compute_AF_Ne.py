import argparse
from argparse import RawTextHelpFormatter
from operator import itemgetter
import heapq
import operator
import numpy as np

# Author: Neda Barghi

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""description""")
parser.add_argument('--input', required=True, dest='Input', type=str, help="Input file contains polymorphic sites in sync format.")
parser.add_argument('--output', required=True, dest='Output', type=str, help="Output file contains the allele frequency and coverage of polymorphic sites. Allele frequencies are polarized based on the rising alleles. Each SNP will have computed allele frequencies followed by the coverage as specified in the header on the file.")

args = parser.parse_args()

#-------------------
# function
#-------------------

def get_sum_of_positions(cols_of_reps,position):
    return sum([int(rep.split(':')[position]) for rep in cols_of_reps])

#-------------------
# main
#-------------------

InputFile=open(args.Input,"r")
OutputFile=open(args.Output,"w")


generations = [0,60]
replicates = 10
first_column_alleleFreq = 4

AF_label = ['F'+str(g)+'R'+str(r) for g in generations for r in range(1,11)]
cov_label = ['F'+str(g)+'R'+str(r)+'_cov' for g in generations for r in range(1,11)]

OutputFile.write('\t'.join(['chr','pos','base']+AF_label+cov_label)+'\n')

index_sync = {'A':0, 'T':1, 'C':2, 'G':3, '*':5 }
for line in InputFile:
    line=line.rstrip()
    cols=line.split('\t')
    sum_allele_freq = [get_sum_of_positions(cols[first_column_alleleFreq-1:first_column_alleleFreq-1+replicates],i) for i in [0,1,2,3,5]] 
    major_allele_index,minor_allele_index =(sum_allele_freq.index(heapq.nlargest(2, sum_allele_freq)[0]),sum_allele_freq.index(heapq.nlargest(2, sum_allele_freq)[1]))
    if major_allele_index == minor_allele_index: 
        sorted_sum_allele_freq = sorted(enumerate(sum_allele_freq), reverse=True, key=itemgetter(1)) 
        major_allele_index,minor_allele_index = (sorted_sum_allele_freq[0][0],sorted_sum_allele_freq[1][0])
    if major_allele_index == 4: major_allele_index=5
    if minor_allele_index == 4: minor_allele_index=5
    allele_count = [int(rep.split(':')[major_allele_index]) + int(rep.split(':')[minor_allele_index]) for rep in cols[first_column_alleleFreq-1:first_column_alleleFreq-1+(len(generations)*replicates)]]
    allele_freq = []
    for rep in cols[first_column_alleleFreq-1:first_column_alleleFreq-1+(len(generations)*replicates)]:
        if float(rep.split(':')[major_allele_index]) > 0.0 or float(rep.split(':')[minor_allele_index]) > 0.0:
            allele_freq.append(float(rep.split(':')[minor_allele_index])/ (float(rep.split(':')[major_allele_index]) + float(rep.split(':')[minor_allele_index])))
        if float(rep.split(':')[major_allele_index]) == 0.0 and float(rep.split(':')[minor_allele_index]) == 0.0:
            allele_freq.append('NA')             
    modified_AF=[0 if i == 'NA' else i for i in allele_freq]    
    AFC_alltimepoints = [map(operator.sub, modified_AF[replicates*i:replicates*(i+1)], modified_AF[0: replicates]) for i in range(1,len(generations))]    
    ave_all_reps=map(np.mean,[[rep[i] for rep in AFC_alltimepoints] for i in range(0,replicates)])    
    if np.mean(ave_all_reps) >= 0: #minor rising
        OutputFile.write('\t'.join(cols[0:first_column_alleleFreq-1])+'\t'+'\t'.join([str(i) for i in allele_freq])+'\t'+'\t'.join([str(i) for i in allele_count])+'\n')
    if np.mean(ave_all_reps) < 0: #major rising
        allele_freq_pol = ['NA' if item == 'NA' else 1-item for item in allele_freq]
        OutputFile.write('\t'.join(cols[0:first_column_alleleFreq-1])+'\t'+'\t'.join([str(i) for i in allele_freq_pol])+'\t'+'\t'.join([str(i) for i in allele_count])+'\n')
OutputFile.flush()
InputFile.close()
