import os
file = open('SARS_COV2.fasta', 'r')
comments = file.readline()
SARS_COV2_genome = file.read()
file.close()
SARS_COV2_genome = SARS_COV2_genome.replace('\n','')

#extract N gene from 28274:29533
N_gene = SARS_COV2_genome[28273:29533]

cDNA_N = N_gene.replace('A', 'X') # replace A with X
cDNA_N = cDNA_N.replace('T', 'A')
cDNA_N = cDNA_N.replace('T', 'A')
cDNA_N = cDNA_N.replace('T', 'A')
cDNA_N = cDNA_N.replace('T', 'A')
cDNA_N = cDNA_N.replace('T', 'A')

ccDNA_N = N_gene #complementary strand of the cDNA strand

DNA_N = (cDNA_N, ccDNA_N)

primers = getPrimers(DNA_N)

PCR_products = PCR(DNA_N, primers, fall_off_rate)