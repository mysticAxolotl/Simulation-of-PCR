import matplotlib.pyplot as plt
from collections import defaultdict
import pandas
import random
import numpy
import time
import os

# param: a single strand 5' to 3' dna
# return: a single strand 5' to 3' dna that is reverse complement to the input strand
def reverse_compliment(dna5_3):  # input is a strand from 5' to 3'
    rc_dna5_3 = dna5_3.replace('A', 'X')  # replace A by X
    rc_dna5_3 = rc_dna5_3.replace('T', 'A')  # replace T by A
    rc_dna5_3 = rc_dna5_3.replace('X', 'T')  # replace X by A
    rc_dna5_3 = rc_dna5_3.replace('C', 'X')  # replace C by X
    rc_dna5_3 = rc_dna5_3.replace('G', 'C')
    rc_dna5_3 = rc_dna5_3.replace('X', 'G')
    rc_dna5_3 = rc_dna5_3[::-1]  # reverse the complementary strand to have a strand from 5' to 3'
    return rc_dna5_3

# param: a single strand 5' to 3' dna
# return: a single strand 5' to 3' dna that is reverse complement to the input strand
def compliment(dna5_3):  # input is a strand from 5' to 3'
    rc_dna5_3 = dna5_3.replace('A', 'X')  # replace A by X
    rc_dna5_3 = rc_dna5_3.replace('T', 'A')  # replace T by A
    rc_dna5_3 = rc_dna5_3.replace('X', 'T')  # replace X by A
    rc_dna5_3 = rc_dna5_3.replace('C', 'X')  # replace C by X
    rc_dna5_3 = rc_dna5_3.replace('G', 'C')
    rc_dna5_3 = rc_dna5_3.replace('X', 'G')
    return rc_dna5_3

# param: a list of tuples of 2 strs, representing double stranded dna segments
# return: a list of single strand dna segments
def denaturation(dna_segments):
    singleStrandDNAs = []
    singleStrandDNAs.append(dna_segments[0])
    singleStrandDNAs.append(dna_segments[1])
    return singleStrandDNAs

# param: a double strand dna, a tuple of 2 strings, representing 2 segments of dna from 5' to 3'
# return: a tuple of 2 strs representing the pair of primers (5' -> 3', GC content > 40%, bases btw the 2 primers: ~200)
def getPrimers():
    # primers: (sequence, start, end, GC content %), Primer 3 chosen
    f_primer = ('TGGACCCCAAAATCAGCGAA', 12, 31, 50)
    r_primer = ('GTGAGAGCGGTGAACCAAGA', 170, 151, 55)
    return f_primer, r_primer
# param: a list of tuples of 2 strings representing double stranded DNA segments
# return: null
def getStats(PCR_products, list_of_lengths):
    list_of_fragments = []
    forChart = defaultdict(int)
    for k in range(1, len(PCR_products)):
        if( PCR_products[k][0] != '' ):
            forChart[len( PCR_products[k][0])] += 1
            list_of_fragments.append(PCR_products[k][0])
        if( PCR_products[k][1] != '' ):
            forChart[len( PCR_products[k][1])] += 1
            list_of_fragments.append(PCR_products[k][1])
        
    # Print out the number of DNA fragments in PCR products
    # make an array, store lengths in it using a for loop
    print('The number of DNA fragments are: ', (len(list_of_fragments)))
    # Prints the maximum length of a DNA strand in PCR products
    print('The maximum length of the DNA fragments are: ', max(list_of_lengths))
    # Prints the minimum length of a DNA strand in PCR products
    print('The minimum length of the DNA fragments are: ', min(list_of_lengths))
    # Calculates the average length of all the DNA fragments in PCR products
    average_DNA = numpy.mean(list_of_lengths)
    print('The average length of the DNA fragments are: ', average_DNA)
    
    # Distribution of lengths of DNA fragments
    plt.xlabel("Sizes")
    plt.ylabel("Counts")
    plt.title("Size Distributions")
    plt.bar(forChart.keys(), forChart.values())
    plt.show()

    #Convert all of the ATGC to upper case then search for GC content over the total length of PCR products
    C_content = 0
    G_content = 0
    for i in list_of_fragments:
        C_content += i.count('C')
        G_content += i.count('G')
    base_average_gc = (C_content + G_content)
    average_gc = base_average_gc / len(PCR_products)
    print('The Average GC Content is: ', average_gc)
    return

# Run only once to generate n_gene.txt
if not(os.path.exists('./n_gene.txt')):
   file = open('SARS_COV2.fasta', 'r')
   comments = file.readline()
   SARS_COV2_genome = file.read()
   file.close()
   SARS_COV2_genome = SARS_COV2_genome.replace('\n','')

   #extract N gene from 28274:29533
   N_gene = SARS_COV2_genome[28273:29533]

   # Save n_gene to file to save time:
   file = open('n_gene.txt', 'w')
   file.write(N_gene)
   file.close()

# Load in the n_gene
file = open('n_gene.txt', 'r')
tDNA_N = file.read() # template DNA strand for N gene (lecture 7, 32:35)
file.close()

# tDNA_N has our entire DNA strand from 3' -> 5'
# If the above is true, then I believe we next generate the reverse compliment strand/coding strand:
cDNA_N = reverse_compliment(tDNA_N) # This is in 5' -> 3'

# Per her recomendation, keep the other one in 5' -> 3'
tDNA_N = tDNA_N[::-1]

# Store these as a tuple, 5' -> 3'
DNA_N = (cDNA_N, tDNA_N)

f_primer, r_primer = getPrimers()
f_compliment = compliment(f_primer[0])
r_compliment = compliment(r_primer[0])

initial_strands = [DNA_N]
random.seed(time.time())

start_time = time.time()
for i in range(4):
    copies = []
    for strand in initial_strands:
        temp = []
        single_segments = denaturation(strand)
        for j in range(0, len(single_segments)):
            falloff = random.randint(-50, 50) + 178
            if r_compliment in single_segments[j]:
                start_index = single_segments[j].find(r_compliment) + len(r_primer[0])
                end_index = start_index + falloff
                temp1 = compliment(single_segments[j][start_index:end_index])
                temp1 = r_primer[0] + temp1
                temp.append(temp1)
            elif f_compliment in single_segments[j][::-1]:
                start_index = single_segments[j][::-1].find(f_compliment)
                end_index = start_index + falloff
                temp1 = compliment(single_segments[j][::-1][start_index:end_index])[::-1]
                temp1 = f_primer[0] + temp1
                temp.append(temp1)
            else:
                temp.append('')
        copies.append(tuple(temp))
    initial_strands.extend(copies)

print(str((time.time() - start_time)/60) + " Minutes")

segment_lengths = []
for pair in initial_strands:
       for strand in pair:
           if strand != '':
               segment_lengths.append(len(strand))
print(len(segment_lengths))
getStats(initial_strands, segment_lengths)
