
import matplotlib
import pandas
import numpy
import os

#
#     Given Functions/Code
#

# param: a single strand 5" to 3" dna
# return: a single strand 5" to 3" dna that is reverse complement to the input strand
def reverse_complement(dna5_3):  # input is a strand from 5" to 3"
    rc_dna5_3 = dna5_3.replace('A', 'X')  # replace A by X
    rc_dna5_3 = rc_dna5_3.replace('T', 'A')  # replace T by A
    rc_dna5_3 = rc_dna5_3.replace('X', 'T')  # replace X by A
    rc_dna5_3 = rc_dna5_3.replace('C', 'X')  # replace C by X
    rc_dna5_3 = rc_dna5_3.replace('G', 'C')
    rc_dna5_3 = rc_dna5_3.replace('X', 'G')
    rc_dna5_3 = rc_dna5_3[::-1]  # reverse the complementary strand to have a strand from 5" to 3"
    return rc_dna5_3

# param: a list of tuples of 2 strs, representing double stranded dna segments
# return: a list of single strand dna segments
def denaturation(dna_segments):
    singleStrandDNAs = []
    return singleStrandDNAs

### POTENETIALL DO NOT NEED
# param: a double strand dna, a tuple of 2 strings, representing 2 segments of dna from 5" to 3"
# return: a tuple of 2 strs representing the pair of primers (5" -> 3", GC content > 40%, bases btw the 2 primers: ~200)
def getPrimers():
    # primers: (sequence, start, end, GC content %), Primer 3 chosen
    f_primer = ('TGGACCCCAAAATCAGCGAA', 12, 31, 50)
    r_primer = ('GTGAGAGCGGTGAACCAAGA', 170, 151, 55)
    return f_primer, r_primer


# param: a list of single strand dna segments, each segment is from 5" to 3"
# return: a list of tuples of 2 strs (2 dna segments from 5" to 3")
def annealing_elongation(singleStrandDNAs, primers, fall_of_rate):
    # ...
    doubleStrandedDNAs = [('a','a'),('a','a')]  # get your sequence of dnas
    return doubleStrandedDNAs


# param: gene to be copied (a tuple of 2 strs, 5' -> 3'), fall of rate of DNA polymerase (int), and num_cycles to run PCR (int)
# return: a list of double stranded dna segments
def PCR(dna_segment_to_be_copied, fall_of_rate, num_cycles):
    # ....

    PCRproducts = [dna_segment_to_be_copied]
    for i in range(num_cycles):
        singleStrandDNAs = denaturation(PCRproducts)
        PCRproducts = annealing_elongation(singleStrandDNAs, getPrimers(), fall_of_rate)

    return PCRproducts


# param: a list of tuples of 2 strings representing double stranded DNA segments
# return: null
def getStats(PCR_products):
    # ...
    dna_segment_lengths = []
    dna_gc_contents = []
    dna_segment_lengths.append(len(pair[0]))
    dna_segment_lengths.append(len(pair[1]))
    # Print out the number of DNA fragments in PCR products
    print('The number of DNA fragements are: ', (2 * len(PCR_products))
    # Prints the maximum lenth of a DNA strand in PCR products
    print('The maximum length of the DNA fragments are: ', max(PCR_products))
    # Prints the minimum length of a DNA strand in PCR products
    print('The minimum length of the DNA fragments are: ', min(PCR_products))
    # Calculates the average length of all the DNA fragments in PCR products
    temp = (PCR_products)
    average_DNA = (float(sum(temp) / len(temp)))
    print('The average length of the DNA fragments are: ', average_DNA)
    # Distribution of lengths of DNA fragments, uses temp as it is the length of PCR products
    hist = plt.hist(PCR_products) #elements a
    plt.figure()
    plt.xlabel('Length range')
    plt.ylabel('Frequency')
    plt.title('Distribution of Lengths of DNA strands')
    #Convert all of the ATGC to upper case then search for GC content over the total length of PCR products
    average_gc_temp = PCR_products.upper()
    dna_c_content = pair[1].count('C')
    dna_g_content = pair[1].count('G')
    average_gc = (dna_c_content + dna_g_content)
    print('The Average GC Content is: ', average_gc)
    return


""" ccDNA_N = N_gene #complementary strand of the cDNA strand

   DNA_N = (cDNA_N, ccDNA_N)

   primers = getPrimers(DNA_N)

   PCR_products = PCR(DNA_N, primers, fall_off_rate)

   # the gene needs to be amplified
   S_gene = "ATG..."
   cDNA_S = reverse_complement(S_gene)
   rccDNA_S = S_gene

   # It is the double strand DNA {cDNA_S, rccDNA_S} that can be amplified.
   DNA_S = (cDNA_S, rccDNA_S)

   fall_of_rate = 50  # between -50 to 50
   num_cycles = 20  # PCR cycles

   # call PCR function to simulation the PCR process
   PCR_products = PCR(DNA_S, fall_of_rate, num_cycles)

   # print stats of PCR_products
   getStats(PCR_products)
"""

# Potential New Code, TODO: Fit into the above functions

import matplotlib
import pandas
import numpy
import os

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

# TODO: Correct me if I am wrong, tDNA_N should have our entire DNA strand from 3' -> 5'
# If the above is true, then I believe we next generate the reverse compliment strand/coding strand:
cDNA_N = reverse_complement(tDNA_N) # This is in 5' -> 3'

# Per her recomendation, keep the other one in 5' -> 3'
tDNA_N = tDNA_N[::-1]

# Store these as a tuple, 5' -> 3'
DNA_N = (cDNA_N, tDNA_N)

# TODO: For the annealing_elongation portion, translating what she said, I think the forward primer
# goes on the start of tDNA_N after tDNA_N has been reversed. The reverse primer goes on the
# end of cDNA_N as it is now.

f_primer, r_primer = getPrimers()
primerStrands = ((DNA_N[1][::-1])[f_primer[1] - 1:r_primer[1]], (DNA_N[0][::-1][f_primer[1] - 1:r_primer[1]])[::-1])
print(primerStrands[0] + '\n')
print(primerStrands[1] + '\n')
