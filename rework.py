import matplotlib.pyplot as plt
from collections import defaultdict
import pandas
import random
import numpy
import time
import os

random.seed(time.time())

def compliment( dna ):
    compDict = { 'A':'T', 'T':'A', 'G':'C', 'C':'G' }
    comp = "";
    for i in dna:
        comp += compDict[ i ]
    return comp


def stats( fragments, lengths ):
    minL = 1259
    maxL = 0
    totalL = 0
    totalGC = 0
    for i in fragments:
        toTest = DNA_N[ i[0]: i[1] ]
        rtoTest = rDNA_N[ i[0]: i[1] ]
        length = len( toTest )
        rlength = len( rtoTest )
        totalL += length + rlength
        if length < minL:
            minL = length
        elif length > maxL:
            maxL = length
        
        if rlength < minL:
            minL = rlength
        elif rlength > maxL:
            maxL = rlength
        totalGC += toTest.count('C') + toTest.count('G') + rtoTest.count('C') + rtoTest.count('G') 
    
    print( "The number of DNA fragments are: ", len( fragments ) )
    print( "The maximun length of the DNA fragments are: ", maxL )
    print( "The minimum length of the DNA fragments are: ", minL )
    print( "The average length of the DNA fragments are: ", totalL / ( len( fragments ) * 2 ) )
    
    # Distribution of lengths of DNA fragments
    plt.xlabel("Sizes")
    plt.ylabel("Counts")
    plt.title("Size Distributions")
    plt.bar( lengths.keys(), lengths.values())
    plt.show()

    #Convert all of the ATGC to upper case then search for GC content over the total length of PCR products
    averageGC = totalGC / len( fragments ) * 2
    print( "The Average GC Content is: ", averageGC)
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
DNA_N = file.read() # template DNA strand for N gene (lecture 7, 32:35)
file.close()

rDNA_N = compliment( DNA_N )
# primers: (sequence, start, end, GC content %), Primer 3 chosen
primers = [ ('TGGACCCCAAAATCAGCGAA', 12, 31, 50), ('GTGAGAGCGGTGAACCAAGA', 170, 151, 55) ]
f_primer = primers[0][0]
r_primer = primers[1][0][::-1]
fc_primer = compliment(f_primer)
rc_primer = compliment(r_primer)

# all dna we've copied and collected 
tDNA = [ ( 0, 1259 ) ]
lengthDict = defaultdict(int)
for i in range( 2 ):
    start = 0
    end = 0
    temp = []
    for strand in tDNA:
        falloff = random.randint(-50, 50) + 178
        toTest = rDNA_N[ strand[0]: strand[1] ] 
        start = toTest.find(fc_primer) 
        if start != -1:
            end = start + falloff + 20
            lengthDict[ end - start ] += 1
            temp.append( ( start, end ) )

        falloff = random.randint(-50, 50) + 178
        toTest = DNA_N[ strand[0]:strand[1] ]
        end = toTest.rfind( rc_primer )
        if start != -1:
            end += 20
            start = end - falloff
            if start < 0:
                start = 0
            lengthDict[ end - start ] += 1
            temp.append( ( start, end ) )
    for i in temp:
        tDNA.append( i )
tDNA.remove( ( 0, 1259 ) )
stats( tDNA, lengthDict )