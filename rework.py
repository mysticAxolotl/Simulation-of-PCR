import matplotlib.pyplot as plt
from collections import defaultdict
import random
import time
import os

random.seed(time.time())

def compliment( dna ):
    compDict = { 'A':'T', 'T':'A', 'G':'C', 'C':'G' }
    comp = "";
    for i in dna:
        comp += compDict[ i ]
    return comp

def getTotalGC( dna ):
    totalGC = 0
    for i in dna:
        toTest = DNA_N[ i[0]: i[1] ]
        totalGC += toTest.count('G') + toTest.count('C')
    return totalGC

def stats( fragments, lengths ):

    minL, maxL, totalGC, totalL, totalLF, totalF = 1259, 0, getTotalGC( fragments ), 0, 0,  len( fragments )
    for key, value in lengths.items():
        if key < minL:
            minL = key
        elif key > maxL:
            maxL = key
        totalL += key * value
        totalLF += value

    print( totalLF )    
    print( "The number of DNA fragments are: ", totalF * 2 )
    print( "The maximun length of the DNA fragments are: ", maxL )
    print( "The minimum length of the DNA fragments are: ", minL )
    print( "The average length of the DNA fragments are: ", totalL / ( totalF * 2 ) )
    
    # Distribution of lengths of DNA fragments
    plt.xlabel("Sizes")
    plt.ylabel("Counts")
    plt.title("Size Distributions")
    plt.bar( lengths.keys(), lengths.values())
    plt.show()

    #Convert all of the ATGC to upper case then search for GC content over the total length of PCR products
    print( "The Average GC Content is: ", totalGC / totalF)
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
f_primer, r_primer = compliment( "TGGACCCCAAAATCAGCGAA" ), compliment ( "GTGAGAGCGGTGAACCAAGA" )[::-1]
# primers = [ ('TGGACCCCAAAATCAGCGAA', 12, 31, 50), ('GTGAGAGCGGTGAACCAAGA', 170, 151, 55) ]

# all dna we've copied and collected 
tDNA = []
lengthDict = defaultdict(int)
lastCycle = [ ( 0, 1259, 0, 1259 ) ]
for i in range( 10 ):
    start, end, rstart, rend = 0, 0, 0, 0
    temp = []
    for strand in lastCycle:
        falloff = random.randint(-50, 50) + 178
        rfalloff = random.randint( -50, 50 ) + 178
        start = rDNA_N[ strand[2]: strand[3] ].find( f_primer )
        rend = DNA_N[ strand[0]:strand[1] ].find( r_primer )
        if start != -1:
            end = start + falloff + 20
            if start < strand[2]:
                start = strand[2]
            temp.append( ( start, end, strand[2], strand[3] ) )
            lengthDict[ end - start ] += 1
            lengthDict[ strand[3] - strand[2] ] += 1

        if rend != -1:
            rstart = rend - rfalloff
            end += 20
            if rstart < strand[0]:
                rstart = strand[0]
            temp.append( ( strand[0], strand[1], rstart, rend ) )
            lengthDict[ rend - rstart ] += 1
            lengthDict[ strand[1] - strand[0] ] += 1

    tDNA.extend( lastCycle )
    lastCycle.clear()
    lastCycle.extend( temp )
lengthDict.pop( 1259 )
stats( tDNA, lengthDict )