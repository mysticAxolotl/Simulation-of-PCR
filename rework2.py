import matplotlib.pyplot as plt
from collections import defaultdict
import random
import time
import os
import sys

random.seed(time.time())

def averageGC_totalLength_max_min( dna, segments, lengths ):
    GCDict = { 'A': 0, 'T': 0, 'G': 0, 'C':0 }
    GCtotal, Ltotal, mx, mn = 0, 0, 0, 1259
    for i in segments:
        temp = dna[ i[0]:i[1] ] + dna[ i[2]:i[3] ]
        GCtotal += temp.count('G') + temp.count('C')
    
    for key, value in lengths.items():
        Ltotal += key * value
        if key < mn:
            mn = key
        if key > mx:
            mx = key
    
    return GCtotal / Ltotal * 100, Ltotal, mn, mx

def stats( dna, segments, lengths, cycles ):
    Ftotal, Stotal = len( segments ) * 2, len( segments )
    averageGC, Ltotal, Lmin, Lmax = averageGC_totalLength_max_min( dna, segments, lengths )

    print( "The number of DNA segments are:", Stotal )
    print( "The number of DNA fragments are:", Ftotal )
    print( "The maximun length of the DNA fragments are:", Lmax )
    print( "The minimum length of the DNA fragments are:", Lmin )
    print( "The average length of the DNA fragments are:", format( Ltotal / Ftotal, ".4f" ) )
    print( "Segment loss is:", 2 ** cycles - Stotal )
    print( "Segment loss rate is", format( ( 2 ** cycles - Stotal ) / 2 ** cycles * 100, ".4f" ), "percent" )
    print( "The average GC content is:", format( averageGC, ".4f" ), "percent" )
    
    # Distribution of lengths of DNA fragments
    # plt.xlabel("Sizes")
    # plt.ylabel("Counts")
    # plt.title("Size Distributions")
    # plt.bar( lengths.keys(), lengths.values())
    # plt.show()
    return

def compliment( dna ):
    compDict = { 'A':'T', 'T':'A', 'G':'C', 'C':'G' }
    comp = ""
    for i in dna:
        comp += compDict[ i ]
    return comp

def PCR( cycles ):
    totalDNA, lengths, lastCycle = [ ( 0, 1259, 0, 1259 ) ], defaultdict(int), [ ( 0, 1259, 0, 1259 ) ]
    fstart, fend, rstart, rend = 0, 0, 0, 0
    for cycles in range( cycles ):
        temp = []
        for i in lastCycle:
            fstart, rend = 11, 170
            if not fstart < i[2]:
                fend = fstart + 178 + random.randint( -25, 25 )
                if fend > i[3]:
                    fend = i[3]
                temp.append( ( fstart, fend, i[2], i[3] ) )
                lengths[ fend - fstart ] += 1
            
            if not rend > i[1]:
                rstart = rend - 20 - 178 + random.randint( -25, 25 )
                if rstart < i[0]:
                    rstart = i[0]
                temp.append( ( i[0], i[1], rstart, rend ) )
                lengths[ rend - rstart ] += 1
        totalDNA.extend( lastCycle )
        lastCycle.clear()
        lastCycle.extend( temp )
    return totalDNA, lengths

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
file = open( 'n_gene.txt', 'r' )
dna = file.read() # template DNA strand for N gene (lecture 7, 32:35)
file.close()
# Primer indexes within the dna string
# We do not need to store the primer strings
# Forward: 11:31
# Reverse: 50:70

cycles = int( sys.argv[1] )
segments, lengths = PCR(  cycles )
stats( dna, segments, lengths, cycles )