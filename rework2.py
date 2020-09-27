import matplotlib.pyplot as plt
from collections import defaultdict
import random
import time
import os

random.seed(time.time())

def averageGC_totalLength( dna, segments, lengths ):
    GCtotal = 0
    Ltotal = 0
    for i in segments:
        temp = dna[ i[0]:i[1] ]
        GCtotal += temp.count( 'G' ) + temp.count( 'C' )
        # GCtotal += temp.count( 'G' or 'C' )
    for key, value in lengths.items():
        Ltotal += key * value
    return ( GCtotal * 2 / Ltotal ) * 100, Ltotal

def stats( dna, segments, lengths ):
    Lmax, Lmin, Ftotal = 0, 1259, len( segments ) * 2
    averageGC, Ltotal = averageGC_totalLength( dna, segments, lengths )
    
    for key, value in lengths.items():
        if key < Lmin:
            Lmin = key
        elif key > Lmax:
            Lmax = key
    print( "The number of DNA segments are:", Ftotal / 2 )
    print( "The number of DNA fragments are:", Ftotal )
    print( "The maximun length of the DNA fragments are:", Lmax )
    print( "The minimum length of the DNA fragments are:", Lmin )
    print( "The average length of the DNA fragments are:", Ltotal / Ftotal )
    print( "Segment loss rate is:", ( 2 ** 10 - ( Ftotal / 2 ) ) / 2 ** 10 * 100 )
    print( "The average GC content is:", averageGC, "percent" )
    
    # Distribution of lengths of DNA fragments
    plt.xlabel("Sizes")
    plt.ylabel("Counts")
    plt.title("Size Distributions")
    plt.bar( lengths.keys(), lengths.values())
    plt.show()
    return

def compliment( dna ):
    compDict = { 'A':'T', 'T':'A', 'G':'C', 'C':'G' }
    comp = ""
    for i in dna:
        comp += compDict[ i ]
    return comp

def PCR( dna, dnaCompliment, fPrimer, rPrimer ):
    totalDNA, lengths, lastCycle = [], defaultdict(int), [ ( 0, 1259, 0, 1259 ) ]
    fstart, fend, rstart, rend = 0, 0, 0, 0
    for cycles in range( 11 ):
        temp = []
        for i in lastCycle:
            fstart = dnaCompliment[ i[0]: i[1] ].find( fPrimer )
            rend = dna[ i[2]: i[3] ].find( rPrimer )
            if fstart != -1:
                fend = fstart + 188 + random.randint(-50, 50)
                if fend > i[3]:
                    fend = i[3]
                temp.append( ( fstart, fend, i[2], i[3] ) )
                lengths[ fend - fstart ] += 1
            
            if rend != -1:
                rstart = rend - 188 + random.randint(-50, 50)
                rend += 20
                if rstart < i[0]:
                    rstart = i[0]
                temp.append( ( i[0], i[1], rstart, rend) )
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
dnaCompliment = compliment( dna )
f_primer, r_primer = compliment( "TGGACCCCAAAATCAGCGAA" ), compliment ( "GTGAGAGCGGTGAACCAAGA" )[::-1] 

segments, lengths = PCR( dna, compliment( dna ), compliment( "TGGACCCCAAAATCAGCGAA" ), compliment( "GTGAGAGCGGTGAACCAAGA" )[::-1] )
stats( dna, segments, lengths )