import matplotlib.pyplot as plt
from collections import defaultdict
import random
import time
import os
import sys

random.seed(time.time())

#-------------------------------------------------------------
# Gets the N gene from SARS_COV_2 and returns the string
#
def getNGene():
    # Run only once to generate n_gene.txt
    if not os.path.exists( "./n_gene.txt" ):
        file = open( "SARS_COV2.fasta", 'r' )
        file.readline()
        SARS_COV2_genome = file.read()
        file.close()
        SARS_COV2_genome = SARS_COV2_genome.replace( '\n', '' )

        N_gene = SARS_COV2_genome[ 28273: 29533 ]

        file = open( "n_gene.txt", 'w' )
        file.write( N_gene )
        file.close()

    file = open( "n_gene.txt", 'r' )
    dna = file.read()
    file.close()

    return dna

#-------------------------------------------------------------
# Finds average GC, total length of fragments, and the min/max
# lengths of all fragments
#   $dna: string
#       The strand of dna we copied in PCR 
#   $segments: array of tuples 
#       All the dna segments we copied
#
def averageGC_totalLength_max_min( dna, segments ):

    # $GCtotal: total # of 'G' and 'C' in all fragments
    # $Ltotal: total length of all fragments
    # $mx, $mn: max/min length of all fragments
    GCtotal, Ltotal, Stotal, mx, mn = 0, 0, 0, 0, 1259

    lengths = defaultdict(int)
    
    # Loops through our dna segments and counts the GC content of each
    # fragment, and populates the lengths dict
    for i in segments:
        lengths[ i[1] - i[0] ] += 1
        lengths[ i[3] - i[2] ] += 1
        temp = dna[ i[0]:i[1] ] + dna[ i[2]:i[3] ]
        GCtotal += temp.count('G') + temp.count('C')
    lengths.pop( 1259 )

    # Loops through the length dictionary and adds all lengths
    # finds min/max
    for key, value in lengths.items():
        Stotal += value
        Ltotal += key * value
        if key < mn:
            mn = key
        if key > mx:
            mx = key

    # Average gc is in a percentage
    return GCtotal / Ltotal * 100, Ltotal / Stotal, mn, mx, lengths

#-------------------------------------------------------------
# Prints the stats of the segments we've collected in pcr, and
# the histogram of lengths
#   $dna: string
#       The strand of dna we copied in PCR 
#   $segments: array of tuples 
#       All the dna segments we copied
#   $cycles: int
#       number of pcr cycles that were run
#
def stats( dna, segments, cycles ):
    # Ftotal: Number of fragments
    # Stotal: Number of segments
    # averageGC, averageL: Average GC content/length
    # Lmin, Lmax: Max/min of lenghts
    # Lengths: Dict of lengths and how often they occur
    #          needed for the histogram
    # precision: How many decimal places we want to display when 
    #            printing stats
    Ftotal, Stotal = len( segments ) * 2, len( segments )
    averageGC, averageL, Lmin, Lmax, lengths = averageGC_totalLength_max_min( dna, segments )
    precision = ".3f"

    print(    "The number of DNA segments are:", Stotal, \
            "\nThe number of DNA fragments are:", Ftotal, \
            "\nThe maximun length of the DNA fragments are:", Lmax, \
            "\nThe minimum length of the DNA fragments are:", Lmin, \
            "\nThe average length of the DNA fragments are:", format( averageL, precision ), \
            "\nSegment loss is:", 2 ** cycles - Stotal, \
            "\nSegment loss rate is", format( ( 2 ** cycles - Stotal ) / 2 ** cycles * 100, precision ), "percent", \
            "\nThe average GC content is:", format( averageGC, precision ), "percent" )
    
    # Distribution of lengths of DNA fragments
    plt.xlabel("Size")
    plt.ylabel("Count")
    plt.title("Size Distributions")
    plt.bar( lengths.keys(), lengths.values())
    plt.show()
    return

#-------------------------------------------------------------
# Runs the entire PCR process
#   $cycles: int
#       number of PCR cycles to run
#
def PCR( cycles, falloff = ( 178, 25, 35 ) ):
    # totalDNA: all the dna segments we copy
    # lastCycle: The copies from the last cycle we ran
    # fstart/fend: start/end points where the forward primer binds to the dna
    # rstart/rend: start/end points where the reverse primer binds to the dna
    totalDNA, lastCycle  = [ ( 0, 1259, 0, 1259 ) ], [ ( 0, 1259, 0, 1259 ) ]
    fstart, fend, rstart, rend = 0, 0, 0, 0
    print( "Fall off is: ", falloff[0], " + rand( -", falloff[1], ", ", falloff[2],")", sep='' )

    for cycles in range( cycles ):
        # stores all copies we've made this cycle
        temp = []
        for i in lastCycle:
            # defines where the first primer begins and the reverse ends 
            # in the dna string
            fstart, rend = 11, 170

            # if fstart < i[2] means there isn't enough dna to bind to properly
            if not fstart < i[2]:
                fend = fstart + falloff[0] + random.randint( -falloff[1], falloff[2] )
                # fend cannot go past the dna we're copying since it doesn't exist
                if fend > i[3]:
                    fend = i[3]
                temp.append( ( fstart, fend, i[2], i[3] ) )
            
            # if rend is > i[1] again means there isn't enough dna to bind to properly
            if not rend > i[1]:
                rstart = rend - 20 - falloff[0] + random.randint( -falloff[1], falloff[2] )
                if rstart < i[0]:
                    rstart = i[0]
                temp.append( ( i[0], i[1], rstart, rend ) )
        
        # adds last cycle to total DNA, clears last cycle, then adds the copies we 
        # got in the current cycle
        totalDNA.extend( lastCycle )
        lastCycle.clear()
        lastCycle.extend( temp )
    return totalDNA

# Primer indexes within the dna string: 
# We do not need to store the primer strings
# or their compliments
#   Forward: 11:31      TGGACCCCAAAATCAGCGAA
#   Reverse: 50:70      GTGAGAGCGGTGAACCAAGA

cycles = int( sys.argv[1] )
if len( sys.argv ) > 2:
    stats( getNGene(), PCR( cycles, ( int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) ) ), cycles ) 
else:
    stats( getNGene(), PCR( cycles ), cycles )