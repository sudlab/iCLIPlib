''' This a modeule that holds functions and classes useful for analysing iCLIP data '''

import CGAT.Experiment as E
import numpy as np
import pandas as pd

def find_first_deletion(cigar):

    position = 0
    for operation, length in cigar:

        if operation == 2:
            return position
        else:
            position += length

    return position



def countChr(reads, chr_len, dtype = 'uint16'):
    ''' Counts the crosslinked bases for each read in the pysam rowiterator
    reads and saves them in pandas Series: those on the positive strand
    and those on the negative strand. The Series are indexed on genome position,
    and are sparse.

    Cross linked bases are definated as in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the iCLIP cDNAs mapped by Bowtie was used to
        define the cross link sites identified by truncated cDNAs.

        [For reads with deletions] The deleted nucleotide in CLIP and iCLIP
        cDNAs mapped by Novoalign was used to define the cross-link sites
        identified by read-through cDNAs. If a cDNA had more than one deletion,
        we selected the one closest to the beginning of the read.
    
    The dtype to use internally for storage can be specified. Large types
    reduce the chance of overflow, but require more memory. With 'uint16'
    the largest count that can be handled is 255. Data is stored sparse,
    so memory is less of a problem. Overflow will cause a ValueError.

    returns a tuple of pandas Series objects, with the positive and negative
    strand arrays and also a counter object that contains the counts for each
    type of site. '''

    pos_depths = pd.Series(dtype=dtype)
    neg_depths = pd.Series(dtype=dtype)

    counter = E.Counter()

    for read in reads:
        
        if not 'D' in read.cigarstring:
            if read.is_reverse:
                series = neg_depths
                pos = read.aend
                counter["truncated_neg"] += 1
            else:
                series = pos_depths
                pos = read.pos - 1
                counter["truncated_pos"] += 1

        else:
            if read.is_reverse:
                cigar = reversed(read.cigar)
                position = find_first_deletion(cigar)
                pos = read.aend - position - 1
                series = neg_depths
                counter["deletion_neg"] += 1
            else:
                position = find_first_deletion(read.cigar)
                pos = read.pos + position
                series = pos_depths
                counter["deletion_pos"] += 1
        
        try:
            series.loc[pos] += 1
        except KeyError:
            series.loc[pos] = 1
        
    #check for integer overflow: counter sum should add up to array sum
    counter_sum = sum([y for x,y in counter.iteritems()])
    array_sum = pos_depths.sum() + neg_depths.sum()
    if not counter_sum == array_sum:
        raise (ValueError,
               "Sum of depths is not equal to number of "
               "reads counted, possibly dtype %s not large enough" % dtype)
    
    E.debug("Counted %i truncated on positive strand, %i on negative"
            % (counter.truncated_pos, counter.truncated_neg))
    E.debug("and %i deletion reads on positive strand, %i on negative"
            % (counter.deletion_pos, counter.deletion_neg))
    
    return (pos_depths, neg_depths, counter)
