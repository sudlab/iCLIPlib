''' This a modeule that holds functions and classes useful for analysing iCLIP data '''

import CGAT.Experiment as E
import numpy as np
import pandas as pd
import CGAT.GTF as GTF
import collections

def find_first_deletion(cigar):
    '''Find the position of the the first deletion in a 
    read from the cigar string, will return 0 if no deletion 
    found '''

    position = 0
    for operation, length in cigar:

        if operation == 2:
            return position
        else:
            position += length

    return position


class TranscriptCoordInterconverter:
    ''' A class to interconvert between genome co-ordinates
    and transcript co-ordinates. Implemented as a class because
    there are expected to be many calls against the same transcript,
    so time can be saved by precomputation 

    TranscriptCoordInterconverter.genome2transcript should be the 
    interverse of TranscriptCoordInterconverter.transcript2genome.

    That is 

    if myConverter = TranscriptCoordInterverter(transcript)
    
    then
    
    myConverter.genome2transcript(myConverter.transcript2genome(x)) == x
    
    and

    myConverter.transcript2genome(myConverter.genome2transcript(x)) == x'''

    

    def __init__(self, transcript, introns=False):
        ''' Pre compute the conversions for each exon '''

        if not introns:
            intervals = GTF.asRanges(transcript, feature="exon")
        else:
            intervals = GTF.toIntronIntervals(transcript)
        
        # get strand
        self.strand = transcript[0].strand

        # store transcript_id
        self.transcript_id = transcript[0].transcript_id

        # sort the exons into "transcript" order
        if self.strand == "-":
            intervals.sort(reverse=True)
            intervals = [(y-1, x-1) for x, y in intervals]
        else:
            intervals.sort(reverse=False)

        self.offset = intervals[0][0]
        self.genome_intervals = [map(abs, (x-self.offset, y-self.offset))
                                 for x, y in intervals]

        interval_sizes = [abs(y-x) for x, y in intervals]

        total = 0
        transcript_intervals = [None]*len(interval_sizes)

        for i in range(len(interval_sizes)):
            transcript_intervals[i] = (total,
                                       interval_sizes[i] + total)
            total += interval_sizes[i]
        
        self.transcript_intervals = transcript_intervals
        self.length = transcript_intervals[-1][1]

    def genome2transcript(self, pos):
        ''' Convert genome coordinate into transcript coordinates.
        pos can be a single value or a nunpy array like object.
        Passing an array ensures that the transcript is only
        searched once, ensuring O(n) performance rather than
        O(nlogn)'''

        if len(pos) == 0:
            return np.array([])

        try:
            relative_pos = pos - self.offset
        except TypeError:
            relative_pos = np.array(pos) - self.offset
        
        if self.strand == "-":
            relative_pos = relative_pos * -1

        ordering = np.argsort(relative_pos)
        relative_pos = np.sort(relative_pos)

        # pre allocate results list for speed
        try:
            results = np.zeros(len(relative_pos))
        except TypeError:
            relative_pos = np.array([relative_pos])
            results = np.zeros(1)

        i = 0
        i_max = len(relative_pos)

        # do linear search for correct exon
        for exon, interval in enumerate(self.genome_intervals):

            if relative_pos[i] < interval[0]:
              
                
                raise ValueError("Position %i is not in transcript %s" %
                                 (pos[i], self.transcript_id) )
            
            while relative_pos[i] < interval[1]:
                
                pos_within_exon = relative_pos[i]-interval[0]
                transcript_exon = self.transcript_intervals[exon]
                transcript_position = transcript_exon[0] + pos_within_exon
                
                results[i] = transcript_position
                i += 1
                if i == i_max:
                    return results[ordering]

        # exon has not been found
       
        raise ValueError("Position %i is not in transcript %s" %
                         (pos, self.transcript_id))

    def transcript2genome(self, pos):
        ''' Convert transcript coodinate into genome coordinate,
        pos can be a single value or a nunpy array like object.
        Passing an array ensures that the transcript is only
        searched once, ensuring O(n) performance rather than
        O(nlogn)'''
    
        if len(pos) == 0:
            return np.array([])
        # Converting a list is only efficient if the list is ordered
        # however want to be able to return list in the same order it
        # arrived, so remember the order and then sort.
        ordering = np.argsort(pos)
        pos = np.sort(pos)

        # pre allocate results list for speed
        try:
            results = np.zeros(len(pos))
        except TypeError:
            pos = np.array([pos])
            results = np.zeros(1)

        i = 0
        i_max = len(pos)

        # do linear search for correct exon
        for exon, interval in enumerate(self.transcript_intervals):

            while pos[i] < interval[1]:
                pos_within_exon = pos[i] - interval[0]
                genome_exon = self.genome_intervals[exon]
                relative_genome_position = genome_exon[0] + pos_within_exon

                if self.strand == "-":
                    results[i] = (self.offset - relative_genome_position)
                    i += 1
                else:
                    results[i] = (self.offset + relative_genome_position)
                    i += 1

                if i == i_max:
                    return results[ordering]
  
        # beyond the end of the transcript
        ValueError("Transcript postion %i outside of transcript %s" %
                   (pos[i], self.transcript_id))

    
def getCrosslink(read):
    ''' Finds the crosslinked base from a pysam read.

    Cross linked bases are definated as in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the iCLIP cDNAs mapped by Bowtie was used to
        define the cross link sites identified by truncated cDNAs.

        [For reads with deletions] The deleted nucleotide in CLIP and iCLIP
        cDNAs mapped by Novoalign was used to define the cross-link sites
        identified by read-through cDNAs. If a cDNA had more than one deletion,
        we selected the one closest to the beginning of the read.

    returns a tuple with the position of the read and one of the following
    categories:

        * truncated_neg
  
        * truncated_pos
      
        * deletion_neg

        * deletion_pos


    to record whether the position came from a truncation or a deletion '''


    if not 'D' in read.cigarstring:
        if read.is_reverse:
            pos = read.aend
            cat = "truncated_neg"
        else:
            pos = read.pos - 1
            cat = "truncated_pos"

    else:
        if read.is_reverse:
            cigar = reversed(read.cigar)
            position = find_first_deletion(cigar)
            pos = read.aend - position - 1
            cat = "deletion_neg"
        else:
            position = find_first_deletion(read.cigar)
            pos = read.pos + position
            cat = "deletion_pos"

    return (pos,cat)
    

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

    pos_depths = collections.defaultdict(int)
    neg_depths = collections.defaultdict(int)

 #   counter = E.Counter()

    for read in reads:
        
        (pos, cat) = getCrosslink(read)
 #       counter[cat] += 1

        if read.is_reverse:
            neg_depths[float(pos)] += 1
        else:
            pos_depths[float(pos)] += 1

    pos_depths = pd.Series(pos_depths, dtype=dtype)
    neg_depths = pd.Series(neg_depths, dtype=dtype)

    # check for integer overflow: counter sum should add up to array sum
    counter_sum = sum([y for x, y in counter.iteritems()])
    array_sum = pos_depths.sum() + neg_depths.sum()
    if not counter_sum == array_sum:
        raise (ValueError,
               "Sum of depths is not equal to number of "
               "reads counted, possibly dtype %s not large enough" % dtype)
    
  
        E.debug("Counted %i truncated on positive strand, %i on negative"
                % (counter.truncated_pos, counter.truncated_neg))
        E.debug("and %i deletion reads on positive strand, %i on negative"
                   % (counter.deletion_pos, counter.deletion_neg))
    
    return (pos_depths, neg_depths)


def count_intervals(bam, intervals, contig, strand=".", dtype='uint16'):
    ''' Count the crosslinked bases accross a transcript '''

    chr_len = bam.lengths[bam.gettid(contig)]
    exon_counts = []
    for exon in intervals:

        reads = bam.fetch(reference=contig,
                          start=exon[0],
                          end=exon[1],
                          reopen=False)

        count_results = countChr(reads, chr_len, dtype)

        # fetch pulls back any reads that *overlap* the specified coordinates
        # exlude Xlinked bases outside the interval (prevents double counting)

        if strand == "+":
            if len(count_results[0]) > 0:
                exon_counts.append(count_results[0].sort_index().loc[
                    float(exon[0]):float(exon[1]-1)])
        elif strand == "-":
            if len(count_results[1]) > 0:
                exon_counts.append(count_results[1].sort_index().loc[
                    float(exon[0]):float(exon[1]-1)])
            
        else:
            sum_counts = count_results[0].loc[(count_results[0].index >= exon[0]) &
                                              (count_results[0].index < exon[1])] + \
                count_results[1].loc[(count_results[1].index >= exon[0]) &
                                     (count_results[1].index < exon[1])]
            exon_counts.append(sum_counts)

    if len(exon_counts) == 0:
        transcript_counts = pd.Series()
    else:
        transcript_counts = pd.concat(exon_counts)
    # transcript_counts = transcript_counts.sort_index()
    return transcript_counts
        
