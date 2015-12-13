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
       
        raise ValueError("Position %i (%i relative) is not in transcript %s\n exons are %s" %
                         (pos[i], relative_pos[i], self.transcript_id, self.genome_intervals))

    def transcript2genome(self, pos):
        ''' Convert transcript coodinate into genome coordinate,
        pos can be a single value or a nunpy array like object.
        Passing an array ensures that the transcript is only
        searched once, ensuring O(n) performance rather than
        O(nlogn)'''
    
        try:
            if len(pos) == 0:
                return np.array([])
        except TypeError:
            pos = np.array([pos])

        # Converting a list is only efficient if the list is ordered
        # however want to be able to return list in the same order it
        # arrived, so remember the order and then sort.
        ordering = np.argsort(pos)
        pos = np.sort(pos)

        # pre allocate results list for speed
       
        results = np.zeros(len(pos))
        
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

    def transcript_interval2genome_intervals(self, interval):
        '''Take an interval in transcript coordinates and returns
        a list of intervals in genome coordinates representing the
        interval on the genome '''

        outlist = []
        for exon in self.transcript_intervals:
            if interval[0] < exon[1]:
                start = interval[0]
                if interval[1] <= exon[1]:
                    outlist.append((start, interval[1]))
                    break
                else:
                    outlist.append((start, exon[1]))
                    interval = (exon[1], interval[1])
       
        genome_list = [tuple(self.transcript2genome((x, y-1))) for
                       x, y in outlist]
        
        # these intervals are zero based-closed. Need to make half open

        if self.strand == "+":
            genome_list = [(x, y+1) for x, y in genome_list]
        else:
            genome_list = [(y, x+1) for x, y in genome_list]

        return sorted(genome_list)



    
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

    counter = E.Counter()

    for read in reads:
        
        (pos, cat) = getCrosslink(read)
        counter[cat] += 1

        if read.is_reverse:
            neg_depths[float(pos)] += 1
        else:
            pos_depths[float(pos)] += 1

    try:
        pos_depths = pd.Series(pos_depths, dtype=dtype)
    except ValueError:
        pos_depths = pd.Series({}, dtype=dtype)

    try:
        neg_depths = pd.Series(neg_depths, dtype=dtype)
    except ValueError:
        neg_depths = pd.Series({}, dtype=dtype)
 
    # check for integer overflow: counter sum should add up to array sum
    counter_sum = sum([y for x, y in counter.iteritems()])
    array_sum = pos_depths.sum() + neg_depths.sum()
    if not counter_sum == array_sum:
        raise (ValueError,
               "Sum of depths is not equal to number of "
               "reads counted, possibly dtype %s not large enough" % dtype)
    
#    E.debug("Counted %i truncated on positive strand, %i on negative"
#            % (counter.truncated_pos, counter.truncated_neg))
#    E.debug("and %i deletion reads on positive strand, %i on negative"
#            % (counter.deletion_pos, counter.deletion_neg))
    
    return (pos_depths, neg_depths, counter)


def count_intervals(bam, intervals, contig, strand=".", dtype='uint16'):
    ''' Count the crosslinked bases accross a transcript '''

    chr_len = bam.lengths[bam.gettid(contig)]
    exon_counts = []
    for exon in intervals:
        
        # X-linked position is first base before read: need to pull back
        # reads that might be one base out. Extra bases will be filtered out
        # later.
        try:
            reads = bam.fetch(reference=contig,
                              start=max(0,exon[0]-1),
                              end=exon[1]+1)
        except ValueError as e:
            E.debug(e)
            E.warning("Skipping intervals on contig %s as not present in bam"
                      % contig)
            return pd.Series()

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
        
def calcAverageDistance(profile1, profile2):
    ''' This function calculates the average distance of all
    pairwise distances in two profiles'''

    def _cartesian(x, y):
        return np.transpose([np.tile(x, len(y)), np.repeat(y, len(x))])

    positions = _cartesian(profile1.index.values, profile2.index.values)
    counts = _cartesian(profile1.values, profile2.values)
    counts = np.prod(counts, axis=1)
    distances = np.abs(positions[:, 0] - positions[:, 1])
    mean_distance = (distances.astype("float64") * counts).sum() / np.sum(counts)

    return mean_distance


def findMinDistance(profile1, profile2):
    '''Finds mean distance between each read in profile1
    and a read in profile2'''

    locations1 = profile1.index.values

    locations2 = profile2.index.values # .astype("int16")

    mat1 = np.repeat(locations1, locations2.size).reshape(
        (locations1.size, locations2.size))
    mat2 = np.tile(locations2, locations1.size).reshape(
        (locations1.size, locations2.size))

    distances = np.abs(mat1-mat2).min(axis=1)

    return distances.mean()

def randomiseSites(profile, start, end, keep_dist=True):
    '''Randomise clipped sites within an interval (between start and end)
    if keep_dist is true, then reads on the same base are kept togehter'''

    if keep_dist:

        profile = profile.copy()
        profile.index = np.random.choice(
            np.arange(start, end), profile.size, replace=False)
        profile = profile.sort_index()
        return profile

    else:
        randomised = np.random.choice(
            np.arange(start, end), profile.sum(), replace=True)
        randomised = pd.Series(randomised).value_counts().sort_index()
        return randomised

def spread(profile, bases, reindex=True):
   
    start = int(profile.index[0] - 2*bases)
    end = int(profile.index[-1] + 2*bases+1)
    
    if reindex:
        profile = profile.reindex(range(start, end))
        profile = profile.fillna(0)

    return pd.rolling_sum(profile, window=2*bases+1, center=True).dropna()

def corr_profile(profile1, profile2, nspread, profile2_ready=False):
    
    profile1 = profile1.reindex(
                    range(int(profile1.index.values.min())-1,
                          int(profile1.index.values.max())+1)).fillna(0)
    profile1 = spread(profile1,nspread, False)
        
    if not profile2_ready:
       profile2 = profile2.reindex(
                       range(int(profile2.index.values.min()),
                             int(profile2.index.values.max()))).fillna(0)
       profile2 = spread(profile2, nspread, False)
        
    return profile1.corr(profile2, method="spearman")

def rand_apply(profile, exon, n, func, keep_dist=False, *args, **kwargs):
    dummy = pd.Series(range(n))
    def _inner_func(x):
        rand = randomiseSites(profile, exon.start, exon.end,
                                    keep_dist=keep_dist)
        return func(rand, *args, **kwargs)
    return dummy.apply(_inner_func)



