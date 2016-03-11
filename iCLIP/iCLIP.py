''' This a modeule that holds functions and classes useful for analysing iCLIP data '''

import CGAT.Experiment as E
import numpy as np
import pandas as pd
import CGAT.GTF as GTF
from CGAT.IndexedFasta import IndexedFasta
import collections
import itertools
import re
import pysam


AMBIGUITY_CODES = {'M': 'AC',
                   'R': 'AG',
                   'W': 'AT',
                   'S': 'CG',
                   'Y': 'CT',
                   'K': 'GT',
                   'V': 'ACG',
                   'H': 'ACT',
                   'D': 'AGT',
                   'B': 'CGT',
                   'N': 'CGAT'}

def IUPAC2Regex(sequence):

    for code, regex in AMBIGUITY_CODES.iteritems():
        sequence = re.sub(code, '[%s]' % regex, sequence)

    return sequence

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

        else:
            pos = read.pos - 1

    else:
        if read.is_reverse:
            cigar = reversed(read.cigar)
            position = find_first_deletion(cigar)
            pos = read.aend - position - 1

        else:
            position = find_first_deletion(read.cigar)
            pos = read.pos + position


    return pos
    

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

    counter = 0

    for read in reads:
        
        pos = getCrosslink(read)
        counter += 1

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
    array_sum = pos_depths.sum() + neg_depths.sum()
    if not counter == array_sum:
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


##################################################
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


##################################################
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


##################################################
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
        idx, counts = np.unique(randomised, return_counts=True)
        randomised = pd.Series(counts, index=idx).sort_index()
        return randomised


##################################################
def spread(profile, bases, reindex=True, right_bases=None):
       
    if right_bases:
        window = bases+right_bases
    else:
        window = 2*bases  + 1
        right_bases = bases

    if reindex:
        start = int(profile.index.min() - window)
        end = int(profile.index.max() + window+1)
        profile = profile.reindex(np.arange(start, end), fill_value=0)

    result = pd.rolling_sum(profile, window=window, center=False).dropna()
    result.index = result.index - (right_bases + 1)
    return result


##################################################
def corr_profile(profile1, profile2, nspread, profile2_ready=False):
    
    profile1 = profile1.reindex(
                    range(int(profile1.index.values.min())-1,
                          int(profile1.index.values.max())+1)).fillna(0)
    profile1 = spread(profile1, nspread, False)
        
    if not profile2_ready:
        profile2 = profile2.reindex(
                       range(int(profile2.index.values.min()),
                             int(profile2.index.values.max()))).fillna(0)
        profile2 = spread(profile2, nspread, False)
        
    return profile1.corr(profile2, method="spearman")


##################################################
def rand_apply(profile, exon, n, func, keep_dist=False,
               *args, **kwargs):
    '''Randomise a profile multiple times and apply a function
    to each randomised profile.

        :param profile: a profile with the number of reads at each base
        :type profile: pandas.Series
        :param exon: a GTF entry specifiying the boundaries to randomiseSites
                     between.
        :type exon: CGAT.GTF.Entry
        :param n: The number of randomisations to apply
        :param func: The function to apply to each randomisation
        :param keep_dist: Keep the read-per-base distribution constant

        :rtype: pandas.Series or pandas.DataFrame

    '''

    dummy = pd.Series(range(n))
 
    def _inner_func(x):
        rand = randomiseSites(profile, exon.start, exon.end,
                              keep_dist=keep_dist)
        return func(rand, *args, **kwargs)

    return dummy.apply(_inner_func)


##################################################
def bin_counts(counts, length, nbins):

    bins = np.linspace(0, length, num=nbins+1, endpoint=True)
    if len(counts.index.levels) == 2:
        bases = counts.index.droplevel()
    else:
        bases = counts.index

    binned_counts = counts.groupby(
        list(pd.cut(bases,
                    bins=bins,
                    labels=range(nbins),
                    include_lowest=True))).sum()

    binned_counts.index.name = "base"
    binned_counts = binned_counts.reindex(range(nbins), fill_value=0)
    return binned_counts


##################################################
def count_transcript(transcript, bam, flanks=0):
    '''Count clip sites from a bam and return a Series which transcript 
    domain coordinates.

        :type transcript: list of CGAT.GTF.Entry or similar
        :type bam: pysam.AlignmentFile
        :param flanks: If specified, flanks of this length are counted
                       and returned to the flank3 and flank5 section
                       of the now multiindexed return series
        :rtype: pandas.Series

    '''

    import pandas
    exons = GTF.asRanges(transcript, "exon")
    
    counts = count_intervals(bam, exons,
                             contig=transcript[0].contig,
                             strand=transcript[0].strand)
    
    coords_translator = TranscriptCoordInterconverter(transcript)
    
    counts.index = coords_translator.genome2transcript(counts.index)

    if flanks == 0:
        return counts
    else:

        if counts.sum() > 0:
            counts.index = pandas.MultiIndex.from_tuples(
                zip(*(['exons']*len(counts), counts.index)),
                names=["region", "base"])
        else:
            # deal with empty case
            counts.index = pandas.MultiIndex(levels=[["flank5", "exon", "flank3"],
                                                     []],
                                             labels=[[], []],
                                             names=["region", "base"])
 
        transcript_min = min([a for a, b in exons])
        transcript_max = max([b for a, b in exons])

        E.debug("Counting 5 flank: (%i, %i)" % (transcript_min - flanks,
                                                transcript_min))
        flank5_counts = count_intervals(bam, [(transcript_min - flanks,
                                              transcript_min)],
                                        contig=transcript[0].contig,
                                        strand=transcript[0].strand)
        E.debug(flank5_counts)
        E.debug("Counting 3 flank: (%i, %i)" % (transcript_max,
                                                transcript_max + flanks))
        flank3_counts = count_intervals(bam, [(transcript_max,
                                              transcript_max + flanks)],
                                        contig=transcript[0].contig,
                                        strand=transcript[0].strand)
        E.debug(flank3_counts)

        if transcript[0].strand == "-":
            index5 = ['flank3']
            index3 = ['flank5']
            bases5 = transcript_min - flank5_counts.index.values
            bases3 = transcript_max + flanks - flank3_counts.index.values

        else:
            index5 = ['flank5']
            index3 = ['flank3']
            bases5 = flank5_counts.index.values - transcript_min + flanks
            bases3 = flank3_counts.index.values - transcript_max

        if flank3_counts.sum() > 0:
            flank3_counts.index = pandas.MultiIndex.from_tuples(
                zip(*(index3*len(bases3), bases3)), names=["region", "base"])

        if flank5_counts.sum() > 0:
            flank5_counts.index = pandas.MultiIndex.from_tuples(
                zip(*(index5*len(bases5), bases5)), names=["region", "base"])

        E.debug("After direction detection and reindexing:")
        E.debug(flank3_counts)
        E.debug(flank5_counts)

        result = pandas.concat([flank5_counts, counts, flank3_counts])
        return result


##################################################
def meta_gene(gtf_filelike, bam, bins=[10, 100, 10], flanks=100,
              output_matrix=False, calculate_flanks=False, 
              pseudo_count=0):
    ''' Produce a metagene profile accross the :param gtf_file: from the reads
    in :param bam_file:.

        :type gtf_filelike: file or buffer
        :type bam: pysam.AlignmentFile
        :param output_matrix: if true, matrix of binned else None is returned

        :rtype: tuple of (pandas.Series, pandas.DataFrame)
    counts for each transcript is returned
    '''
    counts_collector = []
    regions = ["flank5", "exons", "flank3"]
    

    try:
        if flanks == 0:
            bins = bins[1]

        nbins_lookup = dict(zip(*(regions, bins)))
        length_lookup = dict(zip(*(regions, [flanks, 0, flanks])))

    except TypeError:
        if flanks > 0:
            nbins_lookup = {x: bins for x in bins}
   
    for transcript in GTF.transcript_iterator(
            GTF.iterator(gtf_filelike)):
      
        length = sum([x.end - x.start for x in transcript if x.feature == "exon"])

        if calculate_flanks:
            flanks = length * float(bins[0])/bins[1]
            length_lookup = dict(zip(*(regions, [flanks, 0, flanks])))

        counts = count_transcript(transcript, bam, flanks=flanks)

        if flanks > 0:
            length_lookup["exons"] = length
            binned_counts = counts.groupby(level=0).apply(
                lambda x: bin_counts(x, length_lookup[x.name], nbins_lookup[x.name]))
        else:
            binned_counts = bin_counts(counts, length, bins)

        binned_counts.name = transcript[0].transcript_id

        counts_collector.append(binned_counts)

    counts_matrix = pd.concat(counts_collector, axis=1)
    counts_matrix = counts_matrix.transpose()
    counts_matrix = (counts_matrix.fillna(0) + pseudo_count).div(counts_matrix.sum(axis=1), axis=0)
    summed_matrix = counts_matrix.sum()
    summed_matrix.name = "density"

    if output_matrix:
        return summed_matrix, counts_matrix
    else:
        return summed_matrix, None


##################################################
def find_all_matches(sequence, regexes):

    def _find_regex(regex):
        matches = regex.finditer(sequence)
        mStarts = [m.start() for m in matches]
        return np.asarray(mStarts)

    return regexes.apply(_find_regex)


def pentamer_frequency(profile, length, regex_matches, nSpread=15):
    '''Calculate the frequency of the each of a collection
    of sequence regexes on the provided read profile, and the
    coresponding sequence

        :param profile: A profile of the number of reads at each base
        :type profile: pandas.Series
        :param length: Length of the sequence represented by profile
        :param regex_matches: A pandas.Series of the locations of hits
                              for a set of regexs, as returned by 
                              find_all_matches
        :param nSpread: How far either side of each read to consider

        :rtype: pandas.Series with the count of each regex'''

    kmer = len(regex_matches.index.values[0])
    profile = profile.reindex(np.arange(-nSpread, length+nSpread-kmer), fill_value=0)
    profile = spread(profile, nSpread, False, nSpread - kmer)
    profile = profile.values


    def _count_regex(hits):

        if hits.size:
            return profile[hits].sum()
        else:
            return 0

    return regex_matches.map(_count_regex)
    

##################################################
LiteExon = collections.namedtuple('LiteExon', "start, end")


##################################################
def pentamer_enrichment(gtf_chunk_iterator, bam, fasta, kmer_length=5,
                        randomisations=100,
                        seperate_UTRs=False, spread=15,
                        pool=None):
    ''' This function calculates the enrichment of pentamers around
    CLIP'd sites

        :type gtf_chunk_iterator: an iterator that returns lists of CGAT.GTF
                                  Entries.
        :type bam: pysam.AlignmentFile
        :param seperate_UTRs: Treat UTRs as seperate areas from the rest of the
                              gene. Requires GTF file to contain CDS entries.
        :param spread: Number of bp around each base to use for sequence.
        :type fasta: CGAT.IndexedFasta
        :param pool: If present work will be parallelize across the worker pool
        :type pool: multiprocessing.Pool

        :rtype: pandas.Series with z-values for each pentamer
    '''

    bases = AMBIGUITY_CODES.keys()
    bases.remove("N")
    kmers = itertools.product('CGAT', repeat=kmer_length)
    kmers = ["".join(kmer) for kmer in kmers]
    regexs = [re.compile(kmer) for kmer in kmers]
    regexs = pd.Series(regexs, index=kmers)
    
    observed_kmer_counts = pd.Series(0, index=kmers)
    randomised_kmer_counts = pd.DataFrame(0,
                                          index=np.arange(randomisations),
                                          columns=kmers)

    args = ((profile, sequence,
             regexs, spread, randomisations)
            for profile, sequence in
            _get_counts_and_sequence(gtf_chunk_iterator, bam, fasta))



    if pool:
        results_iterator = pool.imap_unordered(_par_call, args)
    else:
        results_iterator = (_get_regex_frequencies(*arg) for arg in args)

    for observed, rands in results_iterator:

        observed_kmer_counts += observed
        randomised_kmer_counts += rands

    randomised_kmer_counts = randomised_kmer_counts.transpose()
    means = randomised_kmer_counts.mean(axis=1)
    sd = randomised_kmer_counts.std(axis=1)

    means.name = "mean"
    sd.name = "sd"
    observed_kmer_counts.name = "count"

    results = pd.concat([observed_kmer_counts, means, sd], axis=1)
    return results.apply(lambda x: (x['count'] - x['mean'])/x['sd'], 1)


def _get_counts_and_sequence(gtf_iterator, bam, fasta,
                             seperate_UTRs=False):
    '''Called by pentamer_enrichment. This function will return an iterator
    that yeilds tuples of profiles accross transcripts or introns and the
    sequence for which the profile is determined'''

    for transcript in gtf_iterator:

        E.debug("Counting transcript %s" % transcript[0].transcript_id)
        contig, strand = transcript[0].contig, transcript[0].strand

        # exons
        exons = GTF.asRanges(transcript, "exon")
        sequence = "".join(fasta.getSequence(contig, strand, exon[0], exon[1])
                           for exon in exons)
        exon_counts = count_transcript(transcript, bam)
        yield (exon_counts, sequence)

        # introns
        intron_intervals = GTF.toIntronIntervals(transcript)
        intron_counts = count_intervals(bam, intron_intervals, contig, strand)

        if intron_counts.sum() == 0:
            continue

        for intron in intron_intervals:
            
            seq = fasta.getSequence(contig, strand, intron[0], intron[1])
            profile = intron_counts.loc[float(intron[0]):float(intron[1])]
            profile.index = profile.index - intron[0]
            yield (profile, seq)

                
def _get_regex_frequencies(profile, sequence, regexs,
                           spread, randomisations):
    '''Called by pentamer_enrichment to get the frequencies
    counts for observed and randomisations. Allows
    parallelisation'''

    regex_matches = find_all_matches(sequence, regexs)
    length = len(sequence)
    boundaries = LiteExon(0, length)
    observed_kmer_counts = pentamer_frequency(profile,
                                              length,
                                              regex_matches,
                                              spread)
    randomised_kmer_counts = rand_apply(profile, boundaries,
                                        randomisations,
                                        pentamer_frequency,
                                        length=length,
                                        regex_matches=regex_matches,
                                        nSpread=spread)

    return observed_kmer_counts, randomised_kmer_counts


def _par_call(args):
    return _get_regex_frequencies(*args)
