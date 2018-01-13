''' These functions deal with calculating the read counts within a genomic region. 
in general they return a pandas Series of read counts. '''


import pandas as pd
import numpy as np
import collections
from functools import partial
from bx.bbi.bigwig_file import BigWigFile
import pysam

import CGAT.Experiment as E
import CGAT.GTF as GTF

from utils import TranscriptCoordInterconverter


def find_first_deletion(cigar):
    '''Find the position of the the first deletion in a
    cigar string.

    Parameters
    ----------
    cigar : tuple of tuple of int
            a cigar string as returned by pysam.AlignedSegment.cigartuples

    Returns
    -------
    int
         The position of the first deletion in the cigar, returns 0 if no
         deletion found.

    '''

    position = 0
    for operation, length in cigar:

        if operation == 2:
            return position
        else:
            position += length

    return position


##################################################
def getCrosslink(read, centre=False):
    '''Finds the crosslinked base from a pysam read.

    Parameters
    ----------
    read : pysam.AlignedSegment
    centre : bool
             Use centre of read rathter than usual Sugimoto rules
             
    Returns
    -------
    int 
         Position of crosslink in 0-based genome coordinates.
         
    Notes
    -----
    
    If the read contains at least one deletion base, as identified
    by the cigar string, then the position of the 5' most deleted
    base is used.
    
    If no deletion is present then cross-linked bases are defined in two
    ways:
        
    1.  As in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the start of the read is 
        used. Read strand is taken into account, such that for reads
        with `is_reverse=True`, the base after `aend` is used rather 
        than the base before `pos`.

    2.  (if `centre=True`) returns the centre base of the read,
        accounting for splicing etc

    '''

    if 'D' not in read.cigarstring:

        if centre:
            reference_bases = read.get_reference_positions(full_length=True)
            i = len(reference_bases)/2
            while reference_bases[i] is None and i > 0:
                i = i -1
            return reference_bases[i]

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


##################################################
def countChr(reads, chr_len, dtype='uint16', centre=False):
    ''' Counts the crosslinked bases for each provided read.

    Scans through the provided pysam.rowiterator reads and tallys the
    number of cross links on each of the positive and negative strand. 
    
    Parameters
    ----------
    reads : iter of reads
        Reads to get cross-links from. Most often a pysam.rowiterator
        returned by pysam.AlignmentFile.fetch
    chr_len : int
        Length of contig/region being counted
    dtype : str
        dtype to use to store the counts. 
    centre : bool
        Use centre of read for cross-link rather than base preceding 5'
        end
        
    Returns
    -------
    positive_counts, negative_counts : pandas.Series
        Counts of cross-linked bases on positive and negative strands. The
        Series are indexed on genome position and are sparse (see below)
    count : 
        total count of cross-linked bases
    
    See also
    --------
    getCrosslink : gets cross-link position for pysam.AlignedSegment
    
    Notes
    -----  
    The `dtype` to use internally for storage can be specified. Large types
    reduce the chance of overflow, but require more memory. With 'uint16'
    the largest count that can be handled is 255. Data is stored sparse,
    so memory is less of a problem. Overflow will cause a ValueError.
    
    See :func:`getCrosslink` for a description of how the cross-link 
    position is determined from each read.

    '''

    pos_depths = collections.defaultdict(int)
    neg_depths = collections.defaultdict(int)

    counter = 0

    for read in reads:
        
        pos = getCrosslink(read, centre)
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


##################################################
def wig_getter(plus_wig, minus_wig, contig, start=0, end=None,
               dtype="uint16", strand="."):
    ''' Get depth over intervals from a bigWig file. Requests for
    intervals starting < 0 will be truncated to 0'''

    if start < 0:
            start = 0
            E.warning("Truncating start of interval  %s:%i-%i"
                      % (contig, start, end))

    if strand == "+" or minus_wig is None:
        counts = plus_wig.get_as_array(contig, start, end)
        result = pd.Series(
            counts, index=np.arange(start, end, dtype="float")).dropna()
        return result
 
    elif strand == "-":
        counts = -1 * minus_wig.get_as_array(contig, start, end)
        result = pd.Series(
            counts, index=np.arange(start, end, dtype="float")).dropna()
        return result

    elif strand == ".":
        plus_counts = plus_wig.get_as_array(contig, start, end)
        minus_counts = -1 * minus_wig.get_as_array(contig, start, end)
        plus_result = pd.Series(
            plus_counts, index=np.arange(start, end, dtype="float")).dropna()
        minus_result = pd.Series(
            minus_counts, index=range(start, end)).dropna()
        result = plus_result.add(minus_result, fill_value=0)
        return result


##################################################    
def bam_getter(bamfile, contig, start=0, end=None, strand=".", dtype="uint16",
               centre=False):
    '''A function to get iCLIP coverage across an interval from a BAM file'''
    chr_len = bamfile.lengths[bamfile.gettid(contig)]
    if end is None:
        end = chr_len

    try:
        reads = bamfile.fetch(reference=contig,
                              start=start-1,
                              end=end+1)
    except ValueError as e:
        E.debug(e)
        E.warning("Skipping intervals on contig %s as not present in bam"
                  % contig)
        return pd.Series()

    counts = countChr(reads, chr_len, dtype, centre)

    # Two sets of extrainous reads to exlucde: firstly we have pull back
    # reads with a 1bp extra window. Second fetch pulls back overlapping
    # reads, not fully contained reads.

    counts = list(counts)
    if strand != "-" and len(counts[0]) > 0:
        counts[0] = counts[0].sort_index().loc[
            float(start):float(end-1)]
    if strand != "+" and len(counts[1]) > 0:
        counts[1] = counts[1].sort_index().loc[
            float(start):float(end-1)]
    
    if strand == "+":
        return counts[0]
    elif strand == "-":
        return counts[1]
    elif strand == ".":
        return counts[0].add(counts[1], fill_value=0)
    else:
        raise ValueError("Invalid strand")


##################################################
def make_getter(bamfile=None, plus_wig=None, minus_wig=None, centre=False):
    ''' A factory for getter functions '''

    if bamfile is not None:
        if not isinstance(bamfile, pysam.AlignmentFile):
            bamfile = pysam.AlignmentFile(bamfile)
        return partial(bam_getter, bamfile, centre=centre)
    else:
        plus_wig = BigWigFile(open(plus_wig))
        if minus_wig is not None:
            minus_wig = BigWigFile(open(minus_wig))

        return partial(wig_getter, plus_wig=plus_wig, minus_wig=minus_wig)


##################################################
def count_intervals(getter, intervals, contig, strand=".", dtype='uint32',
                    use_centre=False):
    ''' Count the crosslinked bases across a set of intervals

    Parameters
    ----------
    getter : *_getter function or pysam.AlignmentFile
        The getter function to retrieve the cross-links from. Returned from
        make_getter.
    intervals : iter of tuples
        The intervals to retrieve cross-links from as an iterable of
        (start, end) tuples.
    contig : str
    strand : {"+", "-", "."}
        If "." then the sum of counts on both strands is returned.
    dtype : str
        `numpy` dtype used for storing counts. Larger dtypes reduce chance
        of overflow, but use more memory. Defaults to 'uint32'.
    
    Returns
    -------
    pandas.Series
        Series of integer counts, where the index of the Series is the
        genome coordinates. Series is spares (i.e. 0 count bases are 
        excluded) 
        
    See also
    --------
    make_getter : Access to cross-link data for range of file-types
    getCrosslink : find cross-link base from pysam.AlignedSegment
    '''

    if isinstance(getter, pysam.AlignmentFile):
        getter = make_getter(bamfile=getter, centre=use_centre)

    exon_counts = []
    for exon in intervals:
        exon_counts.append(getter(contig=contig,
                                  start=exon[0],
                                  end=exon[1],
                                  strand=strand,
                                  dtype=dtype))

    if len(exon_counts) == 0:
        transcript_counts = pd.Series()
    else:
        transcript_counts = pd.concat(exon_counts)
    # transcript_counts = transcript_counts.sort_index()

    return transcript_counts


##################################################
def count_transcript(transcript, bam, flanks=0):
    '''Count clip cross-link sites across a transcript
    
    Provided with a complete transcript, count_transcript returns all 
    cross-links in the transcript and returns counts in transcript domain
    coordinates.
    
    Parameters
    ----------
    transcript : iter of `CGAT.GTF.Entry` objects
        The transcript to be counted over. Usually as returned by one of
        `CGAT.GTF` iterators such as `CGAT.GTF.transcript_iterator` or
        `CGAT.GTF.flat_gene_iterator`
    bam : *_getter function
        getter function to retrieve crosslinks from, as returned by 
        :func:`make_getter`.
    flanks : int
        Length of 5' and 3' flanks in bp to count on either end of the 
        transcript. Note: changes index type of returned Series.
        
    Returns
    -------
    pandas.Series
        Series of counts. Inner-most level of the index is the coordinates
        in the transcript domain (see below). If `flanks` > 0 then 
        will have MultiIndex (see below).
        
    See also
    --------
    make_getter : Access to cross-link data for range of file-types.
    getCrosslink : find cross-link base from pysam.AlignedSegment.
    TranscriptCoordInterconverter : Conversion between genome and
        transcript coordinates.
    Notes
    -----
    
    Returns coordinates in the transcript domain. That is the base
    corresponding to the TSS of the transcript will be base 0 and the 
    bases are numbered 0 to sum(length of exons). Introns are excluded both
    from the counting and from the coordinates. Strand is also accounted
    for. 
    
    This function can also return flanking regions upstream and downstream
    of the transcript. This is specified by `flanks`. If requested, the
    returned Series will have a `pandas.MultiIndex`. The inner level will
    correspond to the base and the outer level to whether the inner level
    refers to the the 5' flank ('flank5'), the 3' flank ('flank3') or the 
    transcript ('exon')

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
            counts.index = pandas.MultiIndex(
                levels=[["flank5", "exon", "flank3"],
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
        else:
            # deal with empty case
            flank3_counts.index = pandas.MultiIndex(
                levels=[["flank5", "exon", "flank3"],
                        []],
                labels=[[], []],
                names=["region", "base"])

        if flank5_counts.sum() > 0:
            flank5_counts.index = pandas.MultiIndex.from_tuples(
                zip(*(index5*len(bases5), bases5)), names=["region", "base"])
        else:
            # deal with empty case
            flank5_counts.index = pandas.MultiIndex(
                levels=[["flank5", "exon", "flank3"],
                        []],
                labels=[[], []],
                names=["region", "base"])

        E.debug("After direction detection and reindexing:")
        E.debug(flank3_counts)
        E.debug(flank5_counts)

        result = pandas.concat([flank5_counts, counts, flank3_counts])

        return result

