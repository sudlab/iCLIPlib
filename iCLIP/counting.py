''' These functions deal with calculating the read counts within a genomic region. 
in general they return a pandas Series of read counts. '''


import pandas as pd
import numpy as np
import collections
from functools import partial
from bx.bbi.bigwig_file import BigWigFile
import pysam

import cgatcore.experiment as E
import cgat.GTF as GTF

from .utils import TranscriptCoordInterconverter

from .getters import make_getter





##################################################
def count_intervals(getter, intervals, contig, strand=".", dtype='uint32',
                    use_centre=False, get_all_and_mask=True):
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
    get_all_and_mask : bool, optional
        Instead of fetching each exon indevidually from the getter,
        fetch entire chromsome region, and then subset to things
        in the interval.

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

    if len(intervals)==0:
        return pd.Series()
    
    if isinstance(getter, pysam.AlignmentFile):
        getter = make_getter(bamfile=getter, centre=use_centre)

    exon_counts = []
    
    if get_all_and_mask:
        start = min(float(exon[0]) for exon in intervals)
        end = max(float(exon[1]) for exon in intervals)
        
        all_counts = getter(contig=contig,
                            start=start,
                            end=end,
                            strand=strand,
                            dtype=dtype)

        if all_counts.sum() == 0:
            return pd.Series()
        
        for exon in intervals:
            exon_counts.append(all_counts.loc[exon[0]:(exon[1]-1)])

    else:
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

    if len(exons)==0:
        if flanks==0:
            return pandas.Series()
        else:
            return pandas.Series(index=pandas.MultiIndex(
                levels=[["flank5", "exon", "flank3"],
                        []],
                codes=[[], []],
                names=["region", "base"]))

    counts = count_intervals(bam, exons,
                             contig=transcript[0].contig,
                             strand=transcript[0].strand)
    
    coords_translator = TranscriptCoordInterconverter(transcript)
    
    counts.index = coords_translator.genome2transcript(counts.index)

    if flanks == 0:
        return counts
    else:

        if counts.sum() != 0:
            counts.index = pandas.MultiIndex.from_tuples(
                zip(*(['exons']*len(counts), counts.index)),
                names=["region", "base"])

        else:
            # deal with empty case. Cannot just return an empty series with the
            # correct index, as this screws up any group-bys that need to
            # run for every group (even if empty).
            counts = pandas.Series(0)
            counts.index = pandas.MultiIndex.from_tuples(
                [('exons', 1)],
                names=["region", "base"])
            #counts.index = pandas.MultiIndex(
            #    levels=[["flank5", "exon", "flank3"],
            #            []],
            #    codes=[[], []],
            #    names=["region", "base"])
 
        transcript_min = min([a for a, b in exons])
        transcript_max = max([b for a, b in exons])

        E.debug("Counting 5 flank: (%i, %i)" % (transcript_min - flanks,
                                                transcript_min))
        flank5_counts = count_intervals(bam, [(max(transcript_min - flanks,0),
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
            # deal with empty case - see counts case above
            flank3_counts = pandas.Series(0)
            flank3_counts.index = pandas.MultiIndex.from_tuples(
                [(index3[0], 1)],
                names=["region", "base"])
            #flank3_counts.index = pandas.MultiIndex(
            #    levels=[["flank5", "exon", "flank3"],
            #            []],
            #    codes=[[], []],
            #    names=["region", "base"])

        if flank5_counts.sum() > 0:
            flank5_counts.index = pandas.MultiIndex.from_tuples(
                zip(*(index5*len(bases5), bases5)), names=["region", "base"])
        else:
            # deal with empty case - see counts case above. 
            flank5_counts = pandas.Series(0)
            flank5_counts.index = pandas.MultiIndex.from_tuples(
                [(index5[0], 1)],
                names=["region", "base"])
            #flank5_counts.index = pandas.MultiIndex(
            #    levels=[["flank5", "exon", "flank3"],
            #            []],
            #    codes=[[], []],
            #    names=["region", "base"])

        E.debug("After direction detection and reindexing:")
        E.debug(flank3_counts)
        E.debug(flank5_counts)

        result = pandas.concat([flank5_counts, counts, flank3_counts])
        return result

