import pandas as pd
import numpy as np
import collections
from functools import partial
from bx.bbi.bigwig_file import BigWigFile
import pysam

from counting import countChr

def getter(contig, start=0, end=None, strand=".", dtype="uint16"):
    '''Get the profile of crosslinks across a genomic region, returning
    the number of crosslinks on each base.

    This stub defines the abstract interface for a getter function.
    Getters are made by calling the :func:`make_getter` function. The
    returned function will have the interface specified here. 

    Parameters
    ----------
    contig : str
    start : int, optional
    end : int or None, optional
        If None specified then profile will be returned from 0 to end of
        contig. 
    strand : ['.','+','-'], optional
        strand to return counts from. If ``.`` then data from both strands
        is returned.
    dtype : str
       dtype to use for storing data. Usually a type of unsigned int. Larger
       size ints are less likely to overflow, but take more memeory.

    Returns
    -------
    pandas.Series
        Profile of crosslinks. Index is the genomic base, value is the count.
        Series is sparse, that is bases with no crosslinks are excluded from
        the index.

    See also
    --------
    make_getter : genearate a getter function from a bamfile or BigWig
        file(s)

    '''
    raise NotImplementedError("This is an abstract function. Please"
                              " create an instance of a concrete version "
                              "using make_getter")

                              
##################################################
def make_getter(bamfile=None, plus_wig=None, minus_wig=None, centre=False):
    ''' A factory for getter functions

    Parameters
    ----------
    bamfile : str or pysam.AlignmentFile
        A bamfile to use to extract crosslinks
    plus_wig : str
        Name of a BigWig file to use as either the unstranded signal, or 
        the plus stranded signal. If no `minus_wig` is provided signal is
        assumed to be unstranded
    minus_wig : str
        Name of a BigWig file to use as the minus strand signal
    centre : bool
        Use the centre of the read rather than the base 5' of the read end.
        Only applicable if using a bamfile to retrieve the signal.

    Returns
    -------
    func
        A :func:`getter` function that will retrieve the profile of crosslinks
        over a specified genomic interval.

    See Also
    --------
    getter : Abstract interface for returning crosslink profiles over genomic
    regions.

    
    Notes
    -----
    
    At least one or `bamfile` or `plus_wig` must be provided.

    ToDo
    ----
    Implement for fetching from tabix indexed bedGraph files.

'''

    if bamfile is not None:
        if not isinstance(bamfile, pysam.AlignmentFile):
            bamfile = pysam.AlignmentFile(bamfile)
        return partial(_bam_getter, bamfile, centre=centre)
    elif plus_wig is not None:
        plus_wig = BigWigFile(open(plus_wig))
        if minus_wig is not None:
            minus_wig = BigWigFile(open(minus_wig))

        return partial(_wig_getter, plus_wig=plus_wig, minus_wig=minus_wig)
    elif bedfile is not None:
        if not isinstance(bedfile, pysam.Tabix):
            bedfile = pysam.Tabix(bedfile)
        return partial(_bed_getter, bedfile=bedfile)
    else:
        raise ValueError("Please provide either a bamfile, a bigwig or a"
                         "befile")        
        

##################################################
def _wig_getter(plus_wig, minus_wig, contig, start=0, end=None,
               strand=".", dtype="unit16"):
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
def _bam_getter(bamfile, contig, start=0, end=None, strand=".", dtype="uint16",
               centre=False):
    '''A function to get iCLIP coverage across an interval from a BAM file'''
    chr_len = bamfile.lengths[bamfile.gettid(contig)]
    if end is None:
        end = chr_len

    try:
        reads = bamfile.fetch(reference=contig,
                              start=max(start-1, 0),
                              end=end+1)
    except ValueError as e:
        E.warning(e)
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
def _bed_getter(bedfile, contig, start=0, end=None, strand=".", dtype="uint16"):
    '''Get crosslink profiles from tabix indexed bedGraph/Bed'''
    
    # fetch the rercords from the specificed region
    crosslinks = bedfile.fetch(contig, start, end, pysam.asBed())
    
    profile = pandas.Series(dtype=dtype)
    
    for base in crosslinks:
        try:
            correct_strand = strand == "." or base.strand == strand
        except AttributeError:
            correct_strand = True
            
        if correct_strand:
            profile[base] = int(base.score)
            check_sum += int(base.score)

    if not check_sum == profile.sum():
        raise OverflowError("Counts exceed specified dtype. Use bigger dtype")

    return profile
 