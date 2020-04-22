import pandas as pd
import numpy as np
import collections
from functools import partial
from bx.bbi.bigwig_file import BigWigFile
import pysam
profile = collections.namedtuple("profile", ["centre",
                                             "read_end",
                                             "use_deletions",
                                             "reverse_direction",
                                             "offset",
                                             "filter_end"])
profiles = {"iclip": profile(centre=False,
                             read_end=False,
                             use_deletions=True,
                             reverse_direction=False,
                             offset=-1,
                             filter_end="none"),
            "eclip": profile(centre=False,
                             read_end=False,
                             use_deletions=True,
                             reverse_direction=False,
                             offset=-1,
                             filter_end="read2"),
            "iclip-centre":  profile(centre=True,
                                     read_end=False,
                                     use_deletions=True,
                                     reverse_direction=False,
                                     offset=-1,
                                     filter_end="read1"),
            "mNETseq-read1":  profile(centre=False,
                                      read_end=True,
                                      use_deletions=False,
                                      reverse_direction=False,
                                      offset=0,
                                      filter_end="read1"),
            "mNETseq-read2": profile(centre=False,
                                     read_end=False,
                                     use_deletions=False,
                                     reverse_direction=True,
                                     offset=0,
                                     filter_end="read2")}


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
def make_getter(bamfile=None, plus_wig=None, minus_wig=None, bedfile=None,
                profile="iclip", **kwargs):
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
    bedfile: str
        Name of a tabix indexed bed file to use to extract the crosslinks
    profile: str
        Name of a predefined set of options that refer to a particular 
        technique.  Currently implemented are "iclip" (default), "iclip-centre",
        "mNETseq-read1" and "mNETSeq-read2". Only applys when getter is made 
        from BAM.
    
    Common extra keyword arguments
    ------------------------------
    centre : bool
        Use the centre of the read rather than the base 5' of the read end.
        Only applicable if using a bamfile to retrieve the signal.
    read_end : bool
        Use 3' end of read rather than 5' end. 
        Only applicable if using a bamfile to retrieve the signal.
    use_deletions : bool
        Use a deletion in the read as the crosslink base.
        Only applicable if using a bamfile to retrieve the singal.
    reverse_direction : bool
        Reverse the strand of the tag so that reads on the + strand add counts
        to the - strand and vice versa (usually combined with filter_end="read2"
    offset : int
        Report the base offset by this much from the end of the read. E.g. in
        iCLIP report the base BEFORE the start of the read
    filter_end : None or str
        Filter reads so that only read1s or read2s are used.
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
    If kwargs are provided they override the settings from the profile


    ToDo
    ----
    Implement for fetching from tabix indexed bedGraph files.

'''

    if isinstance(profile, basestring):
        profile = profiles[profile]
        
    centre = kwargs.get("centre", profile.centre)
    read_end = kwargs.get("read_end", profile.read_end)
    use_deletions = kwargs.get("use_deletions", profile.use_deletions)
    reverse_direction = kwargs.get("reverse_direction",
                                   profile.reverse_direction)
    offset = kwargs.get("offset", profile.offset)
    filter_end = kwargs.get("filter_end", profile.filter_end)
    
    if bamfile is not None:
        if not isinstance(bamfile, pysam.AlignmentFile):
            bamfile = pysam.AlignmentFile(bamfile)
        return partial(_bam_getter, bamfile, centre=centre,
                       read_end=read_end, use_deletions=use_deletions,
                       reverse_strand=reverse_direction,
                       offset=offset,
                       filter_end=filter_end)
    elif plus_wig is not None:
        plus_wig = BigWigFile(open(plus_wig))
        if minus_wig is not None:
            minus_wig = BigWigFile(open(minus_wig))

        return partial(_wig_getter, plus_wig=plus_wig, minus_wig=minus_wig)
    elif bedfile is not None:
        if not isinstance(bedfile, pysam.TabixFile):
            bedfile = pysam.TabixFile(bedfile)
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
        counts = plus_wig.get_as_array(contig, int(start), int(end))
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
                centre=False, read_end=False, use_deletions=True,
                reverse_strand=False, offset=-1, filter_end=None):
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

    if filter_end == "read1":
        reads = (r for r in reads if not r.is_read2)
    elif filter_end == "read2":
        reads = (r for r in reads if not r.is_read1)
        
    counts = countChr(reads, chr_len, dtype, centre, read_end, use_deletions,
                      offset)

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

    if reverse_strand:
        positive = "-"
        negative = "+"
    else:
        positive = "+"
        negative = "-"
        
    if strand == positive:
        return counts[0]
    elif strand == negative:
        return counts[1]
    elif strand == ".":
        return counts[0].add(counts[1], fill_value=0)
    elif strand == "both":
        if reverse_strand:
            return (counts[1], counts[0], counts[2])
        else:
            return counts
    else:
        raise ValueError("Invalid strand")

        
##################################################    
def _bed_getter(bedfile, contig, start=0, end=None, strand=".", dtype="uint16"):
    '''Get crosslink profiles from tabix indexed bedGraph/Bed'''

    # check the file contains some data for the requested contig
    if not contig in bedfile.contigs:
        #print "%s not in bedfile" % contig
        return pd.Series(dict(), dtype=dtype)
    
    # fetch the rercords from the specificed region
    crosslinks = bedfile.fetch(contig, start, end, parser=pysam.asBed())
    
    profile = dict()

    check_sum = 0
    
    for base in crosslinks:
        try:
            correct_strand = strand == "." or base.strand == strand
        except AttributeError:
            correct_strand = True
        except KeyError:
            correct_strand = True
            
        if correct_strand:

            try:
                profile[float(base.start)] = int(base.score)
                check_sum += int(base.score)
            except AttributeError:
                profile[float(base.start)] = 1
                check_sum += 1
            except KeyError:
                profile[float(base.start)] = 1
                check_sum += 1

    if len(profile.keys())==0:
        profile = pd.Series(profile, dtype=dtype, index=pd.Index([], dtype="float"))
    else:
        profile = pd.Series(dict(profile), dtype=dtype)

            
    #if not check_sum == profile.sum():
    #    raise OverflowError("Check sum failed (%i = %i). Possibly counts exceed specified dtype. Use bigger dtype"
    #                        % (check_sum, profile.sum()))

    return profile
 
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
def getCrosslink(read, centre=False, read_end=False, use_deletions=True, offset=-1):
    '''Finds the crosslinked base from a pysam read.

    Parameters
    ----------
    read : pysam.AlignedSegment
    centre : bool
             Use centre of read rathter than usual Sugimoto rules
    read_end : bool
        Use the end of the read rather than the start
    use_deletions : bool
        If a deletion is present in the gene, use it as the crosslinks
        base. 
             
    Returns
    -------
    int 
         Position of crosslink in 0-based genome coordinates.
         
    Notes
    -----
    
    If the read contains at least one deletion base (and use_deletions
    is True), as identified by the cigar string, then the position of
    the 5' most deleted base is used.
    
    If no deletion is present then cross-linked bases are defined in three
    ways:
        
    1.  As in Sugimoto et al, Genome Biology 2012

        The nucleotide preceding the start of the read is 
        used. Read strand is taken into account, such that for reads
        with `is_reverse=True`, the base after `aend` is used rather 
        than the base before `pos`.

    2.  (if `read_end==True`) As 1, but using the end of the read rather
        than the start.

    3.  (if `centre=True`) returns the centre base of the read,
        accounting for splicing etc

    '''

    reverse_direction = (read.is_reverse and not read_end) or \
                        (not read.is_reverse and read_end)
    
    if  not use_deletions or 'D' not in read.cigarstring:

        if centre:
            reference_bases = read.get_reference_positions(full_length=True)
            i = len(reference_bases)/2
            while reference_bases[i] is None and i > 0:
                i = i -1
            return reference_bases[i]

        if reverse_direction:
            pos = read.aend - 1 - offset

        else:
            pos = read.pos + offset

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
def countChr(reads, chr_len, dtype='uint16', centre=False,
             read_end=False, use_deletions=True, offset=-1):
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
        
        pos = getCrosslink(read, centre, read_end, use_deletions, offset)
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
        raise (ValueError(
               "Sum of depths is not equal to number of "
               "reads counted, possibly dtype %s not large enough" % dtype))
    
#    E.debug("Counted %i truncated on positive strand, %i on negative"
#            % (counter.truncated_pos, counter.truncated_neg))
#    E.debug("and %i deletion reads on positive strand, %i on negative"
#            % (counter.deletion_pos, counter.deletion_neg))
    
    return (pos_depths, neg_depths, counter)
