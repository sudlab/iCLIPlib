import re
import numpy as np
import pandas as pd
import cgat.GTF as GTF


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


##################################################
def IUPAC2Regex(sequence):

    for code, regex in AMBIGUITY_CODES.items():
        sequence = re.sub(code, '[%s]' % regex, sequence)

    return sequence


##################################################
class TranscriptCoordInterconverter:
    '''Interconvert between genome domain and transcript domain 
    coordinates.
    
    As there are expected to be many calls against the same transcript,
    time can be saved by precomputation. Overlapping exons are merged.

    Parameters
    ----------
    transcript : sequence of `CGAT.GTF.Entry`
        Set of GTF entires representing a transcript.
    introns : bool, optional
        Use introns instead of exons (see below).

    Attributes
    ----------
    transcript_id : str
        Value of the transcript_id field of the transcript in quesiton.
        Taken from the transcript_id field of the first entry in 
        `transcript`.
    strand : str
        Strand of the transcript. Taken from the strand field of the first
        entry in `transcript`.
    offset : int
        Position of the start of the transcript in genome coordinates
    genome_intervals : list of tuples of int
        Coordinates of exons (or introns) in genome space as the difference
        from `offset`. These are sorted in transcript order (see below)
    transcript_intervals : list of tuples of int  
        Coordinates of exons (or introns) in transcript space. That is
        absolute distance from transcription start site after splicing.
    length : int
        Total length of intervals (exons or introns) in the transcript

    Notes
    -----

    Imagine the following transcript::

        chr1  protein_coding  exon  100  108   .  -  .  transcript_id "t1"; gene_id "g1";
        chr1  protein_coding  exon  112  119   .  -  .  transcript_id "t1"; gene_id "g1";
        chr1  protein_coding  exon  100  108   .  -  .  transcript_id "t1"; gene_id "g1";      

    We can visualise the relationship between the different coordinate
    domains as below::

        Genome coordinates:    1         1         1         1
                               0         1         2         3       
                               0123456789012345678901234567890
        Transcript:            |<<<<<<|----|<<<<<|-----|<<<<<|
        Transcript Coords:      2             1              0
                               10987654    3210987     6543210
          with `introns=True`:         8765       43210   

    Thus the intervals representing the exons in the transcript domain are
    (0, 7), (7,14), (14, 22), and the genome base 115 corresponds to
    transcript base 10. 

    

    TranscriptCoordInterconverter.genome2transcript should be the
    interverse of TranscriptCoordInterconverter.transcript2genome.

    That is if::

        myConverter = TranscriptCoordInterverter(transcript)
    
    then::
    
        myConverter.genome2transcript(myConverter.transcript2genome(x)) == x
    
    and::

        myConverter.transcript2genome(myConverter.genome2transcript(x)) == x
    '''

    def __init__(self, transcript, introns=False):
        ''' Pre compute the conversions for each exon '''

        if not introns:
            intervals = GTF.asRanges(transcript, feature="exon")
        else:
            intervals = GTF.toIntronIntervals(transcript)
        
        # get strand
        self.strand = transcript[0].strand

        # store transcript_id
        try:
            self.transcript_id = transcript[0].transcript_id
        except AttributeError:
            self.transcript_id = transcript[0].gene_id
            
        # sort the exons into "transcript" order
        if self.strand == "-":
            intervals.sort(reverse=True)
            intervals = [(y-1, x-1) for x, y in intervals]
        else:
            intervals.sort(reverse=False)

        self.offset = intervals[0][0]
        self.genome_intervals = [list(map(abs, (x-self.offset, y-self.offset)))
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
        '''Convert genome coordinate into transcript coordinates.

        Can be a single coordinate or a array-like of coordinates

        Parameters
        ----------
        pos : int or array-like or int
            positions, in the genome domain, to be converted
        
        Returns
        -------
        int or numpy.array
            The position, or positions converted into the transcript
            domain.

        Raises
        ------
        ValueError
            If the supplied genome position is not in the transcript.
            This could be because it falls into one of the introns
            (or exons if the converter was created with ``introns=True``,
            or because the requested coordinates are before the start or
            after the end of the transcript.

        See Also
        --------
        transcript2genome : The inverse of this operation
        genome_interval2transcript : Convert intervals rather than single 
            positions

        Notes
        -----

        A key point to be aware of is that this function converts positions
        not intervals. Because of the zero-based half-open nature of python
        coordinates, you can't just convert the start and end positions
        into the transcript space. Use :method:`genome_interval2transcript`
        for this.

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
                                 (pos[i], self.transcript_id))
            
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
        '''Convert transcript coordinates into genome coordinates.

        Can be a single coordinate or a array-like of coordinates

        Parameters
        ----------
        pos : int or array-like or int
            positions, in the transcript domain, to be converted
        
        Returns
        -------
        int or numpy.array
            The position, or positions converted into the genome
            domain.

        Raises
        ------
        ValueError
            If the supplied genome position is not in the transcript.
            This would generatlly be because the supplied genome is
            either negative or greater than the length of the transcript.

        See Also
        --------
        genome2transcript : The inverse of this operation
        transcript_interval2genome_intervals : Convert intervals rather than
        single positions

        Notes
        -----

        A key point to be aware of is that this function converts positions
        not intervals. Because of the zero-based half-open nature of python
        coordinates, you can't just convert the start and end positions
        into the transcript space. Use 
        :method:`transcript_interval2genome_intervals` for this.

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
        '''Convert an interval in transcript coordinates and convert to
        a list of intervals in genome coordinates.

        Returns a list because a single interval might cross several
        exons, and the returned list of intervals would exclude the 
        intronic sequence.

        Parameters
        ----------
        interval : tuple of int
            Tuple with zero-based half open interval of form (start, end)
            in the transcript domain.

        Returns
        -------
        list of tuples
            List of zero-based, half open (start, end) tuples that encompas
            the same bases as described by `interval` but in the genome
            domain.

        Raises
        ------
        ValueError
            If the supplied strart or end is not in the transcript.

        See Also
        --------
        transcript2genome : Does the actaul conversion
        genome_interval2transcript : Almost the inverse of this

        Notes
        -----

        Because this method returns possibly discontinous intervals
        that will always only contain exonic regions (or only intronic if
        the converter was created with ``introns=True``), it is not quite
        guarenteed to be the inverse of genome_interval2transcript. This 
        is because a single interval that spans both introns and exons 
        can be passed to that method (as long as both start and end are
        exonic), and will be converted to a single interval that covers
        the exon bases, but converting the same interval back would give
        to intervals covering only the exonic parts.
        '''
        

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

    def genome_interval2transcript(self, interval):
        '''Convert an interval in genomic coordinates into an interval
        in transcript-coordinates.

        Parameters
        ----------
        interval : tuple of int
            Tuple with zero-based half open interval of form 
        (``start``, ``end``) in the genome domain. ``start`` < ``end``

        Returns
        -------
        tuple of int
            Half open (start, end) tuple starts and ends at the same bases
            as the interval described in `interval`, but in transcript 
            domain coordinates. 

        Raises
        ------
        ValueError
            If the supplied strart or end is not in the transcript.

        See Also
        --------
        genome2transcript : Does the actaul conversion
        genome_interval2transcript : Almost the inverse of this

        Notes
        -----

        This is not quite the inverse of 
        :method:`transcript_interval2genome_intervals`, because of how
        intervals that split across introns are handled. See
        :method:`transcript_interval2genome_intervals` for details of the
        difference.

        .. warning::
            This method current has a known issue where if the interval
            includes the last base of the transcripts (one the +ve strand)
            or the first base (on the -ve strand), an error will occur. I
            aim to fix this before release. 

        '''
        
        transcript_list = self.genome2transcript(interval)

        if self.strand == "+":
            return sorted(transcript_list)
        else:
            # half open
            transcript_list = transcript_list[0] + 1, transcript_list[1] + 1
            return sorted(transcript_list)
        
        
##################################################
def randomiseSites(profile, start, end, keep_dist=True):
    '''Randomise clipped sites within an interval.

    If keep_dist is true, then reads on the same base are kept togehter,
    preverving the distribution of per-base read counts.

    Parameters
    ----------
    profile : pandas.Series
        Cross-link profile to be randomised. Index is base, the value is
        the number of cross-links at that base.
    start : int
        Lowest bases at which allow randomised bases to fall.
    end : int
        Highest base at which to allow randomised bases to fall.
    keep_dist : bool, optional
        Maintain the count-per-base distribution by moving cross-links on
        the same base together. 

    Returns
    -------
    pandas.Series
        Randomised profile, in sparse format.

    
'''

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
    '''Extend cross-link sites in both directions.

    By default the cross-links are "spread" by the same amount in both
    directions. However, providing `right_bases` will make the spreading
    happen `bases` in the left direction and `right_bases` in the right
    direction. 

    Parameters
    ----------
    profile : pandas.Series
        Cross-link profile with positions to be spread
    bases : int
        Number of bases to spread the cross-link sites by. operates both
        left and right unless `right_bases` is provided.
    reindex : bool, optional
        Should the profile be reindexed before spreading. Only set
        ``False`` if the profile has already been reindexed (see below).
    right_bases : int, optional
        If set with extend `bases` in the left direction and `right_bases`
        in the right direction. 

    Returns
    -------
    pandas.Series
        The profile of the spread crosslinked bases. This will not be in 
        sparse form (will contain 0s).

    Notes
    -----

    The `reindex` parameter is provided as a covenience for when this has
    already been done. The following would need to be done if 
    ``reindex=False``:

        start = int(profile.index.min() - window)
        end = int(profile.index.max() + window+1)
        profile = profile.reindex(np.arange(start, end), fill_value=0)

    '''
    
    if right_bases:
        window = bases+right_bases
    else:
        window = 2 * bases + 1
        right_bases = bases

    if reindex:
        start = int(profile.index.min() - window)
        end = int(profile.index.max() + window+1)
        profile = profile.reindex(np.arange(start, end), fill_value=0)

    result = profile.rolling(window=window, center=False).sum().iloc[(window-1):]
    #    print result
    #    result =result.dropna()

    result.index = result.index - (right_bases + 1)
    return result


##################################################
def rand_apply(profile, exon, n, func, keep_dist=False,
               *args, **kwargs):
    '''Randomise a profile multiple times and apply a function
    to each randomised profile.
    
    Parameters
    ----------
    profile : pandas.Series
        profile to randomise, as produced by a *_getter function. 
    exon : CGAT.GTF.Entry or iCLIP.LiteExon
        Object with :attribute:`start` and :attribute:`end` specifying the
        boundaries to randomise between.
    n : int
        The number of randomisations to perform
    func : callable
        A function-like object to call on each randomised profile. The
        profile will be passed as the first argument to the function.
    keep_dist : bool, optional
        Maintain cross-links-per-base distribution when randomising
    *args, **kwargs : 
        Further arguments to `func`
        
    Returns
    -------
    pandas.Series or pandas.DataFrame
        Return type depends on the return type of `func`. Generally a 
        Series with the entries as the same return type as `func`. However,
        if `func` returns Series and all series are of the same length, 
        a pandas.DataFrame may be returned.
        
    See Also
    --------
    randomiseSites : Does the actual work of randomising.
      
    '''

    dummy = pd.Series(range(n))
 
    def _inner_func(x):
        rand = randomiseSites(profile, exon.start, exon.end,
                              keep_dist=keep_dist)
        return func(rand, *args, **kwargs)

    return dummy.apply(_inner_func)


