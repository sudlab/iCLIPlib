import re
import pandas as pd 
import numpy as np
import collections
import itertools

import CGAT.GTF as GTF
import CGAT.Experiment as E

from utils import spread, rand_apply, AMBIGUITY_CODES
from counting import count_transcript, count_intervals


##################################################
def find_all_matches(sequence, regexes):
	'''Return start positions of a collection of regexes

	Parameters
	----------
	sequence : str
		String to search in. Usually an RNA or DNA sequence
	regexs : itr of re.RegEx
		List of regexs to search :param:`sequence` for.
		
	Returns
	-------
	pandas.Series of np.array of int
		Series of arrays where each array is the integer 
		positions matches for one regex. Index is the regex objects.
		
	'''
		
    def _find_regex(regex):
        matches = regex.finditer(sequence)
        mStarts = [m.start() for m in matches]
        return np.asarray(mStarts)

    return regexes.apply(_find_regex)


##################################################
def pentamer_frequency(profile, length, regex_matches, nSpread=15):
    '''Calculate the frequency of overlaps between a list of cross-link
	sites (possibly extended) and a collection of sites
	
	The second collection of site would usually be a Series of ararys of
	positions as returned by :func:`find_all_matches`.

	Parameters
	----------
	profile : pandas.Series
		A profile of the number of cross-links at each base.
	length : int
		Length of the sequence represetnted by `profile`.
	regex_matches : pandas.Series of nd.array of int
		Each array entry represents a single match on the sequence. Each 
		array represents a different thing matched (e.g. a different regex)
		. This structure is would usually be returned by
		:func:`find_all_matches`.
	nSpread : int, optional
		How far either side of a cross-link location to consider when 
		calculating overlaps (defaults to 15).
	
	Returns
	-------
	pandas.Series of int
		Each entry is the number of overlaps between cross-links sites
		in profile and a single entry in `regex_matches`. Each same index
		as `regex_matches`.
		
	'''
	
    kmer = len(regex_matches.index.values[0])
    profile = profile.reindex(
        np.arange(-nSpread, length+nSpread-kmer), fill_value=0)
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
    ''' This function calculates the z-score of enrichment of kmers around
    cross-linked sites.
	
	Parameters
	----------
	gtf_chunk_iterator : iter of list-like of `CGAT.GTF.Entry`
		An iterator that returns a list-like object containing 
		`CGAT.GTF.Entry` objects representing a gene or transcript. 
	bam : *_getter func
		A getter function, as returned by :func:`make_getter`
	fasta : CGAT.IndexedFasta
		Indexed fasta containing the genome sequence to get the transcript
		sequences from. 
	kmer_length : int, optional
		length of kmer to search for enrichment of. All possible kmers of 
		this length will be tested. Defaults to 5.
	randomisations : int, optional
		Number of sequence randomisation to perform in order to determined
		mean and sd of the frequency of each kmer under a the null model.
	seperate_UTRs : bool, optional
		No effect **to be removed**
	spread : int, optional
		Number of bases to consider around a clip site when calculating an
		overlap. Defaults to 15.
	pool : multiprocessing.Pool, optional
		If present, work will be parallelized across the worker pool
		
	Returns
	-------
	pandas.Series of float
		z-values of the frequency of overlaps between cross-link sites and
		all kmers. Series is the kmer in question
		
	Notes
	-----
	
	For each transcript/gene/collection of intervals in
	`gtf_chunk_iterator`, the count of cross-link sites is calculated. Each
	cross-link site is then extended by `spread` bases in both directions. 
	
	For every possible kmer, the start positions in the sequence of the
	transcript is then calculated and the number of times a cross-link
	site overlaps a kmer start position is counted. 
	
	The positions of the cross-link site are then randomized. And the
	above process is repeated on each of the randomized profiles. 
	
	This is performed on introns and exons seperately. 
	
	Once this has been performed on all transcripts/genes in the iterator,
	scores are summed across all transcripts/genes and the mean and
	standard deviation of the counts for each kmer is calculated across
	the randomisations. For each kmer :math:`i \in 0,1,...4^kmer` the
	z-score :math:`z_i` is calculated:
	
	.. math::
		z_i = \frac{observed_counts - \mu_i}{\sigma_i}
		
	where :math:`\mu_i` is the mean frequency of overlaps between
	cross-links and the i*th* kmer across the randomisations and
	:math:`\sigma_i` is the corresponding standard deviation. 
	
	This process can be very time consuming for longer kmers. To help 
	elleviate this, it can be parrellized accross cores of the machine by
	providing a `multiprocessing.Pool` object to the `pool` parameter. The
	task will then be parrellized, gene-wise. 
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
