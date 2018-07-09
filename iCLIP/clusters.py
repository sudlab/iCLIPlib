'''This file contains functions and classes relating to calling iCLIP
clusters'''

import CGAT.Experiment as E
import numpy as np
import pandas as pd
import CGAT.GTF as GTF

from utils import spread, rand_apply, TranscriptCoordInterconverter
from counting import count_transcript, count_intervals
from kmers import LiteExon


def Ph(profile, exon, nspread):
    '''Calculates a Series, Ph, such that Ph[i] is the
    P(X >= i) where X is the height of signal on a base of the
    profile'''

    profile = profile.reindex(np.arange(exon.start-nspread,
                                        exon.end + nspread + 1),
                              fill_value=0)
    profile = spread(profile, nspread, reindex=False)
    profile = profile[profile > 0]
    pdf = profile.value_counts()
    pdf = pdf/pdf.sum()
    pdf = pdf.reindex(np.arange(0, pdf.index.values.max()), fill_value=0)
    cdf = 1 - pdf.cumsum()
    cdf.index = cdf.index + 1
    return cdf


def fdr(profile, exon, nspread, randomizations):
    '''Calculate the FDR of finding a particular heights
    by using randomizations'''

    profile_Ph = Ph(profile, exon, nspread)
    rands = rand_apply(profile, exon,
                       randomizations,
                       Ph,
                       False,
                       exon,
                       nspread)
    rands = rands.reindex(columns=profile_Ph.index)
    rands = rands.fillna(0)
    muh = rands.mean()
    sigmah = rands.std()
    fdr_thresholds = (muh + sigmah) / profile_Ph
    spread_profile = spread(profile, nspread)
    fdrs = spread_profile.map(fdr_thresholds)
    try:
        fdrs = fdrs.loc[profile.index]
    except KeyError:
        print profile.index
        print spread_profile.index
        print fdrs.index
        raise
    fdrs = fdrs.reindex(profile.index)
    return fdrs


def _get_profiles_and_conveter(gtf_iterator, bam):

    for transcript in gtf_iterator:

        transcript = [e for e in transcript if hasattr(e, "transcript_id")]
        if len(transcript) == 0:
            continue
        
        gene_id = transcript[0].gene_id
        #transcript_id = transcript[0].transcript_id
        contig = transcript[0].contig
        strand = transcript[0].strand
         
        E.debug("Crunching gene: %s:"
                % gene_id)
         
        # exons
        profile = count_transcript(transcript, bam)

        if profile.sum() > 0:
            converter = TranscriptCoordInterconverter(transcript)
            yield (profile, converter, LiteExon(0, converter.length),
                   contig, strand)

        # introns
       
        intron_intervals = GTF.toIntronIntervals(transcript)
        intron_counts = count_intervals(bam, intron_intervals,
                                        contig, strand)
        if intron_counts.sum() == 0:
            continue

        converter = TranscriptCoordInterconverter(transcript,
                                                  introns=True)
        intron_counts.index = converter.genome2transcript(
            intron_counts.index.values)

        for intron in intron_intervals:
            intron = (intron[0], intron[1] - 1)
            intron = converter.genome2transcript(intron)
            intron = sorted(intron)
            intron = (intron[0], intron[1] + 1)
            profile = intron_counts.loc[float(intron[0]):float(intron[1])]
            if profile.sum() > 0:
                yield (profile, converter, LiteExon(*intron), contig, strand)


def _get_fdr_for_transcript(profile, exon, nspread, randomizations,
                            converter, contig, strand):

    fdrs = fdr(profile, exon, nspread, randomizations)
    fdrs.name = "fdr"
    profile.name = "Depth"
    fdrs = pd.concat([fdrs, profile], axis=1)
    fdrs.fdr = fdrs.fdr.fillna(1)
    fdrs.Depth = fdrs.Depth.fillna(0)
    fdrs["strand"] = strand
    fdrs.index = pd.MultiIndex.from_tuples(
        [(contig, x) for x in converter.transcript2genome(fdrs.index.values)])
    return fdrs


def _par_get_fdr_for_transcript(args):
    '''Multiprocessing wrapper for _get_fdr_for_transcript'''

    return _get_fdr_for_transcript(*args)
    

def get_crosslink_fdr_by_randomisation(gtf_iterator, bam,
                                       randomisations=100,
                                       nspread=15, pool=None):
    ''' This function will carry out the assessment of crosslink site
    significance using the method outlined in Wang Z et al.

    Breifly the empirical distribution of hieights is calculated
    accross the transcript after calculating the height of any one
    base by summing reads sites within 15 bases.

    FDR is then assessed by comparison to randomized profiles.

    Results will be returned in the order of the gtf_iterator unless
    paralellization is used, in which case not particular order is
    guarenteed.

        :param gtf_iterator: An iterator that returns listss of
                             CGAT.GTF.Entry
        :type bam: pysam.AlignmentFile
        :param pool: If a worker pool is provided work will be
                     parallelised accross the pool

        :rtype: pd.Series with a MultiIndex first level contig,
                second level base'''

    args = ((profile, exon, nspread, randomisations, converter, contig, strand)
            for profile, converter, exon, contig, strand 
            in _get_profiles_and_conveter(gtf_iterator, bam))
    if pool:
        results = pool.imap(_par_get_fdr_for_transcript, args)
    else:
        results = (_get_fdr_for_transcript(*arg) for arg in args)

    results = list(results)
    
    results = pd.concat(results)
        
    return results

    
