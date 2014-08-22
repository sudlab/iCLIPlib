'''
find_significant_bases.py - template for CGAT scripts
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script calculates p values for cross linked positions in transcripts from
an iCLIP experiment. Input (from stdin) is a gtf file containing transcripts
and (specified to the script) a bamfile containing the mapped reads. The output
is a bedgraph file containing the (optionally corrected) pvalues for
significant cross-link clustering around each crosslinked base.

The model is based on several assumptions,

1)  The null probability of a cross-link at any position in a transcript is
    equal. Thus:

2)  The probability of getting x_i read counts at any given position is
    distributed binomially. :math:`X_i ~ Bin(n,1/p)` where n is the total
    number of counts mapping to a region and p is the number of bases in the
    region.

3)  The probability of getting x_win read counts within a window around a
    cross-linked base is distributed :math:`X_win ~ Bin(n-x_i, 1/p_win)` where
    p_win is p minus the window size.

4)  Thus the probability of getting x_i reads or more at position i and X_win
    or more reads around it is:
    
    :math:`P(X_i>=x_i)*P(X_win>=x_win)`

5)  All transcripts from a gene are equally expressed. Where a crosslink
    position is with in more than one transcript the probability for that base
    is. Thus for 2 transcripts:
    
    :math:`P(X>x) = sum (P(X>=x| transcript i)*P(transcript i))`
           for i in transcripts

    as all transcripts are equally expressed this becomes

    :math:`P(X>x) = mean(P(X>=x | transcript i))`
           for i in transcripts.

    i.e. overlapping bases are averaged. This is probably not good. If you
    don't want this please deal with overlapping transcripts before running
    the script.

.. Note:: Input GTF should be sorted such that all lines from the same
         transcript are consecutive, and all transcripts from the same gene are
         consecutive.

.. Warning:: Genes must not overlap. If they do the overlapping bases will
             appear more than once in the output.

.. TODO:: Allow normalisation by rnaseq input
Options
-------

-g, --groupby: How to divide the transcript into regions that will be tested
               together. The default is exons. Here all exons from a transcript
               are concatenated, and treated as a single region, and the
               introns likwise.

               Alternatively, the utrs can also be treated as seperate regions
               (--groupby=utrs). If this is the case the GTF must contain CDS
               entires. Or the whole genomic region from tss to tts can be
               treated as a single region.

-w, --window-size: The script uses windows around the cross-linked base to
               determine clustering. This parameter specifies the distance
               either side of the cross-linked base to use.

-p, --pipeout: Specifies that results should be written to the output as soon
               as generated. This will save memory, and allow downstream
               computations, but may produce a truncated output if there is an
               error and does not allow an FDR correction.

-f, --fdr:     Compute an BH FDR correction on the results.
               Implies not --pipeout.

-t, --dtype:   The numpy dtype to use for storing counts. The default is
               uint32. Smaller types will use less memory, but run the risk of
               integer overflow (detected).

-
Usage
-----

.. Example use case

Example::

   python find_significant_bases.py mybam.bam > outfile.bg < transcripts.gtf

Usage::

   python find_significant_bases.py [OPTIONS] BAMFILE < GTFFILE


Command line options
--------------------

'''

from scipy.stats import binom
import pandas as pd
import sys
import CGAT.Experiment as E
import iCLIP
import CGAT.GTF as GTF
import pysam
from statsmodels.stats.multitest import multipletests


class DeferredOutput:
    ''' This class looks like a file like object, but stores up all
    the objects passed for output, and outputs them all when close
    is called. Optionally performs multiple testing correction '''

    def __init__(self, outfile, correct):

        self.outfile = outfile
        self.output = pd.Series(index=pd.MultiIndex(levels=[[], []],
                                                    labels=[[], []],
                                                    names=["contig",
                                                           "position"]),
                                dtype="float_")
        self.correct = correct

    def write(self, gene_results):

        self.output = self.output.append(gene_results)

    def close(self):

        if self.correct:
            E.info("Correcting p-values using BH ...")
            corrected_pvals = multipletests(self.output, method="fdr_bh")
            self.output = pd.Series(corrected_pvals[1], index=self.output.index)

        E.info("Writing output")
        E.debug("output contains %i entries" % len(self.output))
        self.output.to_csv(self.outfile,
                           sep="\t",
                           header=False)


class InstantOutput:
    ''' This class looks file a file like object but takes pandas Series
    objects and outputs them to a file handle it keeps open '''

    def __init__(self, outfile):
        self.outfile = outfile

    def write(self, gene_results):
        gene_results.to_csv(self.outfile,
                            sep="\t",
                            header=False)

    def close(self):
        pass


def calculateProbabilities(counts, window_size, length, start=0):
    '''Calculates the probablity of observing the counted
    number of reads in windows of "window_size" around each
    cross-linked size.

    Currently assumes that coordinates are in transcript space.

    The length of the transcript must be provided because it
    is needed for calculating the prior and may be outside of the
    provided coordinates.

    Start allows, together with length, for only using part of
    counts.'''

    E.debug("Calculate probabilies called with:")
    E.debug((counts,window_size, length, start))
    # limit to subset
    counts = counts[(counts.index.values >= start) &
                    (counts.index.values < (start + length))]
    total_counts = counts.sum()
    
   
    single_base_ps = 1 - binom.cdf(counts-2, total_counts, 1.0/length)

    window_ps = pd.Series(dtype="float_")
    for base in counts.index.values:

        window_start = max(0, base-window_size)
        window_end = min(start+length, base + window_size)
        try:
            window = counts.index.slice_indexer(start=float(window_start),
                                                end=float(window_end),
                                                step=1)
        except KeyError:
            print (window_start,window_end)
            print counts
  
        height = counts.iloc[window].sum() - counts[base]

        
        p = (window_end - window_start) / length
        E.debug((window_start,window_end))
        window_ps[base] = (1 - binom.cdf(height - 1, total_counts, p))
        E.debug("Heigh %i, total_counts %i, p %f" % (height, total_counts, p))
        E.debug("window p is %f" % window_ps[base])

    return single_base_ps * window_ps


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    grouping_choices = ["exons",
                        "utrs",
                        "all"]
    parser.add_option("-g", "--grouping", dest="grouping", type="choice",
                      choices=grouping_choices,
                      help="How to group transcript regions choices are [%s]"
                            % ",".join(grouping_choices))
    parser.add_option("-p", "--pipeout", dest="pipeout", action="store_true",
                      help="Output continuously to the pipe rather than in a"
                           "chunk at the end")
    parser.add_option("-t", "--dtype", dest="dtype", type="string",
                      default="uint32",
                      help="Numpy dtype for storing counts")
    parser.add_option("-w", "--window-size", dest="window_size",
                      type="int", default=15,
                      help="Size of window either size of crosslinked base to"
                           "consider")
    parser.add_option("-f", "--fdr", dest="fdr", action="store_true",
                      default=False,
                      help="perform BH fdr correction on p-values, implies not"
                           "--pipeout")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # Standard in contains the transcripts
    
    gffs = GTF.gene_iterator(GTF.iterator(options.stdin))

    # bam file is the first positional arguement
    bamfile = pysam.Samfile(args[0])

    if options.fdr and options.pipeout:
        E.warning("--fdr implies not --pipeout, instant output disabled")
        options.pipeout = False

    if options.pipeout:
        output = InstantOutput(options.stdout)
    else:
        output = DeferredOutput(options.stdout, options.fdr)

    E.info("Counting accross transcripts ...")

    for gene in gffs:

        if options.grouping == "all":
            gene = GTF.merged_gene_iterator(gene)

        transcript_ps = {}

        for transcript in gene:
            coords_converter = iCLIP.TranscriptCoordInterconverter(transcript)
            exons = GTF.asRanges(transcript, "exon")
            counts = iCLIP.count_intervals(bamfile,
                                           exons,
                                           strand=transcript[0].strand,
                                           contig=transcript[0].contig,
                                           dtype=options.dtype)

            counts.index = coords_converter.genome2transcript(counts.index.values)
            counts = counts.sort_index()
            cds = GTF.asRanges(transcript, "CDS")

            if options.grouping == "utrs" and len(cds) > 0:
                
                cds_interval = (cds[0][0], cds[-1][1])
                cds_interval = coords_converter.genome2transcript(cds_interval)
                cds_interval.sort()
                cds_length = cds_interval[1] - cds_interval[0]

                p_intervals = [(0, cds_interval[0]),
                               (cds_interval[0], cds_length),
                               (cds_interval[1], coords_converter.length - cds_interval[1])]

            else:  # do not group by cds or there is no cds
                p_intervals = [(0, coords_converter.length)]

            p_values = [calculateProbabilities(counts, options.window_size,
                                              length=length, start=start)
                        for start, length in p_intervals
                        if length > 0]

            if len(p_values) > 1:
                p_values = pd.concat(p_values)
            else:
                p_values = p_values[0]

            p_values.index = coords_converter.transcript2genome(p_values.index.values)

            intron_intervals = GTF.toIntronIntervals(transcript)

            if len(intron_intervals) > 0:
                intron_coords = iCLIP.TranscriptCoordInterconverter(transcript,
                                                                    introns=True)
                intron_counts = iCLIP.count_intervals(bamfile,
                                                      intron_intervals,
                                                      strand=transcript[0].strand,
                                                      contig=transcript[0].contig,
                                                      dtype=options.dtype)
             
                intron_counts.index = intron_coords.genome2transcript(
                    intron_counts.index.values)
                intron_counts = intron_counts.sort_index()
                intron_pvalues = calculateProbabilities(intron_counts,
                                                        options.window_size,
                                                        intron_coords.length)
                                                        
                intron_pvalues.index = intron_coords.transcript2genome(
                    intron_pvalues.index.values)
                p_values = p_values.append(intron_pvalues)
                
            transcript_ps[transcript[0].transcript_id] = p_values

        transcript_df = pd.DataFrame(transcript_ps)

        # convert to 1-based coords
        transcript_df.index = transcript_df.index + 1

        transcript_df.index.rename("position", inplace=True)
        transcript_df["contig"] = gene[0][0].contig
        transcript_df.set_index("contig", append=True, inplace=True)
        transcript_df = transcript_df.reorder_levels(["contig", "position"])

        gene_ps = transcript_df.mean(1)
    
        output.write(gene_ps)

    output.close()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
