'''
iCLIP_bam2geneprofile.py - produce geneprofile of iCLIP sites
============================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

While bam2geneprofile in CGAT is a very flexible tool, it is not
neccesarily suitable for iCLIP as we should only consider the first
base (or any mutant bases).

This script wraps iCLIP.meta_gene to produce metagene profiles of
iCLIP bam files. It can also use bigwig files - provide either unstranded
data to `--plus-wig` or stranded data by also using `--minus-wig`.

Profiles are always normalised to the sum of each window before summing 
accross windoes, but the full matrix may also be output if custom normalisation
is required.

Using single bases means we don't have to worry about over sampling single
reads, so profiles should be less resolution sensitive

Usage
-----

.. Example use case

Example::

   python iCLIP_bam2geneprofile.py -I geneset.gtf.gz mybam.bam

Type::

   python iCLIP_bam2geneprofile.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import CGAT.Experiment as E
import pysam
import CGAT.IOTools as IOTools

sys.path.insert(1, os.path.join(
    os.path.dirname(__file__), ".."))

import iCLIP


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--output-matrix", dest="matrix", type="string",
                      default=None,
                      help="output full matrix to this file")
    parser.add_option("-f", "--flanks", dest="flanks", type="int",
                      default=100,
                      help="number of basepairs to use for gene flanks")
    parser.add_option("-b", "--exon-bins", dest="exon_bins", type="int",
                      default=1000,
                      help="number of bins to divide transcripts into")
    parser.add_option("--flank-bins", dest="flank_bins", type="int",
                      default=10,
                      help="number of bins to divide flanks into")
    parser.add_option("--scale-flanks", dest="scale_flanks", action="store_true",
                      default=False,
                      help="Scale the size of the flank bins to match the size of the"
                      "exon bins for each transcript")
    parser.add_option("--pseudo_count", dest="pseudo_count", type="float",
                      default=0,
                      help="add pseduo count to bins to mitiage effects of low numbers of reads")
    parser.add_option("--normalised_profile", dest="normalize_profile", action="store_true",
                      default=False,
                      help="Normlize profile by profile sum")
    parser.add_option("--plus-wig", dest="plus_wig",
                      default=None,
                      help="Use this wig file instead of a BAM file to get clip density"
                      "may be used as only wig file, or may be provided together with"
                      "--minus-wig for standed computation")
    parser.add_option("--minus_wig", dest="minus_wig", default=None,
                      help="Use this to provide stranded wig data")
    parser.add_option("--bed", dest="bedfile", default=None,
                      help="Use bed file with signal instead of bam")
    parser.add_option("--centre", dest="centre", action="store_true",
                      default=False,
                      help="Use centre of read rather than end")
    parser.add_option("--no-gene-norm", dest="row_norm", action="store_false",
                      default=True,
                      help="Do not normalise profile from each gene")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.plus_wig:
        bam = iCLIP.make_getter(plus_wig=options.plus_wig,
                                minus_wig=options.minus_wig)
    elif options.bedfile:
        bam = iCLIP.make_getter(bedfile=options.bedfile)
    else:
        bam = iCLIP.make_getter(bamfile=args[0], centre=options.centre)

    if options.flanks > 0:
        bins = [options.flank_bins,
                options.exon_bins,
                options.flank_bins]
    else:
        bins = options.exon_bins

        
    summed_matrix, counts_matrix = iCLIP.meta_gene(
        options.stdin,
        bam,
        bins,
        options.flanks, 
        output_matrix=(options.matrix is not None),
        calculate_flanks=options.scale_flanks,
        pseudo_count=options.pseudo_count,
        row_norm=options.row_norm)

    if options.flanks > 0:
        summed_matrix = summed_matrix[["flank5", "exons", "flank3"]]

    summed_matrix = summed_matrix.reset_index()
    if options.normalize_profile:
        summed_matrix["density"] = summed_matrix["density"]/summed_matrix["density"].sum()

    summed_matrix.to_csv(options.stdout, sep="\t",
                         index=True,
                         index_label="bin")

    if options.matrix:
        counts_matrix = counts_matrix.transpose()

        counts_matrix = counts_matrix.loc[["flank5", "exons", "flank3"],:]
        counts_matrix = counts_matrix.reset_index(drop=True)
        counts_matrix = counts_matrix.transpose()

        counts_matrix.to_csv(IOTools.openFile(options.matrix, "w"),
                             sep="\t",
                             index=True,
                             index_label="transcript_id")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
