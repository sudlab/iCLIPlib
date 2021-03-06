'''
cgat_script_template.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Count the number of clip sites in each gene, transcript or exon. GTF file is supplied to stdin.
Bamfile is supplied as the first positional arguement.

Centres can be counted with --use-centres


Usage
-----

python count_clip_sites.py BAMFILE [OPTIONS]


.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Intervals as Intervals
import os
import pysam

sys.path.insert(1, os.path.join(
    os.path.dirname(__file__)
    ,".."))

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

    parser.add_option("-f", "--feature", dest="feature", type="choice",
                      choices=["gene", "transcript", "exon"],
                      default="transcript",
                      help="which feature to use: gene/transcript/exon")
    parser.add_option("--unstranded-bw", dest="unstranded_wig", type="string",
                      help="BigWig with tag counts on both strands")
    parser.add_option("--plus-bw", dest="plus_wig", type="string",
                      help="BigWig with tag counts from plus strand")
    parser.add_option("--minus-bw", dest="minus_wig", type="string",
                      help="BigWig with tag counts from minus strand")
    parser.add_option("--bed", dest="bedfile", type="string",
                      help="tabix indexed bed file with tag counts"),
    parser.add_option("-c", "--use-centre", dest="centre", action="store_true",
                      default=False,
                      help="Use centre of read rather than start")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    iterator = GTF.iterator(options.stdin)

    if options.feature == "gene":
        iterator = GTF.flat_gene_iterator(iterator)
    elif options.feature == "transcript":
        iterator = GTF.transcript_iterator(iterator)
    elif options.feature == "exon":
        def _exon_iterator(gff_iterator):
            for exon in gff_iterator:
                yield [exon]
        iterator = _exon_iterator(iterator)

    if options.unstranded_wig:
        bamfile = iCLIP.make_getter(plus_wig=options.unstranded_wig)
    elif options.plus_wig:
        if not options.minus_wig:
            raise ValueError(
                "Please provide wigs for both strands or use --unstranded_wig")
        bamfile = iCLIP.make_getter(plus_wig=options.plus_wig,
                                    minus_wig=options.minus_wig)
    elif options.bedfile:
        bamfile = iCLIP.make_getter(bedfile=options.bedfile)   
    else:
        bamfile = pysam.AlignmentFile(args[0])
        
    outlines = []
    for feature in iterator:
        exons = GTF.asRanges(feature, "exon")

        exon_counts = iCLIP.count_intervals(bamfile,
                                            exons,
                                            feature[0].contig,
                                            feature[0].strand,
                                            dtype="uint32",
                                            use_centre=options.centre)

        exon_counts = exon_counts.sum()

        introns = Intervals.complement(exons)
        intron_counts = iCLIP.count_intervals(bamfile,
                                              introns,
                                              feature[0].contig,
                                              feature[0].strand,
                                              dtype="uint32",
                                              use_centre=options.centre)

        intron_counts = intron_counts.sum()

        if options.feature == "exon":

            try:
                exon_id = feature[0].exon_id
            except AttributeError:
                try:
                    exon_id = feature[0].exon_number
                except AttributeError:
                    exon_id = "missing"

            gene_id = feature[0].gene_id
            transcript_id = feature[0].transcript_id
            intron_counts = "NA"
        else:
            exon_id = "NA"
            gene_id = feature[0].gene_id
            transcript_id = feature[0].transcript_id
            intron_counts = float(intron_counts)
            
        outlines.append([gene_id,
                         transcript_id,
                         exon_id,
                         str(float(exon_counts)),
                         str(intron_counts)])

    options.stdout.write("\t".join(["gene_id",
                                    "transcript_id",
                                    "exon_id",
                                    "exon_count",
                                    "intron_count"])+"\n")

    outlines = ["\t".join(outline) for outline in outlines]
    outlines = "\n".join(outlines)
    options.stdout.write(outlines + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
