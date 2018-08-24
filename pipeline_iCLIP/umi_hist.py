'''
umi_stats.py - Count frequency of UMIs
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Takes a BAM file produced from mapping a FASTQ created by
extract_umi.py and produces a histogram of their usage.

Options
-------

Will read bamfile off of the stdin or from -I, but 
BAM file must be indexed and so cannot be manipulated on 
stdin.


Usage
-----


python umi_stats.py -I [BAMFILE]



Command line options
--------------------

'''


import sys
import collections
import CGAT.Experiment as E
import pysam


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    out_hist = collections.defaultdict(int)

    if options.stdin == sys.stdin:
        in_bam = pysam.Samfile("-", "rb")
    else:
        fn = options.stdin.name
        options.stdin.close()
        in_bam = pysam.Samfile(fn, "rb")

    for read in in_bam.fetch():
        barcode = read.qname.split("_")[-1]
        out_hist[barcode] += 1

    header = "\t".join(["UMI", "Count"])
    outlines = ["\t".join(map(str, x)) for x in out_hist.iteritems()]
    outlines = header + "\n" + "\n".join(outlines) + "\n"
    options.stdout.write(outlines)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
