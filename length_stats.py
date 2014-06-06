'''
length_stats.py - Count frequency of UMIs
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Take a bam file and calculate alignment length histogram.

Options
-------

-l --length
     The read length. Any aligned length longer than this is
     truncated to this length. 


Usage
-----


python length_stats.py -I [BAMFILE]



Command line options
--------------------

'''

import os
import sys
import re
import optparse
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

    parser.add_option("-l","--length", dest = "length", type="int",
                      help="max read length", default = 100)
    parser.add_option("-p", "--paired", dest="paired", type = "int",
                     help="Data is paired. Use fragment length where aligned"
                          "length is greater than or equal to this length", default=None)
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    out_hist = collections.defaultdict(int)
    
    if options.stdin == sys.stdin:
        in_bam = pysam.Samfile("-", "rb")
    else:
        fn = options.stdin.name
        options.stdin.close()
        in_bam = pysam.Samfile(fn, "rb")

    nreads=0
    nlonger=0
    for read in in_bam.fetch():

        if read.is_read2 or read.is_unmapped or read.mate_is_unmapped:
            continue

        nreads+=1
        length = read.inferred_length - sum([l for o,l in read.cigar if o=="S"])
        if read.inferred_length > length:
            nlonger+=1

        if not (options.paired is None) and length >= options.paired:

            splice_1 = max(0, read.alen-length)
            
            out_hist[abs(read.tlen) - splice_1] += 1
        else:
            out_hist[length] += 1

    header = ["Length", "Count"]

    keys = sorted(out_hist.keys())

    outlines = ["\t".join(map(str, (key,out_hist[key]))) for key in keys]
    outlines = "\t".join(header) + "\n" + "\n".join(outlines) + "\n"
    options.stdout.write(outlines)

    # write footer and output benchmark information.
    E.info("%i reads processed. %i reads have inferred reads length longer than aligned length" % (nreads, nlonger))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
