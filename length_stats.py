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

        nreads+=1
        if read.inferred_length > read.alen:
            nlonger+=1
        if read.inferred_length >= options.length:
            out_hist[options.length] += 1
        else:
            out_hist[read.inferred_length] += 1

    header = ["Length", "Count"]
    outlines = ["\t".join(map(str, x)) for x in out_hist.iteritems()]
    outlines = "\t".join(header) + "\n" + "\n".join(outlines) + "\n"
    options.stdout.write(outlines)

    # write footer and output benchmark information.
    E.info("%i reads processed. %i reads have inferred reads length longer than aligned length" % (nreads, nlonger))
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
