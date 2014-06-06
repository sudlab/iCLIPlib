'''
DedupUMI.py - Deduplicate reads that are coded with a UMI
=========================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The perpose of this script is to deduplicate BAM files based
on the first mapping co-ordinate and the UMI attached to the read.
It is assumed that the FASTQ files were processed with extract_UMI.py
and thus the UMI is the last word of the read name.

Two reads that have the same start co-ordinate and the same UMI and the same
splicing are assumed to be PCR duplicates.

The following criteria are applied to select the read that will be retained:

1. The read with the lowest nubmer of hits
2. The read with the highest mapping quality

Otherwise a read is chosen at random.

The input file must be sorted.

.. note::
   You need to redirect either logging information or output to a file or turn
   it off logging via -v 0 in order to get a valid sam/bam file.

Options
-------

The script assumes that input and output are in BAM format. Use that
-i/--in-sam and -o/--out-sam options to specifiy that the input/output
are SAM format.

The option --subset can be used to randomly select subsets of the reads
for futher processing.

Usage
-----

    python dedup_umi.py -I infile.bam -O deduped.bam

or

    cat infile.bam | python dedup_umi.py > deduped.bam 2> deduped.log

Command line options
--------------------

'''

import sys
import pysam
import CGAT.Experiment as E
import random
import collections

def dedup(insam, outsam, ignore_umi, subset):

    last_pos = 0
    last_chr = ""
    reads_dict = collections.defaultdict(dict)
    read_counts = collections.defaultdict(dict)
    
    for read in insam.fetch():

        if subset:
            if random.random() >= subset:
                continue

        if read.is_unmapped:
            continue

        if read.mate_is_unmapped:
            continue

        if read.is_read2:
            continue

        if read.is_reverse:
            pos = read.aend
        else:
            pos = read.pos
            
        if not read.pos == last_pos or not read.tid == last_chr:

            out_keys = [x for x in reads_dict.keys() if x < read.pos]
            
            for p in out_keys:
                for x in reads_dict[p].itervalues():
                    outsam.write(x)
                del reads_dict[p]
                del read_counts[p]


            
            last_pos = read.pos
            last_chr = read.tid

        if ignore_umi:
            umi = ""
        else:
            umi = read.qname.split("_")[-1]

        if 'N' in read.cigarstring:
            key = umi+"T"+str(read.is_reverse)+str(read.tlen)
        else:
            key = umi+"F"+str(read.is_reverse)+str(read.tlen)

        try:
            if reads_dict[pos][key].mapq > read.mapq:
                continue
        except KeyError:
            reads_dict[pos][key] = read
            read_counts[pos][key] = 0
        else:
            if reads_dict[pos][key].mapq < read.mapq:
                reads_dict[pos][key] = read
                read_counts[pos][key] = 0
                continue

            if reads_dict[pos][key].opt("NH") < read.opt("NH"):
                continue
            elif reads_dict[pos][key].opt("NH") > read.opt("NH"):
                reads_dict[pos][key] = read
                read_counts[pos][key] = 0

            read_counts[pos][key] += 1
            prob = 1.0/read_counts[pos][key]

            if random.random() < prob:
                reads_dict[pos][key] = read

     
def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--in-sam", dest="in_sam", action="store_true",
                      help="Input file is in sam format", default=False)
    parser.add_option("-o", "--out-sam", dest="out_sam", action="store_true",
                      help="Output alignments in sam format", default=False)
    parser.add_option("--ignore-umi", dest="ignore_umi", action="store_true",
                      help="Ignore UMI and dedup only on position",
                      default=False)
    parser.add_option("--subset", dest="subset", type="string",
                      help="Use only a fraction of reads, specified by subset",
                      default=1.1)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        in_name = "-"

    if options.stdout != sys.stdout:
        out_name = options.stdout.name
        options.stdout.close()
    else:
        out_name = "-"

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.out_sam:
        out_mode = "w"
    else:
        out_mode = "wb"

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode,
                            template=infile)

    dedup(infile, outfile, options.ignore_umi, subset=float(options.subset))
        
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
