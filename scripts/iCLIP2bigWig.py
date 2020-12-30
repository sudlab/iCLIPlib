'''
iCLIP2bigWig -- convert iCLIP BAM files to two wig files
============================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Takes a bam file on the stdin and converts it into two
bigWig files, each containing the depth of crosslinked bases
at any given position on each strand.

Output files are named according to the provided template
with _plus and _minus suffixes.

Requires the ucsc wigToBigWig program to work properly.

Options
-------

Will read bamfile off of the stdin or from -I, but 
BAM file must be indexed and so cannot be manipulated on 
stdin.

If wig files are required as output, for example if 
wigToBigWig is not installed, --wig will output wig files
rather than bigWig files.


Usage
-----


python umi_stats.py -I [BAMFILE] [OUT_TEMPLATE]



Command line options
--------------------

'''

import sys
import os
import shutil
import CGAT.Experiment as E
import pysam
import subprocess
import tempfile
import CGAT.IOTools as IOTools
import pandas as pd

sys.path.insert(1, os.path.join(
    os.path.dirname(__file__), ".."))

import iCLIP


def outputToBW(infile, outfile_prefix, chrom_sizes):

    E.debug("Attempting to output %s to bigWig" % outfile_prefix)
    command = ["bedGraphToBigWig",
               infile,
               chrom_sizes,
               outfile_prefix + ".bw"]
    try:
        subprocess.check_call(command)
    except Exception as e:
        E.error("Error on conversion: %s" % e)
        E.info("Outputting to wig")
        shutil.move(infile, outfile_prefix + ".bg")
    else:
        E.debug("Conversion successful")
        os.unlink(infile)
    
def outputToBG(depths,chrom, chrom_size, bgfile):
    '''depths is a pandas series keyed on chromosome position,
    chrom is a chromosome, wigfile is a file to output to.
    This function converts a series of depths into wig formated
    text and writes it to the specified file '''

    depths.index = depths.index
    for row in depths.items():
        row = list(row)
        if row[0] <= 0 or row[0] >= chrom_size-1:
            continue
        
        row = "\t".join(map(str,[chrom,
                                 int(row[0]),
                                 int(row[0])+1,
                                 row[1]])) + "\n"
        bgfile.write(row)

def output2Bed(pos_depth, neg_depth, chrom, outfile):

    pos_depth.name = "+"
    neg_depth.name = "-"
    neg_depth = abs(neg_depth)
    
    both_depth = pd.concat([pos_depth, neg_depth], axis=1)

    def _row_to_bed(row, strand):
        return "\t".join(
            map(str, [chrom,
                      int(row[0]),
                      int(row[0])+1,
                      ".",
                      row[1][strand],
                      strand]))

    for row in both_depth.iterrows():
        if row[1]["+"] > 0:
            outfile.write(_row_to_bed(row, "+") + "\n")

        if row[1]["-"] > 0:
            outfile.write(_row_to_bed(row, "-") + "\n")
                          
    
def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    profiles = list(iCLIP.getters.profiles.keys())
    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-p", "--profile", dest="profile", type="choice",
                      choices=profiles,
                      default="iclip",
                      help="Experiment profile to use. Sets various things"
                      "about obtaining 1-bp position from read. Options are"
                      " %s" % ", ".join(profiles))
    parser.add_option("-c", "--use-centre", dest="centre", action="store_true",
                      default=None,
                      help="Use centre of read rather than frist base."
                      "Overrides profile")
    parser.add_option("-f", "--format", dest="format",
                      choices=["bigWig", "bigwig", "BigWig",
                               "bedGraph", "bg", "bedgraph",
                               "bed", "Bed", "BED"],
                      help="Output format. Either bigWig (2 files, + and - strand)"
                           ", bedGraph (2 files), or bed (1 file, depth in column 5,"
                           "strand in column 6",
                      default="bigWig")
    parser.add_option("-w", "--wig", dest="output_wig", action="store_true",
                      default=False,
                      help="Write output to bedgraph file rather than bigwig")
    parser.add_option("--dtype", dest = "dtype", type="string",
                      default="uint32",
                      help="dtype for storing depths")
    parser.add_option("--cpm", dest="cpm", action="store_true",
                      default=False,
                      help="Normalize output depths to number of mapped reads (in millions) in BAM")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    options.format = options.format.lower()
    if options.format == "bg":
        options.format = "bedgraph"

    profile = iCLIP.getters.profiles[options.profile]
    
    if options.centre is not None:
        centre=True
    else:
        centre=profile.centre
    
    if options.stdin == sys.stdin:
        in_bam = pysam.Samfile("-", "rb")
                                 
    else:
        fn = options.stdin.name
        options.stdin.close()
        in_bam = pysam.Samfile(fn, "rb")
                                  

    getter = iCLIP.make_getter(in_bam, profile=profile, centre=centre)

    if options.cpm:
        scale_factor = sum(contig.mapped for contig in in_bam.get_index_statistics())

        scale_factor = 1000000.0/scale_factor

    if options.format == "bed":
        bedfile = IOTools.openFile(args[0], "w")
    else:
        plus_wig = tempfile.NamedTemporaryFile(delete=False)
        minus_wig = tempfile.NamedTemporaryFile(delete=False)

    contig_sizes = []

 
    for chrom, chrom_length in zip(in_bam.references, in_bam.lengths):

        # get depths over chromosome
        pos_depth, neg_depth, counter = getter(chrom, strand="both", dtype=options.dtype)
        pos_depth_sorted = pos_depth.sort_index()
        del pos_depth
        neg_depth_sorted = neg_depth.sort_index()
        del neg_depth
        neg_depth_sorted = -1*neg_depth_sorted
 
        if options.cpm:
            pos_depth_sorted = pos_depth_sorted * scale_factor
            neg_depth_sorted = neg_depth_sorted * scale_factor

        if options.cpm:
            pos_depth = pos_depth * scale_factor
            neg_depth = neg_depth * scale_factor
            
        # output to temporary wig file
        if options.format == "bed":
            output2Bed(pos_depth_sorted, neg_depth_sorted, chrom, bedfile)
        else:
            outputToBG(pos_depth_sorted, chrom, chrom_length, plus_wig)
            outputToBG(neg_depth_sorted, chrom, chrom_length, minus_wig)
    
        contig_sizes.append([chrom, chrom_length])

        del pos_depth_sorted
        del neg_depth_sorted

    if options.format == "bed":
        bedfile.close()
    else:
        plus_wig_name = plus_wig.name
        minus_wig_name = minus_wig.name
        plus_wig.close()
        minus_wig.close()

    outname_plus = args[0] + "_plus"
    outname_minus = args[0] + "_minus"

    if options.format == "bedgraph":
        E.debug("Outputting to bedGraph")
        shutil.move(plus_wig_name, outname_plus + ".bg")
        shutil.move(minus_wig_name, outname_minus + ".bg")
        
    elif options.format == "bigwig":
        chrom_sizes_file = tempfile.NamedTemporaryFile(delete=False, dir=".")
        contig_sizes = ["\t".join(map(str,row)) for row in contig_sizes]
        contig_sizes = "\n".join(contig_sizes) + "\n"
        chrom_sizes_file.write(contig_sizes)
        chrom_sizes_filename = chrom_sizes_file.name
        chrom_sizes_file.close()

        outputToBW(plus_wig_name, outname_plus, chrom_sizes_filename)
        outputToBW(minus_wig_name, outname_minus, chrom_sizes_filename)


    # write footer and output benchmark information.
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
