'''
call_polyA_sites.py - Use Quant-seq read2s to locate polyA sites
====================================================

:Author: Ian Sudbery, Charlotte Vandermeulen, Mikayla Fernholz
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Based on 3P-seq analysis described here:
http://bartellab.wi.mit.edu/protocols/3PSeqV2.pdf

Usage
-----

python call_polyA_sites.py BAMFILE [OPTIONS]


Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

from email.policy import default
import sys
import cgatcore.experiment as E
import cgat.Bed as Bed
import os
import pysam

import iCLIP

def find_peaks_in_contig(getter,
                         contig,
                         strand,
                         window_width=10,
                         threshold_height=2,
                         threshold_sum=2,
                         dtype="int32"):
    '''This function performs the main algo. It scans the specified contig
    finding the maximum point, records its position and signal, and then
    groups it with the `width` bases either side, and removes them from the
    signal. It returns a list of "peaks", each a list with the  with the
    following fields:

    0: contig
    1: start (0-based)
    2: end (half-open)
    3: peak_height
    4: sum of signal in peak
    5: strand'''

    results = list()
    contig0 = getter(contig, strand=strand, dtype=dtype)
    contig_sum = contig0.sum()
    last_index = contig0.last_valid_index()

    while contig_sum > 0:
        peak_location = contig0.idxmax()
        #E.debug(print(peak_location))
        peak_start = peak_location - window_width
        peak_end = peak_location + window_width + 1
        peak = contig0.loc[peak_start:peak_end]
        peak_height = contig0.loc[peak_location]
        peak_sum = peak.sum()
        contig0[peak_start:peak_end] = 0
        contig_sum = contig_sum - peak_sum

        if (peak_start < 0):
            peak_start = 0
        if (peak_end > last_index):
            E.debug(peak_end)
            peak_end = last_index + 1
            E.debug(peak_end)

        if (peak_height <= threshold_height) or (
            peak_sum <= threshold_sum):
            continue

        results.append([contig, int(peak_start), int(peak_end),
                        peak_height, peak_sum, strand])

    return (results)

def output_peak_to_bed(outfile, peak):
    '''Takes a peak line as output by find_peaks_in_contig and
    converts to bed format and outputs to the provided file'''

    outbed = Bed.Bed()
    outbed.contig = peak[0]
    outbed.start = peak[1]
    outbed.end = peak[2]
    outbed["name"] = peak[3]
    outbed["score"] = peak[4]
    outbed["strand"] = peak[5]

    outfile.write(str(outbed) + "\n")

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

    parser.add_option("-H", "--peak-height", dest="peak_height", type="int",
                      default=2,
                      help="Minimum signal at highest point to call a peak")
    parser.add_option("--peak-sum", dest="peak_sum", type="int",
                      default=2,
                      help="Minimum number of reads in window to call peak")
    parser.add_option("-w", "--peak-width", dest="peak_width", type="int",
                      default=10,
                      help="Number of bases to include either side of highest"
                           "point in peak")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args[0])
    polyA_signal = iCLIP.make_getter(bamfile=bamfile,
                                use_deletions=False,
                                 offset=0,
                                 filter_end="read1")

    contigs = bamfile.references

    out_peaks = list()
    for contig in contigs:
        results_plus = find_peaks_in_contig(polyA_signal, contig, "+",
                                       window_width=options.peak_width,
                                       threshold_height=options.peak_height,
                                       threshold_sum=options.peak_sum)
        results_minus = find_peaks_in_contig(polyA_signal, contig, "-",
                                       window_width=options.peak_width,
                                       threshold_height=options.peak_height,
                                       threshold_sum=options.peak_sum)
        out_peaks.extend(results_plus)
        out_peaks.extend(results_minus)
        #E.debug(
        #    "Done with contig %s, found %i peaks" %
        #    (contig, len(results_plus + results_minus)))

    #E.debug("Starting output....")

    for peak in out_peaks:
        output_peak_to_bed(options.stdout, peak)

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
