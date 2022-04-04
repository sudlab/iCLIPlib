'''
filter_internal_priming.py - Use Quant-seq read2s to locate polyA sites
====================================================
:Author: Ian Sudbery, Mikayla Fernholz
:Release: $Id$
:Date: |today|
:Tags: Python
Purpose
-------
Filter Quant-seq peak calling bam files for internal priming

Usage
-----
python filter_internal_priming.py BAMFILE [OPTIONS]
Type::
   python cgat_script_template.py --help
for command line help.
Command line options
--------------------
'''

import pysam
from cgat.Genomics import reverse_complement
from email.policy import default
import sys
import cgatcore.experiment as E
import cgat.Bed as Bed
import os
import iCLIP



def filter_internal_priming(bamfile,
                            genome,
                            out_name,
                            allow_mismatch = False):
    ''' This function filters a bam file for internal priming as part of polyA
    seq '''

    outfile = pysam.AlignmentFile(out_name,
                                 mode="wb",
                                 template=bamfile)


    last_base_T = 0
    start_not_aligned = 0
    mismatchs = 0
    motif_in_flank = 0
    too_many_As = 0
    run_of_As = 0
    output = 0
    input_reads = 0

    for read in bamfile.fetch(until_eof=True):

        #if input_reads >= 500:
        #    break

        input_reads += 1

        ##Check if last nt is A
        read_sequence = read.get_forward_sequence()

        if read_sequence[0] == "T":
            last_base_T += 1
            continue

        ##Check if last 3nt map perfectly
        if allow_mismatch:
            cigar = read.cigartuples

            first_operation = cigar[0]

            if not (first_operation[0] == 0
                    and first_operation[1] >= 3):
                start_not_aligned += 1
                #print (read)
                continue

            if not read.is_reverse:
                sequence = genome.fetch(read.reference_name,
                                        read.reference_start,
                                        read.reference_start +3)

                if not sequence == read_sequence[0:3]:
                    mismatchs += 1
                    continue

            else:

                sequence = genome.fetch(read.reference_name,
                                        read.reference_end -3,
                                        read.reference_end)

                sequence = reverse_complement(sequence)

                if not sequence == read_sequence[0:3]:
                    mismatchs += 1
                    continue

        ##Get sequence of flanking 10nt
        if not read.is_reverse:
            start = read.reference_start - 10
            end = read.reference_start

            sequence = genome.fetch(reference=read.reference_name,
                                    start = start,
                                    end = end)
            sequence = reverse_complement(sequence)

        else:

            start = read.reference_end
            end = read.reference_end + 10

            sequence = genome.fetch(reference=read.reference_name,
                                    start = start,
                                    end = end)

        ##Check for AAAA, AGAA, AAGA and AAAG
        if sequence.startswith("AAAA"):
            motif_in_flank += 1
            continue

        if sequence.startswith("AGAA"):
            motif_in_flank += 1
            continue

        if sequence.startswith("AAGA"):
            motif_in_flank += 1
            continue

        if sequence.startswith("AAAG"):
            motif_in_flank += 1
            continue

        ##Check total A content
        if sequence.count("A") >= 7:
            too_many_As += 1
            continue

        ##Check < 6 consecurive As
        if "AAAAAA" in sequence:
            run_of_As += 1
            continue

        outfile.write(read)
        output += 1
    outfile.close()

    number_filtered = ["Number of input reads:", str(input_reads),
    ", Number of passed reads:", str(output)]
    #E.debug(number_filtered)
    return(number_filtered)

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

    parser.add_option("-M", "--dont_allow_mismatch", dest="dont_allow_mismatch",
                      action="store_true",
                      default=False,
                      help="Check and filter out reads with mismatches at the last 3nt. Default: Allow mismatches in the last 3nt")
    parser.add_option(
        "-f", "--reference-fasta-file", dest="reference_fasta_file",
        help="reference genomic sequence in fasta format. ")


    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    bamfile = pysam.AlignmentFile(args[0])
    genome = pysam.FastaFile(options.reference_fasta_file)
    name_outfile = os.path.abspath(args[0])[:-4] + ".filtered.bam"
    num_filter = filter_internal_priming(bamfile,
                                        genome = genome,
                                        out_name = name_outfile,
                                        allow_mismatch = options.dont_allow_mismatch)
    #E.debug(num_filter)
    name_summary= os.path.abspath(args[0])[:-4] + ".summary"
    file_number_filtered = open(name_summary,'w')
    for element in num_filter:
        file_number_filtered.write(element)
        file_number_filtered.write(' ')
    file_number_filtered.close()

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
