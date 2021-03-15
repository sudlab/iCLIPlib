'''
iCLIP_bam2heatmap.py - Produce sorted heatmaps of binding activity
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Overall purpose and function of the script>

Usage
-----

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
import os
import cgatcore.experiment as E
from cgat import GTF
from cgatcore import iotools
from cgat import Intervals
import pysam
import pandas

sys.path.insert(1, os.path.join(
    os.path.dirname(__file__), ".."))

import iCLIP

sort_choices = ["length", "first-exon", "3utr", "5utr", "manual", "none"]
align_choices = ["start", "end"]
norm_choices = ["quantile", "sum", "none"]
annotation_choices = ["start", "end", "3utr", "5utr"]


def normalize(matrix, method, quantile=0.99):

    if method == "quantile":
        return iCLIP.quantile_row_norm(matrix, quantile=quantile)
    elif method == "sum":
        return iCLIP.sum_row_norm(matrix)
    else:
        return matrix


def get_matrix(getter, lengths, options):

    if getter is None:
        E.error("No bamfile or wigfile specified")
        print((globals()["__usage__"]))
        return(1)

    f = iotools.open_file(options.gtf)
    if options.feature == "gene":
        gtf_iterator = GTF.flat_gene_iterator(GTF.iterator(f))
    else:
        gtf_iterator = GTF.transcript_iterator(GTF.iterator(f))

    if options.ds_win is None:
        ds_win = lengths.max()
    else:
        ds_win = options.ds_win

    if options.align_at == "start":
        align_at = 0
        us_win, ds_win = options.us_win, ds_win
    elif options.align_at == "end":
        align_at = lengths
        ds_win, us_win = options.us_win, ds_win
            
    if options.rstrand:
        def _it_reverse(gtf):
            for transcript in gtf:
                if transcript[0].strand == "+":
                    transcript[0].strand = "-"
                else:
                    transcript[0].strand = "+"
                yield transcript

        gtf_iterator = _it_reverse(gtf_iterator)
        ds_win, us_win = us_win, ds_win
        align_at = lengths

    raw_matrix = iCLIP.get_binding_matrix(getter, gtf_iterator,
                                          align_at=align_at,
                                          bin_size=options.bin_size,
                                          left_margin=us_win,
                                          right_margin=ds_win)
    if options.rstrand:
        
        raw_matrix.columns = -1 * raw_matrix.columns.values
        raw_matrix = raw_matrix.sort_index(axis=1)

    return raw_matrix

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

    parser.add_option("-g", "--gtf-file", dest="gtf", type="string",
                      help="GTF containing gene annotations")
    parser.add_option("-s", "--sort", dest="sort", type="choice",
                      default="length",
                      choices=sort_choices,
                      help="Property to sort rows by. Choices are %s"
                           % ", ".join(sort_choices))
    parser.add_option("-b", "--bin-size", dest="bin_size", type="int",
                      default=25,
                      help="Size of window over which to sum reads")
    parser.add_option("-u", "--upstream-window", dest="us_win", type="int",
                      default=500,
                      help="Amount of sequence upstream of alignment point (less introns)")
    parser.add_option("-d", "--downstream-window", dest="ds_win", type="int",
                      default=None,
                      help="Amount of sequence downstream of alignment point (default longest segment)")
    parser.add_option("-a", "--align-at", dest="align_at", type="choice",
                      default="start",
                      choices=align_choices,
                      help="Where to align genes/transcripts at. Choices are %s"
                            % ", ".join(align_choices))
    parser.add_option("-H", "--height", dest="height", type="int",
                      default=None,
                      help="Number of rows in output matrix/heigh of plot in px")
    parser.add_option("-w", "--width", dest="width", type="int",
                      default=None,
                      help="Number of columns in output/width of plot in px"
                           "default based on bin size")
    parser.add_option("-n", "--normalize", dest="normalize", type="choice",
                      default="none",
                      choices=norm_choices,
                      help="Row normalization to apply. Choices are: %s"
                           % ", ".join(norm_choices))
    parser.add_option("-r", "--renormalize", dest="renormalize", type="choice",
                      default="none",
                      choices=norm_choices,
                      help="Row normalization to apply after row/column compression")
    parser.add_option("--no-plot", dest="plot", action="store_false",
                      default=True,
                      help="Do not output plot - compute matrix only")
    parser.add_option("--use-matrix", dest="use_matrix", type="string",
                      default=None,
                      help="Use existing matrix")
    parser.add_option("--annotations", dest="annotations", type="choice",
                      action="append",
                      choices=annotation_choices,
                      help="Add annotations to the output plot")
    parser.add_option("--reverse-strand", dest="rstrand", action="store_true",
                      default=False,
                      help="Find reads on reverse strand")
    parser.add_option("-f", "--feature", dest="feature", type="choice",
                      choices=["gene", "transcript"],
                      default="gene",
                      help="use genes or transcripts")
    parser.add_option("--quantile", dest="quantile", type="float",
                      default=0.99,
                      help="Quantile to use in quantile normalization")
    parser.add_option("-o", "--outfile-prefix", dest="outfile_pattern", type="string",
                      default=None,
                      help="base of names for output files")
    parser.add_option("-c", "--crop", dest="crop", type="string",
                      default=None,
                      help="crop view to a certain range on the xaxis. Specify like"
                      "-500:1000")
    parser.add_option("--format", dest="format", type="string",
                      default="png",
                      help="Output format, use valid R graphics device")
    parser.add_option("--plus-wig", dest="plus_wig", type="string",
                      help="Use this wig for plus strand info rather than bam file")
    parser.add_option("--minus-wig", dest="minus_wig", type="string",
                      help="Use this wig for minus strand info rather than bam file")
    parser.add_option("--bed", dest="bed", type="string",
                      help="Use this bed for signal(must be indexed)")
    parser.add_option("--norm-mat", dest="norm_mat", type="string",
                      help="Use this matrix for normalizing (e.g. RNA data")
    parser.add_option("--sort-order-file", dest="sort_file", type="string",
                      default=None,
                      help="Two column file containing gene names in the first"
                           "column and a numeric value to sort on in the second")
    parser.add_option("-p", "--profile", dest="profile", type="choice",
                      choices=profiles,
                      default="iclip",
                      help="Experiment profile to use. Sets various things"
                      "about obtaining 1-bp position from read. Options are"
                      " %s" % ", ".join(profiles))                       

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    if options.plot and (options.height is None):
        options.height = 100

    if options.gtf:
        
        f = iotools.open_file(options.gtf)
        if options.feature == "gene":
            gtf_iterator = GTF.flat_gene_iterator(GTF.iterator(f))
        else:
            gtf_iterator = GTF.transcript_iterator(GTF.iterator(f))

        lengths = dict()
        utr3_lengths = dict()
        utr5_lengths = dict()
        first_exon_lengths = dict()
        for transcript in gtf_iterator:
            lengths[transcript[0].transcript_id] = sum(
                [e[1] - e[0] for e in GTF.asRanges(transcript, "exon")])

            exons = GTF.asRanges(transcript, "exon")
            utrs = GTF.asRanges(transcript, "UTR")
            coding = Intervals.truncate(exons, utrs)
            coding.sort()

            utr5 = [utr for utr in utrs if utr[1] <= coding[0][0]]
            utr3 = [utr for utr in utrs if utr[0] >= coding[-1][-1]]

            if transcript[0].strand == "-":
                utr3, utr5 = utr5, utr3
            
            if transcript[0].strand == "+" or len(exons) == 1:
                first_exon_lengths[transcript[0].transcript_id] = \
                    exons[0][1] - exons[0][0]
            else:
                first_exon_lengths[transcript[0].transcript_id] = \
                    exons[-1][1] - exons[-1][0]

            utr3_lengths[transcript[0].transcript_id] = sum(
                [e[1] - e[0] for e in utr3])

            utr5_lengths[transcript[0].transcript_id] = sum(
                [e[1] - e[0] for e in utr5])

        lengths = pandas.Series(lengths)
        utr3_lengths = pandas.Series(utr3_lengths)
        utr5_lengths = pandas.Series(utr5_lengths)
        first_exon_lengths = pandas.Series(first_exon_lengths)

    else:
        options.sort = "none"
        options.annotations = None

    if options.plus_wig:
        getter = iCLIP.make_getter(plus_wig=options.plus_wig,
                                   minus_wig=options.minus_wig)
    elif options.bed:
        getter = iCLIP.make_getter(bedfile=options.bed)
    else:
        try:
            getter = iCLIP.make_getter(bamfile=args[0], profile=options.profile)
        except IOError:
            E.error("Cannot open bamfile %s" % args[0])
            return(1)
        except IndexError:
            getter = None

    if options.use_matrix:
        raw_matrix = pandas.read_csv(options.use_matrix,
                                     sep="\t",
                                     index_col=0)
        raw_matrix.columns = raw_matrix.columns.astype("int")
    else:
        raw_matrix = get_matrix(getter, lengths, options)

    if options.crop:
        crop_from, crop_to = list(map(int, options.crop.split(":")))
        raw_matrix = raw_matrix.loc[:, crop_from:crop_to]

    if options.norm_mat:
        norm_matrix = pandas.read_csv(options.norm_mat,
                                      sep="\t",
                                      index_col=0)
        norm_matrix.columns = norm_matrix.columns.astype("int")

        if options.crop:
            norm_matrix = norm_matrix.loc[:, crop_from:crop_to]
        
        if all(norm_matrix.columns == raw_matrix.columns) and \
           all(raw_matrix.index.isin(norm_matrix.index.values)):
            norm_matrix = norm_matrix.loc[raw_matrix.index]
            norm_matrix = norm_matrix.replace(
                0, norm_matrix[norm_matrix > 0].min().min())
            raw_matrix = raw_matrix/norm_matrix
            norm_matrix = None

        else:
            raise ValueError("Incompatible normalisation matrix")

    normalized_matrix = normalize(raw_matrix, options.normalize,
                                  quantile=options.quantile)

    if options.sort == "length":
        sorter = lengths
    elif options.sort == "3utr":
        sorter = utr3_lengths
    elif options.sort == "5utr":
        sorter = utr5_lengths
    elif options.sort == "first-exon":
        sorter = first_exon_lengths
    elif options.sort == "manual":
        sorter = pandas.read_csv(options.sort_file, sep="\t",
                                 index_col=0, usecols=[0, 1])
        sorter = sorter[sorter.columns[0]]
    elif options.sort == "none":
        sorter = pandas.Series(range(raw_matrix.shape[0]),
                               index=raw_matrix.index[::-1])

    sorter = sorter[sorter.index.isin(normalized_matrix.index)]
    sorter = sorter.sort_values(ascending=False)
    sorted_matrix = normalized_matrix.loc[sorter.index.values]

    compress_matrix = iCLIP.compress_matrix(sorted_matrix,
                                            ncols=options.width,
                                            nrows=options.height)

    renormalized_matrix = normalize(compress_matrix, options.renormalize,
                                    quantile=options.quantile)

    if renormalized_matrix is raw_matrix and options.use_matrix is not None:
        E.info("Input and output matrices are identical, no matrix output")
    else:
        if options.outfile_pattern:
            mat_outfile = iotools.open_file(
                options.outfile_pattern + ".matrix.tsv.gz", "w")
        else:
            mat_outfile = options.stdout

        renormalized_matrix.to_csv(mat_outfile, sep="\t")

    if options.plot:

        try:
            from rpy2.robjects import r as R
            from rpy2 import robjects as ro
        except ImportError:
            E.info("No rpy2. Not plotting image")
            return(0)

        from rpy2.robjects.numpy2ri import numpy2ri
        ro.conversion.py2ri = numpy2ri
        ro.numpy2ri.activate()

        if options.outfile_pattern:
            plot_outfile = options.outfile_pattern + ".png"
        else:
            plot_outfile = "bam2heatmap_out.png"

        c = R["c"]

        R[options.format](plot_outfile,
                          width=renormalized_matrix.shape[1] + 72,
                          height=renormalized_matrix.shape[0] + 72,
                          unit="px",
                          res=72)
        R.par(mai=c(1, 0.5, 0, 0.5))
        cols = R["colorRampPalette"](c("white", "blue"))(50)
        bases = renormalized_matrix.columns.values.astype("int")
        groups = renormalized_matrix.index.values.astype("int")
        mat = renormalized_matrix.as_matrix()
        mat[mat >= 1] = 1

        R.image(bases, groups, R.t(mat),
                zlim=c(0, 1),
                raster=True,
                col=cols,
                xlab="Base",
                yaxt="n")

        def _sort_and_compress_annotation(anno):
            sorted_anno = anno.loc[sorter.index]
            comp_anno = iCLIP.compress_matrix(
                sorted_anno, renormalized_matrix.shape[0])
            return comp_anno

        if options.annotations:
            ends = _sort_and_compress_annotation(lengths)
            starts = pandas.Series(0, index=renormalized_matrix.index)

            if options.align_at == "end":
                starts, ends = -1 * ends, starts

            if "start" in options.annotations:
                R.lines(starts.values, starts.index.values, col="black", pch=".")
            if "end" in options.annotations:
                R.lines(ends.values, ends.index.values,
                        pch=".", col="black")
            if "5utr" in options.annotations:
                utr5s = _sort_and_compress_annotation(utr5_lengths)
                utr5s = starts + utr5s
                R.lines(utr5s.values, utr5s.index.values, col="orange", pch=".")
            if "3utr" in options.annotations:
                utr3s = _sort_and_compress_annotation(utr3_lengths)
                utr3s = ends - utr3s
                R.lines(utr3s.values, utr3s.index.values, col="orange", pch=".")

        R["dev.off"]()

    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
