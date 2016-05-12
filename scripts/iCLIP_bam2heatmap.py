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
import CGAT.Experiment as E
import iCLIP
import pysam
import pandas

sort_choices = ["length", "3utr", "5utr", "none"]
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


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--gtf-file", dest="gtf", type="string",
                      help="GTF containing gene annotations")
    parser.add_option("-s", "--sort", dest="sort", type="choice",
                      default="none",
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
    parser.add_option("-a", "--align-at", dest="align_at", type=choice,
                      default="start",
                      choices=align_choices,
                      help="Where to align genes/transcripts at. Choices are %s"
                            % ", ".join(align_choices))
    parser.add_option("-h", "--height", dest="height", type="int",
                      default=100,
                      help="Number of rows in output matrix/heigh of plot in px")
    parser.add_option("-w", "--width", dest="width", type="int",
                      default=None,
                      help="Number of columns in output/width of plot in px"
                           "default based on bin size")
    parser.add_option("-n", "--normalize", dest="normalize", type=choice,
                      default="None",
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
                      help="Add annotations to the output plot")
    parser.add_option("--reverse-strand", dest="rstrand", action="store_true",
                      default=False,
                      help="Find reads on reverse strand")
    parser.add_option("-f", "--feature", dest="feature", type=choice,
                      choices=["gene", "transcript"],
                      default="gene",
                      "use genes or transcripts")


    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.gtf:
        
        if options.feature == "gene":
            gtf_iterator = GTF.flat_gene_iterator(GTF.iterator(options.gtf))
        else:
            gtf_iterator = GTF.transcript_iterator(GTF.iterator(options.gtf))

        lengths = dict()
        utr3_lengths = dict()
        utr5_lengths = dict()

        for transcript in gtf_iterator:
            lengths[transcript[0].transcript_id] = sum(
                [e[1] - e[0] for exon in GTF.asRanges(transcript, "exon")])

            utrs = GTF.asRanges(transcript, "UTR")
            coding = Intervals.truncate(GTF.asRanges(transcript, "exon"), utrs)
            coding.sort()

            utr5 = [utr for utr in utrs if utr[1] <= coding[0][0]]
            utr3 = [utr for utr in turs if utr[0] >= coding[-1][-1]]

            if transcript[0].strand == "-":
                utr3, utr5 = utr5, utr3

            utr3_lengths[transcript[0].transcript_id] = sum(
                [e[1] - e[0] for e in utr3])

            utr5_lengths[transcript[0].transcript_id] = sum(
                [e[1] - e[0] for e in utr5])

        lengths = pandas.Series(lengths)
        utr3_lengths = pandas.Series(utr3_lengths)
        utr5_lengths = pandas.Series(utr5_lengths)

    if options.use_matrix:
        raw_matrix = pandas.read_csv(options.use_matrix,
                                     sep="\t",
                                     index_col=1)
    else:

        try:
            bamfile = pysam.AlignmentFile(args[0])
        except IOError:
            E.error("Cannot open bamfile %s" % args[0])
            return(1)
        except IndexError:
            E.error("No bamfile specified")
            print(globals()["__usage__"])
            return(1)

        if options.feature == "gene":
            gtf_iterator = GTF.flat_gene_iterator(GTF.iterator(options.gtf))
        else:
            gtf_iterator = GTF.transcript_iterator(GTF.iterator(options.gtf))

        if options.ds_win is None:
            options.ds_win = lengths.max()

        if options.align_at == "start":
            align_at = 0
        elif options.align_at == "end":
            align_at = lengths
            options.ds_win, options.us_win = options.us_win, options.ds_win
            
        if options.rstrand:
            def _it_reverse(gtf):
                for transcript in gtf:
                    if transcript[0].strand == "+":
                        transcript[0] == "-"
                    else:
                        transcript[0] == "+"
                    yield transcript

            gtf_iterator = it_reverse(gtf_iterator)
            options.ds_win, options.us_win = options.us_win, options.ds_win

        raw_matrix = iCLIP.get_binding_matrix(bamfile, gtf_iterator,
                                              align_at=align_at,
                                              bin_size=options.bin_size,
                                              left_margin=options.us_win,
                                              right_margin=options.ds_win)
        if options.rstrand:
            raw_matrix.columns = -1 * raw_matrix.columns.values

    normalized_matrix = normalize(raw_matrix, options.normalize,
                                  quantile=options.quantile)

    if options.sort == "length":
        sorter = lengths
    elif options.sort == "3utr":
        sorter = utr3_lengths
    elif options.sort == "5utr":
        sorter = utr5_lengths
    elif options.sort == "none":
        sorter = pandas.Series(range(raw_matrix.shape[0]),
                               index=raw_matrix.index)

    sorter = sorter[sorter.index.isin(normalized_matrix.index)]
    sorter = sorter.sort_values()
    sorted_matrix = normalized_matrix.loc[sorter.index.values]

    compress_matrix = iCLIP.compress_matrix(sorted_matrix,
                                            ncols=options.ncols,
                                            nrows=options.nrows)

    renormalized_matrix = normalize(compress_matrix, options.renormalize,
                                    qunatile=options.quantile)

    if renormalized_matrix == raw_matrix and options.use_matrix is not None:
        E.info("Input and output matrices are identical, no matrix output")
    else:
        if options.outfile_pattern:
            mat_outfile = IOTools.openFile(
                options.outfile_pattern % "matrix.tsv.gz", "w")
        else:
            mat_outfile = options.stdout

        renormalized_matrix.to_csv(mat_outfile, sep="\t")

    if options.plot:

        try:
            from rpy2.robjects import r as R
        except:
            E.info("No rpy2. Not plotting image")
            return(0)

        from rpy2.robjects.numpy2ri import numpy2ri
        ro.conversion.py2ri = numpy2ri
        rpy2.robjects.numpy2ri.activate()

        if options.outfile_pattern:
            plot_outfile = options.outfile_pattern % "png", "w"
        else:
            plot_outfile = "bam2heatamp_out.png"

        c = R["c"]

        R["png"](plot_outfile,
                 width=renormalized_matrix.shape[1] + 72,
                 height=renormalized_matrix.shape[0] + 72,
                 res=72)
        R.par(mai=c(1, 0.5, 0, 0.5))
        cols = R["colorRampPallette"](c("white", "blue"))(50)
        bases = renormalized_matrix.columns.values.astype("int")
        groups = renormalized_matrix.index.values.astype("int")
        mat = renormalized_matrix.as_matrix()
        mat[mat > 0.99999] = 0.99999
        R.image(bases, groups, R.t(mat),
                zlim=c(0, 1),
                raster=True,
                col=cols)
        R["dev.off"]()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
