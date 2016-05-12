import numpy as np
import pandas as pd
import CGAT.GTF as GTF

from counting import count_transcript
from counting import count_intervals

##################################################
def bin_counts(counts, length, nbins):

    bins = np.linspace(0, length, num=nbins+1, endpoint=True)
    if len(counts.index.levels) == 2:
        bases = counts.index.droplevel()
    else:
        bases = counts.index

    binned_counts = counts.groupby(
        list(pd.cut(bases,
                    bins=bins,
                    labels=range(nbins),
                    include_lowest=True))).sum()

    binned_counts.index.name = "base"
    binned_counts = binned_counts.reindex(range(nbins), fill_value=0)
    return binned_counts


##################################################
def meta_gene(gtf_filelike, bam, bins=[10, 100, 10], flanks=100,
              output_matrix=False, calculate_flanks=False,
              pseudo_count=0):
    ''' Produce a metagene profile accross the :param gtf_file: from the reads
    in :param bam_file:.

        :type gtf_filelike: file or buffer
        :type bam: pysam.AlignmentFile
        :param output_matrix: if true, matrix of binned else None is returned

        :rtype: tuple of (pandas.Series, pandas.DataFrame)
    counts for each transcript is returned
    '''
    counts_collector = []
    regions = ["flank5", "exons", "flank3"]
    
    try:
        if flanks == 0:
            bins = bins[1]

        nbins_lookup = dict(zip(*(regions, bins)))
        length_lookup = dict(zip(*(regions, [flanks, 0, flanks])))

    except TypeError:
        if flanks > 0:
            nbins_lookup = {x: bins for x in bins}
   
    for transcript in GTF.transcript_iterator(
            GTF.iterator(gtf_filelike)):
      
        length = sum(
            [x.end - x.start for x in transcript if x.feature == "exon"])

        if calculate_flanks:
            flanks = length * float(bins[0])/bins[1]
            length_lookup = dict(zip(*(regions, [flanks, 0, flanks])))

        counts = count_transcript(transcript, bam, flanks=flanks)

        if flanks > 0:
            length_lookup["exons"] = length
            binned_counts = counts.groupby(level=0).apply(
                lambda x: bin_counts(x, length_lookup[x.name],
                                     nbins_lookup[x.name]))
        else:
            binned_counts = bin_counts(counts, length, bins)

        binned_counts.name = transcript[0].transcript_id

        counts_collector.append(binned_counts)

    counts_matrix = pd.concat(counts_collector, axis=1)
    counts_matrix = counts_matrix.transpose()
    counts_matrix = (counts_matrix.fillna(0) + pseudo_count).div(
        counts_matrix.sum(axis=1), axis=0)
    summed_matrix = counts_matrix.sum()
    summed_matrix.name = "density"

    if output_matrix:
        return summed_matrix, counts_matrix
    else:
        return summed_matrix, None


##################################################
def processing_index(interval_iterator, bam, window_size=50):
    '''Calculate the processing index for the speicied sample, using the
    provided interval_iterator to get the cleavage sites. The iterator
    can be GTF or BED, as long as it has end, contig and strand
    attributes. The end attribute will be used to define the
    cleavage site.

    The proccessing index for G genes is defined as:
    
    .. math::

       pi = log_2( \frac{\sum_{i=1}^{G} N_i^{PM}}{\sum_{i=1}^{G} N_i^M})

    after Baejen et al Mol Cell 5(55):745-757. However, Beaejen et al
    normalise this number to the total number of genes, which seems
    wrong to me. '''

    n_pm = 0
    n_m = 0
 
    for site in interval_iterator:
    
        if site.strand == "+":
            pos = site.end
        elif site.strand == "-":
            pos = site.start
        else:
            raise ValueError(
                "processing index not valid for unstranded cleavage points in "
                "entry\n" + str(site)+"\n")

        upstream_interval = (pos - window_size, pos)
        downstream_interval = (pos, pos + window_size)

        counts = count_intervals(bam, [upstream_interval, downstream_interval],
                                 site.contig, site.strand)

        # We are currently in genome cooridinates, not transcript

        if site.strand == "+":
            # pandas indexing is inclusive
            n_up = counts.iloc[:pos-1].sum()
            n_down = counts.iloc[pos:].sum()
        elif site.strand == "-":
            n_up = counts.iloc[pos:].sum()
            n_down = counts.iloc[:pos-1].sum()

        n_pm += n_down
        n_m += n_up-n_down

    pi = np.log2(float(n_pm)/float(max(1, n_m)))

    return pi


##################################################
def get_binding_matrix(bamfile,
                       gtflike_iterator,
                       align_at=0,
                       bin_size=25,
                       left_margin=500,
                       right_margin=10000):
    ''' Returns a matrix containing the binding densities across all the
    genes in the iterator binned into equal size bins of :param bin_size:
    base pairs.  Transcripts/genes are collapsed to just their exons, but
    sequence up and down stream of the ends of the transcript/gene is
    added. Unclipped transcripts are skipped.

      :param bam: BAM file with iCLIP reads
      :param gtflike_iterator: iterator returning bundles of CGAT.GTF
                  objects.  Note overlaping exons will be merged.
      :param align_at: Position in transcript at which to align profiles
                  from each transcript. The 0 default mean the start of the
                  gene/transcript.  Can be int or a Series index on
                  transcript_id with different in alignment positions
                  for each gene/transcript. 
      :param bin_size: Number of bases to include in a single bin. Depth
                  will be summed across bins
      :param left_margin: Number of bases to include to the left of the
                   alignment point
      :param right_margin: Ditto but to the left.

      :rtype: A dataframe with columns being bins and rows being genes. '''

    matrix = []
    for gene in gtflike_iterator:
        length = sum(e[1] - e[0] for e in GTF.asRanges(gene, "exon"))
        
        try:
            align_this_at = align_at[gene[0].transcript_id]
        except TypeError:
            align_this_at = align_at

        flanks = max(left_margin - align_this_at,
                     right_margin - (length-align_this_at))
        counts = count_transcript(gene, bamfile, flanks=flanks)

        # skip unclipped transcripts
        if counts.sum() == 0:
            continue

        counts = counts.reset_index()
    
        flank3 = counts["region"] == "flank3"
        flank5 = counts["region"] == "flank5"
    
        counts.loc[flank5, "base"] = counts.loc[flank5, "base"] - flanks
        counts.loc[flank3, "base"] = counts.loc[flank3, "base"] + length

        counts["base"] = counts["base"] - align_this_at
        counts["base"] = (counts["base"] // bin_size) * bin_size

        counts = counts.groupby("base").sum()[0]

        counts.name = gene[0].transcript_id
        matrix.append(counts)

    matrix = pd.concat(matrix, axis=1)
    matrix = matrix.reindex(range(-1*left_margin, right_margin, bin_size),
                            fill_value = 0)
    matrix = matrix.T.fillna(0)

    return matrix


def quantile_row_norm(matrix, quantile=0.99, min_value='auto'):
    '''Normalise the rows of a matrix so that 1 corresponds to the given
    quantile. :param min_value: allows a value to be provided so that
    the minimum normalisation value will never be smaller than this.
    If set to 'auto' then the smallest non-zero value across the whole
    matrix will be found and used'''

    if min_value=='auto':
        min_value = matrix[matrix > 0].min().min()

    normed_matrix = matrix.apply(lambda x: x/max(min_value, x.quantile(quantile)),
                                 axis=1)

    return normed_matrix


def sum_row_norm(matrix):
    '''Normalise the rows of a matrix by the sum of the row'''

    normed_matrix = matrix.apply(lambda x: x/x.sum())

    return normed_matrix


def compress_matrix(matrix, nrows=None, ncols=None):
    '''Compress a matrix to a new number of rows/columns. Cells
    in the matrix are collapsed by averaging. Assumes matrix is sorted'''

    if ncols:
        groups, bins = pd.cut(matrix.columns.values, ncols,
                              retbins=True, labels=False)
        groups = bins[groups]
        matrix = matrix.groupby(groups, axis=1).mean()

    if nrows:
        groups = pd.cut(range(matrix.shape[0]), nrows, labels=False)
        matrix = matrix.groupby(groups).mean()

    return matrix

