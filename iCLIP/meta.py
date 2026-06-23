"""Functions in this module produce profiles of binding across sets of
genes/transcripts, or provide functions for normalizing/aggregating
such profiles.

"""


import functools

import numpy as np
import pandas as pd
import cgat.GTF as GTF
import time
import logging

from .counting import count_transcript
from .counting import count_intervals


##################################################
def bin_counts(counts, length, nbins):
    """Aggregate counts into specified number of bins spatial bins.

    Parameters
    -----------
    counts : pandas.Series
        Series of counts to be binned. Index is position, value is count.
    length : int
        Total length of sequence to be binned. The range 0...`length` is
        divided into equally sized bins.
    nbins : int
        Number of bins to divide the range 0...`length` into.

    Returns
    -------
    pandas.Series
        Value of returned Series is the summed counts over all bases that
        fall into a particular bin.

    """

    bins = np.linspace(0, length, num=nbins+1, endpoint=True)
    if isinstance(counts.index, pd.core.indexes.multi.MultiIndex):
        bases = counts.index.droplevel()
    else:
        bases = counts.index

    if counts.sum() != 0:

        binned_counts = counts.groupby(
            list(pd.cut(bases,
                        bins=bins,
                        labels=range(nbins),
                        include_lowest=True))).sum()
        binned_counts.index.name = "base"
        binned_counts = binned_counts.reindex(range(nbins), fill_value=0)
    else:
        binned_counts = pd.Series([0]*nbins, index=range(nbins))

    return binned_counts


##################################################
def meta_gene(gtf_filelike, bam, bins=[10, 100, 10], flanks=100,
              output_matrix=False, calculate_flanks=False,
              pseudo_count=0, row_norm=True):
    ''' Produce a metagene across a gtf file from CLIP data, where each
    gene is divided into the same number of bins.

    Parameters
    ----------
    gtffile_like : file or buffer or similar
        Handle to a GTF file containing annotations to count across
    bam : func
        getter function, as created by :func:`make_getter`, from which to
        retrieve cross-link data
    bins : sequence of int, optional
        Length 3 sequence containing the number of bins in the 5' flank,
        the body of each transcripts and the 3' flank.
    flanks : int, optional
        Length, in bp of the flanking regions.
    output_matrix : False, optional
        Output the whole (#transcripts, #bins) matrix alongside the
        averaged profile.
    calculate_flanks : bool, optional
        Calculate the size of the flanks seperately for each transcript so
        that each flank bin is the same width as each transcript body bin.
    pseudo_count : int, optional
        pseudocount to add to each bin, , minimizing the effect of
        transcripts with a single or small number of cross-linked bases in
        a `row_norm`-ed profile.
    row_norm : bool, optional
        Divide counts in each bin by the row-sum for each transcript?

    Returns
    -------
    pandas.Series of float
        Averaged profile over all transcripts. Will have
        `pandas.MultiIndex` with first level corresponding to the region
        ["flank5", "exons", "flank3"]. Second level is the bin within the
        regions.


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

            binned_counts = counts.groupby(level=0, dropna=False).apply(
                lambda x: bin_counts(x, length_lookup[x.name],
                                     nbins_lookup[x.name]))
        else:
            binned_counts = bin_counts(counts, length, bins)

        binned_counts.name = transcript[0].transcript_id

        counts_collector.append(binned_counts)

    counts_matrix = pd.concat(counts_collector, axis=1)

    counts_matrix = counts_matrix.transpose()

    counts_matrix = counts_matrix.fillna(0)
    counts_matrix[counts_matrix.sum(axis=1) > 0] = \
        counts_matrix[counts_matrix.sum(axis=1) > 0] + pseudo_count

    if row_norm:
        counts_matrix = counts_matrix.div(
            counts_matrix.sum(axis=1), axis=0)

    summed_matrix = counts_matrix.sum()
    summed_matrix.name = "density"

    if output_matrix:
        return summed_matrix, counts_matrix
    else:
        return summed_matrix, None


##################################################
def processing_index(interval_iterator, bam, window_size=50):
    '''Calculate the ratio of processed transcripts to non-processed

    Parameters
    ----------
    interval_iterator : CGAT.Bed or CGAT.GTF-like iterator
        The iterator must yeild objects that have a start, end and strand
        attribute. Processing index will be calculated around these.
    bam : *_getter-like function
        A getter function returned by the `make_getter` function, this will
        be used to retrieve cross-link counts.
    window_size : int, optional
        How far up and downstream of the the processing site to consider.

    Returns
    -------
    float
        processing index averaged over all processing sites given.
    '''

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
    '''Get matrix containing the binding counts across all the
    genes in the iterator binned into equal sized bins

    Transcripts/genes are collapsed to just their exons, but
    sequence up and down stream of the ends of the transcript/gene is
    added. Zero count transcripts/genes are excluded.

    Parameters
    ----------
    bam : *_getter func
        Function to access cross-links from, as returned by
        :func:`make_getter`
    gtflike_iterator : iter of CGAT.GTF.Entry
        iterator returning sequences of CGAT.GTF.Entry objects. Note
        overlapping exons will be merged.
    align_at : int or pandas.Series, optional
        Position in transcript at which to align profiles form each
        transcript. The defaults means the start of the gene/transcript.
        If int same position is used for every gene/transcript. If
        `pandas.Series`, index should contain transcript_ids.
    bin_size : int, optional
        Number of bases to include in a single bin. Depth will be summed
        across these bins.
    left_margin, right_margin : int, optional
        Number of bases to include to the left and right of the alignment
        point.

    Returns
    -------
    pandas.Dataframe
        Columns are bins and row are genes.

    '''

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
    matrix = matrix.reindex(range(-1*(left_margin//25)*25, 
                                  (right_margin//25)*25,
                                   bin_size),
                            fill_value = 0)
    matrix = matrix.T.fillna(0)

    return matrix


def quantile_row_norm(matrix, quantile=0.99, min_value='auto'):
    '''Normalise the rows of a matrix so that 1 corresponds to the given
    quantile.

    :param min_value: allows a value to be provided so that
    the minimum normalisation value will never be smaller than this.
    If set to 'auto' then the smallest non-zero value across the whole
    matrix will be found and used

    Parameters
    ----------
    matrix : pandas.DataFrame
        Matrix of data for which the rows are to be normalized
    quantile : float, optional
        Quantile to use as the normalizer, defaults to 0.99 thus removing
        only the most outlying values.
    min_value : float or "auto"
        The minimum value the normalizer can take. This prevents rows being
        normalized by 0. If '``auto``' then for each for the smallest non
        zero value is used.

    Returns
    -------
    pandas.DataFrame
        Same row and column indexes as `matrix`. Rows are now normalized'''

    if min_value=='auto':
        min_value = matrix[matrix > 0].min().min()

    normed_matrix = matrix.apply(lambda x: x/max(min_value, x.quantile(quantile)),
                                 axis=1)

    return normed_matrix


def sum_row_norm(matrix):
    '''Normalise the rows of a matrix by the sum of the row

    Parameters
    ----------
    matrix : pandas.DataFrame
        Matrix over which to normalise rows

    Returns
    -------
    pandas.DataFrame
        Indexes will be identical to `matrix`

    '''

    normed_matrix = matrix.apply(lambda x: x/x.sum(), axis=1)

    return normed_matrix


def compress_matrix(matrix, nrows=None, ncols=None):
    '''Compress a matrix to a new number of rows/columns. Cells
    in the matrix are collapsed by averaging. Assumes matrix is sorted

    Parameters
    ----------
    matrix : pandas.DataFrame
        matrix to be compressed. Index(s) must be numeric.
    nrows, nocls : int
        Number of rows/columns in the output matrix. If ``None`` then will
        be same as input matrix.

    Returns
    -------
    pandas.DataFrame
        If `nrows`/`ncols` is ``None`` then the corresponding index will
        be unchanged. Otherwise, index will be integer 0 to `nrows`/`ncols`
        .
    '''

    if ncols and not ncols == matrix.shape[1]:
        groups, bins = pd.cut(matrix.columns.values, ncols,
                              retbins=True, labels=False)
        groups = bins[groups]
        matrix = matrix.groupby(groups, axis=1).mean()

    if nrows and not nrows == matrix.shape[0]:
        groups = pd.cut(range(matrix.shape[0]), nrows, labels=False)
        matrix = matrix.groupby(groups).mean()

    return matrix

def get_window(profile, position, upstream, downstream):
    '''Return a window around point with an index aligned
    so the specified point is 0

    Parameters
    ----------
    profile : pandas.Series
        Profile to be subset, as generated by a *_getter function
    position : int
        Position with the profile to extract a window around
    upstream : int
        Number of bases upstream of the `position` to extract
    downstream : int
        Number of bases downstream of `position` to extract

    Value
    -----
    pandas.Series
        The part of the input profile from ``position-upstream`` to 
        ``position + downstream``, indexed such that ``position`` is now 0.
    '''

    window = profile.loc[(position-upstream):(position+downstream)]
    window.index = window.index - position
    return window


def transcript_region_meta(transcript, getter, regions, names, bins, length_norm=True, aggregate=False):
    '''Return a profile of binding across a transcript, divided into
    specified regions. Each region is divided into the specified number of bins,
    and the counts in each bin are summed. If `aggregate` is True,
    then each exon within a region is binned separately and counts across all
    exons in that region are summed.

    If a region has no length, it will be excluded from the output profile. If
    all regions have no length, an empty profile with the correct MultiIndex
    will be returned to enable safe concatenation with other profiles.
    '''

    # basic validation
    if not (len(names) == len(regions) == len(bins)):
        raise ValueError("names, regions and bins must be same length")

    logger = logging.getLogger(__name__)
    t_start = time.perf_counter()
    t_count = 0.0
    n_regions_processed = 0
    n_exons_processed = 0

    # get the exon lists for each region
    region_exons = [region_fun(transcript) for region_fun in regions]
    region_lengths = [sum(x.end - x.start for x in r) for r in region_exons]

    # keep only regions with positive length and preserve name/bin alignment
    filtered = [(r, l, n, b) for r, l, n, b in zip(region_exons, region_lengths, names, bins) if l > 0]

    # prepare collectors
    region_binned_counts = []
    valid_names = []

    if aggregate:
        from functools import reduce
        # For each region, fetch counts for all exons with a single call to
        # count_intervals (which itself fetches once) and then split per exon.
        for r, l, n, b in filtered:
            if len(r) == 0:
                continue

            n_regions_processed += 1
            exon_profiles = []

            # Fetch counts for the whole block covering the region exons in one go
            intervals = [(e.start, e.end) for e in r]
            intervals.sort(key=lambda x: x[0])  # sort by start position
            
            t0 = time.perf_counter()
            all_exon_counts = count_intervals(getter, intervals, contig=transcript[0].contig, strand=transcript[0].strand)
            t1 = time.perf_counter()
            t_count += (t1 - t0)

            # Vectorized binning: map genome positions to exon indices and bin indices
            exon_starts = np.array([e.start for e in r], dtype=float)
            exon_ends = np.array([e.end for e in r], dtype=float)
            exon_lens = exon_ends - exon_starts
            n_exons = len(r)

            if all_exon_counts.sum() == 0:
                # no counts; produce zero bins
                region_bins = np.zeros(b, dtype=float)

            else:
                positions = np.asarray(all_exon_counts.index.values, dtype=float)
                values = np.asarray(all_exon_counts.values, dtype=float)

                # map positions to exon indices using pandas IntervalIndex
                intervals_index = pd.IntervalIndex.from_arrays(exon_starts, exon_ends, closed='left')

                # vectorized lookup: returns -1 for positions not in any interval
                exon_idx = intervals_index.get_indexer(positions).astype(int)

                # ensure all positions map to an exon; if not, raise so we know
                # something unexpected occurred
                assert np.all(exon_idx >= 0), f"Position(s) not assigned to any exon for transcript {transcript[0].transcript_id}"

                # compute base offsets within exon and bin indices; account for
                # transcript strand so offsets are measured from the transcript
                # 5' end of the exon.
                if transcript[0].strand == "-":
                    offsets = exon_ends[exon_idx] - positions - 1
                else:
                    offsets = positions - exon_starts[exon_idx]

                lengths = exon_lens[exon_idx]
                # avoid division by zero
                lengths_nonzero = np.where(lengths == 0, 1, lengths)
                bin_idx = np.floor(offsets * b / lengths_nonzero).astype(int)
                bin_idx = np.clip(bin_idx, 0, b-1)

                # composite index to aggregate per-exon per-bin
                composite = exon_idx * b + bin_idx
                minlength = n_exons * b
                agg = np.bincount(composite, weights=values, minlength=minlength)
                agg = agg.reshape(n_exons, b)

                # sum across exons to get region bins
                region_bins = agg.sum(axis=0)

                # apply length normalization at the region level to match
                # non-aggregate behaviour (scale by b / region_length)
                if length_norm:
                    region_bins = region_bins * (float(b) / float(l))

            # Use the region-level binned counts (already length-normalized if requested)
            combined = pd.Series(region_bins, index=range(b))
            region_binned_counts.append(combined)
            valid_names.append(n)

    else:
        # non-aggregate: bin each region as a whole
        if len(filtered) == 0:
            empty_index = pd.MultiIndex(levels=[names, []], codes=[[], []], names=["region", "region_bin"])
            return pd.Series(dtype=float, index=empty_index)

        region_exons_f, region_lengths_f, valid_names, bins_f = zip(*filtered)

        region_counts = [count_transcript(t, getter) for t in region_exons_f]
        region_binned_counts = [bin_counts(c, l, b) for c, l, b in
                                zip(region_counts, region_lengths_f, bins_f)]

        if length_norm:
            region_binned_counts = [x * (float(b) / float(l)) for x, l, b in
                                    zip(region_binned_counts, region_lengths_f, bins_f)]

    # build final profile
    if len(region_binned_counts) == 0:
        empty_index = pd.MultiIndex(levels=[names, []], codes=[[], []], names=["region", "region_bin"])
        return pd.Series(dtype=float, index=empty_index)

    profile = pd.concat(region_binned_counts, keys=valid_names, names=["region", "region_bin"])

    t_end = time.perf_counter()
    logger.debug("transcript_region_meta: total_time=%.6fs, count_time=%.6fs, regions=%d, exons=%d", \
                 (t_end - t_start), t_count, n_regions_processed, n_exons_processed)

    return profile
