import numpy as np
import pandas as pd
import CGAT.GTF as GTF


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
      
        length = sum([x.end - x.start for x in transcript if x.feature == "exon"])

        if calculate_flanks:
            flanks = length * float(bins[0])/bins[1]
            length_lookup = dict(zip(*(regions, [flanks, 0, flanks])))

        counts = count_transcript(transcript, bam, flanks=flanks)

        if flanks > 0:
            length_lookup["exons"] = length
            binned_counts = counts.groupby(level=0).apply(
                lambda x: bin_counts(x, length_lookup[x.name], nbins_lookup[x.name]))
        else:
            binned_counts = bin_counts(counts, length, bins)

        binned_counts.name = transcript[0].transcript_id

        counts_collector.append(binned_counts)

    counts_matrix = pd.concat(counts_collector, axis=1)
    counts_matrix = counts_matrix.transpose()
    counts_matrix = (counts_matrix.fillna(0) + pseudo_count).div(counts_matrix.sum(axis=1), axis=0)
    summed_matrix = counts_matrix.sum()
    summed_matrix.name = "density"

    if output_matrix:
        return summed_matrix, counts_matrix
    else:
        return summed_matrix, None


