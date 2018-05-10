import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGATPipelines.Pipeline as P
import CGAT.FastaIterator as FastaIterator
import CGAT.Experiment as E
from CGATPipelines.Pipeline import cluster_runnable
from CGATPipelines.PipelineMapping import SequenceCollectionProcessor
import pandas
import os
import re
import pysam

# The PARAMS dictionary must be provided by the importing
# code

class SemiProcessedCollector(SequenceCollectionProcessor):
    '''This class takes local or remote, single ended files that have been
    semi-processed, i.e. they have been filtered, trimmed and
    demuxed. The UMI may not be in the expected place, but can be
    moved using pre-processed patterns in the [preprocessed] section
    of the ini

    '''

    def build(self, infiles, outfile):

        cmd_preprocess, mapfiles = self.preprocess(infiles, outfile)
        
        statement = [cmd_preprocess]
        assert len(mapfiles)==1
        infile = mapfiles[0][0]
        statement.append('''
                zcat %(infile)s |
                sed 's/ /_/g' |
                sed -E 's/%%(preprocess_in_pattern)s/%%(preprocess_out_pattern)s/g;n;n;n' |
                gzip > %(outfile)s;
           ''' % locals())

        return " checkpoint; ".join(statement)

    
def checkParams():

    if not len(PARAMS) > 0:
        raise ValueError(
                "Please set PARAMS dictionary in PipelineiCLIP module")


def getBarcodeCG(table, outfile):
    ''' Annotate barcode use statistics with %GC '''

    statement = " SELECT * FROM %(table)s" % locals()

    umi_stats = PUtils.fetch_DataFrame(statement)

    def _GC(x):
        return float(x.count("G") + x.count("G"))/len(x)

    barcode_gc = umi_stats.Barcode.apply(_GC)
    sample_gc = umi_stats.Sample.apply(_GC)
    umi_gc = umi_stats.UMI.apply(_GC)

    gc_stats= pandas.DataFrame({"Barcode":umi_stats.Barcode,
                                "barcode_gc": barcode_gc,
                                "sample_gc": sample_gc,
                                "umi_gc": umi_gc})

    gc_stats.to_csv(IOTools.openFile(outfile,"w"), sep="\t", index=False)


###################################################################
def callClusters(bamfile, gtffile, outfiles,
                 window_size=None,
                 pthresh=None):
    ''' Wrapper around find_reproducible_clusters.
    If no window_size is specified, it is taken from
    pipeline.ini'''

    checkParams()

    bedGraph, bed12 = outfiles
    logfile = P.snip(bed12, ".bed.gz")

    if window_size:
        options = "--window-size=%i" % window_size
    else:
        options = "--window-size=%s" % PARAMS["clusters_window_size"]

    if PARAMS["clusters_fdr"]:
        options += " --fdr"
    if PARAMS["clusters_grouping"]:
        options += " --grouping=%s" % PARAMS["clusters_grouping"]

    if pthresh:
        options += " -t %s" % str(pthresh)
    else:
        options += " -t %s" % PARAMS["clusters_pthresh"]

    
    job_options = "-l mem_free=10G"
    statement = '''python %(scriptsdir)s/gtf2gtf.py -L %(logfile)s.log
                           -I %(gtffile)s
                          --method=sort --sort-order=gene+transcript
                 | python %(scriptsdir)s/gtf2gtf.py -L %(logfile)s.log 
                          --method=set-transcript-to-gene
                 | python %(project_src)s/find_significant_bases.py
                   %(bamfile)s
                   %(options)s
                   --output-both=%(bed12)s
                  -L %(logfile)s.log
                | gzip -c > %(bedGraph)s '''

    P.run()


###################################################################
def callReproducibleClusters(infiles, outfile, min_overlap):
    '''Find clusters that appear in more than one replicate'''

    checkParams()

    merge_template = '''<( zcat %s
                          | sort -k1,1 -k2,2n
                          | cgat bed2bed
                          --method=merge
                          --merge-and-resolve-blocks
                          --merge-stranded
                           -L /dev/null
                            ) '''
    infiles = " ".join(
        [merge_template % infile
         for infile in infiles])
    logfile = P.snip(outfile, ".bed.gz")
    statement = ''' cat %(infiles)s
                  | sort -k1,1 -k2,2n -k3,3n
                  | cgat bed2bed
                          --method=merge
                          --merge-and-resolve-blocks
                          --merge-min-intervals=%(min_overlap)s
                          --merge-stranded
                           -L %(logfile)s.log
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
def removeInputOverlappingClusters(sample, control, outfile):
    '''Remove reproducible clusters that overlap with reproducible
    input clusters '''

    statement = ''' bedtools intersect -v -a %(sample)s -b %(control)s
                   | gzip > %(outfile)s '''
    P.run()


###################################################################
def clustersToBigBed(infile, outfile):
    '''Convert beds to bigbed '''

    checkParams()

    tmp = P.getTempFilename()
    genome_file = os.path.join(PARAMS["annotations_dir"],
                               PARAMS_ANNOTATIONS["interface_contigs_tsv"])
    statement = ''' zcat %(infile)s | sort -k1,1 -k2,2n 
                    | awk 'BEGIN{OFS="\\t"} $5=1' > %(tmp)s;
                    checkpoint;
                    bedToBigBed %(tmp)s %(genome_file)s %(outfile)s;
                    checkpoint;
                    rm %(tmp)s'''
    P.run()


###################################################################
def makeClustersUCSC(infiles, outfile, group, label):
    '''Compile UCSC track file for bigbeds from list of cluster files'''

    template = '''
         track %(track_name)s
         parent %(group)s %(visible)s
         shortLabel %(short_label)s
         longLabel %(long_label)s
         bigDataUrl %(big_data_url)s
         type bigBed 12'''

    outlines = []

    for infile in infiles:

        big_data_url = os.path.basename(infile)
        

        if "reproducible" in infile:
            visible = "on"
            track_name = group + "_" + re.match(
                ".+/(.+).reproducible.*bigBed", infile).groups()[0]
            long_label = "Clusters from %s appearing in at least %s replicates" \
                         % (track_name, PARAMS["clusters_min_reproducible"])
        else:
            visible = "off"
            track_name = group + "_" + re.match(
                ".+/(.+)\.bigBed", infile).groups()[0]
            long_label = "Clusters from track %s" % track_name

        short_label = "%s clusters" % track_name
        outlines.append(template % locals())

    composite_stanaz = '''

    track %(group)s
    superTrack on
    shortLabel %(label)s
    longLabel %(label)s
    '''
    outlines = [composite_stanaz % locals()] + outlines

    outlines = "\n".join(outlines)

    with IOTools.openFile(outfile, "w") as outf:
        outf.write(outlines+"\n\n")


