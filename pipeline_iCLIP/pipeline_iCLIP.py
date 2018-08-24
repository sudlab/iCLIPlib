###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
===========================
Pipeline iCLIP
===========================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

A pipeline template.

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use
CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

Input
-----

The inputs should be a fastq file. The pipeline expects these to be raw fastq
files. That is that they contain the UMIs and the barcodes still on the 5' end
of the reads. It is expected that these will not be demultiplexed, although
demultiplex fastqs will also work.

In addition to the fastq files, a table of barcodes and samples is required as
sample_table.tsv.

It has four columns:

The first contains the barcode including UMI bases, marked as Xs.
The second contains the barcode sequence without the UMI bases.
The third contains the sample name you'd like to use
The fourth contains the fastq files that contain reads from this sample

e.g.

NNNGGTTNN	GGTT	FlipIn-FLAG-R1	hiseq7,hiseq8,miseq

Means that the sample FlipIn-FLAG-R1 should have reads in the fastq files
hiseq7, hiseq8 and miseq, is marked by the barcode GGTT and is embeded in the
UMI as NNNGGTTNN.

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|CGAPipelines        |                   |Pipelining infrastructure, mapping pipeline     |
+--------------------+-------------------+------------------------------------------------+
|CGAT                | >=0.2.4           |Various                                         |
+--------------------+-------------------+------------------------------------------------+
|Bowtie              | >=1.1.2           |Filtering out PhiX reads                        |
+--------------------+-------------------+------------------------------------------------+
|FastQC              | >=0.11.2          |Quality Control of demuxed reads                |
+--------------------+-------------------+------------------------------------------------+
|STAR or HiSat       |                   |Spliced mapping of reads                        |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |Interval manipulation                           |
+--------------------+-------------------+------------------------------------------------+
|samtools            |                   |Read manipulation                               |
+--------------------+-------------------+------------------------------------------------+
|subread             |                   |FeatureCounts read quantifiacation              |
+--------------------+-------------------+------------------------------------------------+
|MEME                |>=4.9.1            |Motif finding                                   |
+--------------------+-------------------+------------------------------------------------+
|DREME               |>=4.9.1            |Motif finding                                   |
+--------------------+-------------------+------------------------------------------------+
|bedGraphToBigWig    |                   |Converstion of results to BigWig                |
+--------------------+-------------------+------------------------------------------------+
|reaper              | 13-100            |Used for demuxing and clipping reads            |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

As well as the report, clusters, as BED files are in the clusters.dir directory,
and traces as bigWig files are in the bigwig directory. Both of these are exported
as a UCSU genome browser track hub in the export directory. 

Example
=======

Example data and configuration is avaiable in example_data.tar.gz


Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
from ruffus.combinatorics import *

import sys, glob, gzip, os, itertools, re, math
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineMapping as PipelineMapping
import PipelineiCLIP
###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
import CGATPipelines.Pipeline as P
P.getParameters(
    [os.path.join(os.path.dirname(__file__), "configuration",
                  "pipeline.ini"),
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

PipelineiCLIP.PARAMS = PARAMS
PipelineiCLIP.PARAMS_ANNOTATIONS = PARAMS_ANNOTATIONS
PARAMS["project_src"] = os.path.join(os.path.dirname(__file__),
                           "..")

###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample3)
for line in IOTools.openFile("sample_table.tsv"):
    track = line.split("\t")[2]
    TRACKS.tracks.append(PipelineTracks.Sample3(filename=track))


###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect(PARAMS["database"])
    statement = '''ATTACH DATABASE '%s' as annotations''' \
                % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


###################################################################
###################################################################
# worker tasks
###################################################################
@transform(os.path.join(PARAMS["genome_dir"],
                        PARAMS["genome"]),
           formatter(),
           "contigs.tsv")
def getContigSizes(infile, outfile):

    from CGAT import IndexedFasta
    
    try:
        prefix = P.snip(infile, ".fasta")
    except ValueError:
        prefix = P.snip(infile, ".fa")
        
    fasta = IndexedFasta.IndexedFasta(prefix)
    outs = IOTools.openFile(outfile, "w")

    for contig, size in fasta.getContigSizes(with_synonyms=False).items():
        outs.write("%s\t%i\n" % (contig, size))

    outs.close()


    
@transform("*.fastq.gz", regex("(.+).fastq.gz"),
           add_inputs(os.path.join(PARAMS["bowtie_index_dir"],
                                   PARAMS["phix_genome"]+".fa")),
           r"\1.fastq.clean.gz")
def filterPhiX(infiles, outfile):
    ''' Use mapping to bowtie to remove any phiX mapping reads '''

    infile, reffile = infiles
    outfile = P.snip(outfile, ".gz")
    bam_out = P.snip(infile, ".fastq.gz") + ".phix.bam"

    job_threads = PARAMS["phix_bowtie_threads"]
    job_memory = PARAMS["phix_bowtie_memory"]
    options = PARAMS["phix_bowtie_options"] + " --un %s" % outfile
    genome = PARAMS["phix_genome"]
    bowtie_threads = PARAMS["phix_bowtie_threads"]

    m = PipelineMapping.Bowtie(executable=PARAMS["phix_bowtie_exe"],
                               strip_sequence=False,
                               remove_non_unique=False,
                               tool_options=options)

    statement = m.build((infile,), bam_out)
    statement += "checkpoint; gzip %(outfile)s"

    P.run()


@transform("sample_table.tsv", suffix(".tsv"), ".load")
def loadSampleInfo(infile, outfile):

    P.load(infile, outfile,
           options="--header-names=format,barcode,track,lanes -i barcode -i track")


###################################################################
@follows(mkdir("demux_fq"))
@transform(filterPhiX, regex("(.+).fastq.clean.gz"),
           r"demux_fq/\1.fastq.umi_trimmed.gz")
def extractUMI(infile, outfile):
    ''' Remove UMI from the start of each read and add to the read
    name to allow later deconvolving of PCR duplicates '''

    statement='''umi_tools extract
                        -I %(infile)s
                        --bc-pattern=%(reads_bc_pattern)s
                        -L %(outfile)s.log
                        -S %(outfile)s '''

    P.run()


###################################################################
@transform(extractUMI, suffix(".fastq.umi_trimmed.gz"),
           ".umi_stats.load")
def loadUMIStats(infile, outfile):
    ''' load stats on UMI usage from the extract_umi log into the
    database '''

    infile = infile + ".log"
    P.load(infile, outfile, "-i sample -i barcode -i UMI")


###################################################################
@transform(filterPhiX,
           regex("(.+).fastq.clean.gz"),
           add_inputs("sample_table.tsv"),
           r"\1_reaper_metadata.tsv")
def generateReaperMetaData(infile, outfile):
    '''Take the sample_table and use it to generate a metadata table
    for reaper in the correct format '''

    adaptor_5prime = PARAMS["reads_5prime_adapt"]
    adaptor_3prime = PARAMS["reads_3prime_adapt"]

    outlines = []
    lane = P.snip(infile[0], ".fastq.clean.gz")
    for line in IOTools.openFile(infile[1]):
        fields = line.split("\t")
        barcode = fields[1]
        lanes = fields[-1].strip().split(",")
        if lane in lanes:
            outlines.append([barcode, adaptor_3prime, adaptor_5prime, "-"])

    header = ["barcode", "3p-ad", "tabu", "5p-si"]
    IOTools.writeLines(outfile, outlines, header)


###################################################################
@follows(loadUMIStats, generateReaperMetaData)
@subdivide(extractUMI, regex(".+/(.+).fastq.umi_trimmed.gz"),
       add_inputs(r"\1_reaper_metadata.tsv", "sample_table.tsv"),
       r"demux_fq/*_\1.fastq.gz")
def demux_fastq(infiles, outfiles):
    '''Demultiplex each fastq file into a seperate file for each
    barcode/UMI combination'''

    infile, meta, samples = infiles
    track = re.match(".+/(.+).fastq.umi_trimmed.gz", infile).groups()[0]
    
    statement = '''reaper -geom 5p-bc
                          -meta %(meta)s
                          -i %(infile)s
                          --noqc
                          %(reads_reaper_options)s
                          -basename demux_fq/%(track)s_
                          -clean-length %(reads_min_length)s > %(track)s_reapear.log;
                   checkpoint;
                   rename _. _ demux_fq/*clean.gz;
                 '''

    for line in IOTools.openFile(samples):
        line = line.split("\t")
        bc, name, lanes = line[1:]
        name = name.strip()
        if PARAMS["reads_paired"]:
            ext = "fastq.1.gz"
        else:
            ext = "fastq.gz"
        if track in lanes.strip().split(","):
            statement += '''checkpoint;
                         mv demux_fq/%(track)s_%(bc)s.clean.gz
                            demux_fq/%(name)s_%(track)s.%(ext)s; ''' % locals()

    P.run()

    
###################################################################
@follows(mkdir("demux_fq"))
@transform("*.remote", formatter(),
           "demux_fq/{basename[0]}.fastq.gz")
def get_remote_reads(infile, outfile):
    '''Get remote reads for reprocessing'''

    m = PipelineiCLIP.SemiProcessedCollector()
    statement = m.build((infile,), outfile)
    P.run()



###################################################################
@follows(mkdir("fastqc"))
@transform([demux_fastq,
            get_remote_reads],
           regex(".+/(.+).fastq.gz"),
           r"fastqc/\1.fastqc")
def qcDemuxedReads(infile, outfile):
    ''' Run fastqc on the post demuxing and trimmed reads'''

    m = PipelineMapping.FastQc(nogroup=False, outdir="fastqc")
    statement = m.build((infile,),outfile)
    exportdir = "fastqc"
    P.run()


###################################################################
@transform(qcDemuxedReads, regex("(.+)/(.+)\.fastqc"),
           inputs(r"\1/\2_fastqc/fastqc_data.txt"),
           r"\1/\2_length_distribution.tsv")
def getLengthDistribution(infile, outfile):
    ''' Parse length distribution out of the fastqc results '''

    statement = '''
       sed -e '/>>Sequence Length Distribution/,/>>END_MODULE/!d' %(infile)s
     | grep -P '^[0-9]+'
     | sed '1istart\\tend\\tcount'
     | sed 's/-/\\t/' > %(outfile)s '''

    P.run()


###################################################################
@merge(getLengthDistribution, "read_length_distribution.load")
def loadLengthDistribution(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_length_distribution.tsv",
                         options="-i start -i end")



###################################################################
@follows(demux_fastq,qcDemuxedReads, loadUMIStats, loadSampleInfo,
         loadLengthDistribution)
def PrepareReads():
    pass


###################################################################
# Mapping
###################################################################
def mapping_files():

    infiles = glob.glob("%s/*.fastq*gz" % PARAMS["input"])
    infiles = [infile for infile in infiles if ".fastq.umi_trimmed.gz" not in infile]
    outfiles = set(
        ["mapping.dir/%(mapper)s.dir/%(track)s.%(mapper)s.bam"
         % {"mapper": PARAMS["mappers"],
            "track": re.match("%s/(.+).fastq.*gz"
                              % PARAMS["input"], infile).groups()[0]}
         for infile in infiles
         if "umi_trimmed" not in infile])

    yield (infiles, outfiles)


###################################################################
@follows(mkdir("mapping.dir"), demux_fastq, get_remote_reads)
@transform([demux_fastq,
            get_remote_reads],
           regex(".+/(.+).fastq.gz"),
           r"mapping.dir/\1.bam")
def run_mapping(infile, outfile):
    ''' Map reads with the specified read mapper '''

    if PARAMS["mapper"] == "star":
        job_threads=PARAMS["star_threads"]
        job_memory=PARAMS["star_memory"]
        star_mapping_genome = PARAMS["star_genome"] or PARAMS["genome"]
        m = PipelineMapping.STAR(
            executable=P.substituteParameters(**locals())["star_executable"],
            strip_sequence=0)
        
    elif PARAMS["mapper"] == "bowtie":
         job_threads = PARAMS["bowtie_threads"]
         job_memory = PARAMS["bowtie_memory"]

         m = PipelineMapping.Bowtie(
             executable="bowtie",
             tool_options=PARAMS["bowtie_options"],
             strip_sequence=0)

         genome = PARAMS["bowtie_genome"]
         reffile = os.path.join(PARAMS["bowtie_index_dir"],
                                PARAMS["bowtie_genome"] + ".fa")

    statement = m.build((infile,), outfile)

    P.run()


###################################################################
@collate(run_mapping,
         regex("(.+)/([^_]+\-.+\-[^_]+)_(.+)\.bam"),
         r"\1/merged_\2.bam")
def mergeBAMFiles(infiles, outfile):
    '''Merge reads from the same library, run on different lanes '''

    if len(infiles) == 1:
        P.clone(infiles[0], outfile)
        P.clone(infiles[0]+".bai", outfile+".bai")
        return

    infiles = " ".join(infiles)
    statement = ''' samtools merge %(outfile)s %(infiles)s >& %(outfile)s.log'''
    P.run()


###################################################################
@transform(mergeBAMFiles, suffix(".bam"), ".bam.bai")
def indexMergedBAMs(infile, outfile):
    ''' Index the newly merged BAM files '''

    statement = ''' samtools index %(infile)s '''
    P.run()


###################################################################

@transform(PARAMS["annotations_geneset"],
           regex("(?:.+/)?([^/]+).gtf.gz"),
           add_inputs(getContigSizes),
           r"\1.context.bed.gz")
def generateContextBed(infiles, outfile):
    ''' Generate full length primary transcript annotations to count
    mapping contexts '''

    infile, genome = infiles
    
    statement = ''' zcat %(infile)s
                  | awk '$3=="exon"'
                  | cgat gtf2gtf
                    --method=exons2introns
                    
                     -L %(outfile)s.log
                  | awk 'BEGIN{FS="\\t";OFS="\\t"} {$2="intron"; print}'
                  | gzip > %(outfile)s.tmp.gtf.gz;

                  checkpoint;

                  zcat %(infile)s %(outfile)s.tmp.gtf.gz
                  | awk '$3=="exon" || $3=="intron"'
                  | cgat gff2bed
                    --set-name=source
                     -L %(outfile)s.log
                  | bedtools sort -i - -faidx %(genome)s
                    > %(outfile)s.tmp.bed;
                    
                    checkpoint;
                    
                    bedtools complement -i %(outfile)s.tmp.bed -g %(genome)s
                  | awk 'BEGIN{OFS="\\t"} {print ($0 "\\tnone\\t0\\t.")}'
                    > %(outfile)s.tmp2.bed;

                    checkpoint;

                    cat %(outfile)s.tmp.bed %(outfile)s.tmp2.bed
                  | sort -k1,1 -k2,2n
                  | gzip > %(outfile)s;
                  
                    checkpoint;
                  
                    rm %(outfile)s.tmp.bed %(outfile)s.tmp2.bed %(outfile)s.tmp.gtf.gz '''

    P.run()


###################################################################
@transform(generateContextBed, suffix(".context.bed.gz"),
           ".context_interval_stats.tsv.gz")
def getContextIntervalStats(infile, outfile):
    ''' Generate length stastics on context file '''

    statement = ''' cgat bed2stats
                            --aggregate-by=name
                            -I %(infile)s
                    | gzip > %(outfile)s '''

    P.run()



###################################################################
@transform(getContextIntervalStats, suffix(".tsv.gz"),
           ".load")
def loadContextIntervalStats(infile, outfile):

    P.load(infile, outfile)


###################################################################
@follows(loadContextIntervalStats )
def mapping():
    pass


###################################################################
# Deduping, Counting, etc
###################################################################
@follows(mkdir("deduped.dir"), run_mapping)
@transform(indexMergedBAMs, regex("(.+)/merged_(.+).bam.bai"),
           inputs(r"\1/merged_\2.bam"),
           r"deduped.dir/\2.bam")
def dedup_alignments(infile, outfile):
    ''' Deduplicate reads, taking UMIs into account'''

    outfile = P.snip(outfile, ".bam")

    job_memory="7G"
    statement = ''' umi_tools dedup
                    %(dedup_options)s
                    -I %(infile)s
                    -S %(outfile)s.bam
                    -L %(outfile)s.log;
 
                    checkpoint;

                    samtools index %(outfile)s.bam ;''' 


    P.run()


@collate(dedup_alignments,
         regex("(.+-.+)-(.+).bam"),
         r"\1-union.bam")
def get_union_bams(infiles, outfile):
    '''Merge replicates (as defined by thrid slot in name) to createView
    "union" tracks '''

    if len(infiles) == 1:
        infile = infiles[0]
        infile = os.path.abspath(infile)
        statement = ''' ln -sf %(infile)s %(outfile)s;
                        checkpoint;
                        ln -sf %(infile)s.bai %(outfile)s.bai; '''
    else:
        infiles = " ".join(infiles)
        statement = ''' samtools merge %(outfile)s %(infiles)s;
                        checkpoint;
                        samtools index %(outfile)s'''

    P.run()

    
###################################################################
@transform([get_union_bams,
            dedup_alignments],
           suffix(".bam"),
           ".bed.gz")
def get_indexed_bed(infile, outfile):
    '''Convert BAMs of reads into beds of signal'''

    outfile = P.snip(outfile, ".gz")
    statement = '''python %(project_src)s/scripts/iCLIP2bigWig.py 
                     -I %(infile)s
                     %(outfile)s
                     --format=bed;

                     checkpoint;

                     bgzip -f %(outfile)s

                    checkpoint;

                     tabix -p bed %(outfile)s.gz'''

    P.run()

    
###################################################################
@transform([dedup_alignments,indexMergedBAMs,
            get_union_bams], 
           regex("(?:merged_)?(.+).bam(?:.bai)?"),
           r"\1.frag_length.tsv")
def getFragLengths(infile, outfile):
    ''' estimate fragment length distribution from read lengths'''

    intrack = re.match("(.+).bam(?:.bai)?", infile).groups()[0]

    curdir = os.path.dirname(__file__)
    statement = ''' python %(curdir)s/length_stats.py
                           -I %(intrack)s.bam
                           -S %(outfile)s
                           -L %(outfile)s.log
               '''

    P.run()


###################################################################
@collate(getFragLengths,
         regex("(mapping|deduped).dir/.+\.frag_length.tsv"),
         r"\1.frag_lengths.load")
def loadFragLengths(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+\-.+\-.+).frag_length.tsv",
                         options=" -i Length")


###################################################################
@transform(dedup_alignments,
           suffix(".bam"), ".bam_stats.tsv")
def dedupedBamStats(infile, outfile):
    ''' Calculate statistics on the dedeupped bams '''

    statement = '''cgat bam2stats
                         --force-output
                          < %(infile)s > %(outfile)s '''

    P.run()


###################################################################
@merge(dedupedBamStats, "deduped_bam_stats.load")
def loadDedupedBamStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).bam_stats.tsv")


###################################################################
@transform((dedup_alignments,
            get_union_bams),
           suffix(".bam"),
           ".nspliced.txt")
def getNspliced(infile, outfile):
    ''' Calculate the number of reads spliced by grepping the
    cigar string for Ns'''

    statement = ''' samtools view %(infile)s
                  | awk '$6 ~ /N/'
                  | wc -l > %(outfile)s '''
    P.run()


###################################################################
@merge(getNspliced, "deduped_nspliced.load")
def loadNspliced(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).nspliced.txt",
                         cat="track",
                         has_titles=False,
                         header="track,nspliced",)


###################################################################
@transform(dedup_alignments, suffix(".bam"), ".umi_stats.tsv.gz")
def deduped_umi_stats(infile, outfile):
    ''' calculate histograms of umi frequencies '''

    curdir = os.path.dirname(__file__)
    statement = '''python %(curdir)s/umi_hist.py
                           -I %(infile)s
                           -L %(outfile)s.log
                  | gzip > %(outfile)s '''

    P.run()


###################################################################
@merge(deduped_umi_stats, "dedup_umi_stats.load")
def loadDedupedUMIStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).umi_stats.tsv.gz",
                         cat="track",
                         options="-i track -i UMI")


###################################################################
@transform([indexMergedBAMs, dedup_alignments,
            get_union_bams],
           regex("(?:merged_)?(.+).bam(?:.bai)?"),
           add_inputs(generateContextBed),
           r"\1.reference_context.tsv")
def buildContextStats(infiles, outfile):
    ''' Find context from reads '''

    infile, reffile = infiles
    infile = re.match("(.+.bam)(?:.bai)?", infile).groups()[0]
    statement = ''' cgat bam_vs_bed
                   --min-overlap=0.5
                   --log=%(outfile)s.log
                   %(infile)s %(reffile)s
                > %(outfile)s '''

    job_options = "-l mem_free=4G"
    P.run()


###################################################################
@collate(buildContextStats,
         regex("(mapping|deduped|saturation).dir/(?:[^/]+.dir/)?(.+).tsv"),
         r"\1_context_stats.load")
def loadContextStats(infiles, outfile):

    if "saturation" in infiles[0]:
        regex_filename = ".+/(.+-.+-.+)\.([0-9]+\.[0-9]+).reference_context.tsv"
        cat = "track,subset"
    else:
        regex_filename = ".+/(.+).reference_context.tsv"
        cat = "track"

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=regex_filename,
                         cat=cat)


###################################################################
@follows(loadContextStats,
         loadDedupedBamStats,
         loadFragLengths,
         loadNspliced,
         loadDedupedUMIStats)
def MappingStats():
    pass

         
###################################################################
# Quality control and reproducibility
###################################################################
@follows(mkdir("reproducibility.dir"))
@collate(dedup_alignments, regex(".+/(.+\-.+)\-.+.bam"),
         r"reproducibility.dir/\1-agg.reproducibility.tsv.gz")
def calculateReproducibility(infiles, outfile):
    ''' Calculate cross-link reproducibility as defined by 
    Sugimoto et al, Genome Biology 2012 '''

    job_options = "-l mem_free=1G"
    infiles = " ".join(infiles)

    statement = '''python %(project_src)s/scripts/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                 | gzip > %(outfile)s '''
    P.run()


###################################################################
@follows(mkdir("reproducibility.dir"))
@merge(dedup_alignments,
       r"reproducibility.dir/agg-agg-agg.reproducibility.tsv.gz")
def reproducibilityAll(infiles, outfile):
    ''' Test weather sites from one file reproduce in other of different factors '''

    job_options = "-l mem_free=10G"
    infiles = " ".join(infiles)

    statement = '''python %(project_src)s/scripts/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                 | gzip > %(outfile)s '''
    P.run()


###################################################################
@merge(calculateReproducibility,
       r"reproducibility.dir/experiment_reproducibility.load")
def loadReproducibility(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile, cat="Experiment",
                         regex_filename=".+/(.+)-agg.reproducibility.tsv.gz",
                         options="-i Track -i fold -i level")


###################################################################
@transform(reproducibilityAll, regex("(.+)"),
           "reproducibility.dir/all_reproducibility.load")
def loadReproducibilityAll(infile, outfile):
    P.load(infile, outfile, "-i Track -i fold -i level")


###################################################################
@permutations(dedup_alignments, formatter(".+/(?P<TRACK>.+).bam"),
              2,
              "reproducibility.dir/{TRACK[0][0]}_vs_{TRACK[1][0]}.tsv.gz")
def computeDistances(infiles, outfile):
    ''' Compute the reproduciblity between each indevidual pair of samples
    this can then be readily converted to a distance measure'''

    track = infiles[0]
    infiles = " ".join(infiles)

    job_options="-l mem_free=2G"

    statement = '''python %(project_src)s/scripts/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                   -t %(track)s
                   -m 1
                 | gzip > %(outfile)s '''

    P.run()


###################################################################
@merge(computeDistances,
       "reproducibility.dir/reproducibility_distance.load")
def loadDistances(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+)_vs_(.+).tsv.gz",
                         cat="File1,File2",
                         options="-i Track, -i File2")


###################################################################
@follows(loadReproducibility,
         loadReproducibilityAll,
         loadDistances)
def reproducibility():
    pass


###################################################################
@follows(mkdir("counts.dir"))
@transform(get_indexed_bed,
           regex(".+/(.+).bed.gz"),
           add_inputs(PARAMS["annotations_geneset"]),
           r"counts.dir/\1.tsv.gz")
def countReadsOverGenes(infiles, outfile):
    ''' use feature counts to quantify the number of tags on each gene'''

    bamfile, annotations = infiles
    if PARAMS["reads_use_centre"]==1:
        use_centre = "--use-centre"
        
    statement = '''python %(project_src)s/scripts/count_clip_sites.py
                   -I %(annotations)s
                   --bed=%(bamfile)s
                   -f gene
                   %(use_centre)s
                   -S %(outfile)s '''
    P.run()


###################################################################
@merge(countReadsOverGenes,
       "counts.dir/track_counts.tsv.gz")
def loadCounts(infiles, outfile):
    '''Merge feature counts data into one table'''

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.+).tsv.gz",
                         options="-i track,gene_id")


###################################################################
# Analysis
###################################################################
@follows(mkdir("gene_profiles.dir"))
@transform(get_indexed_bed,
           regex(".+/(.+).bed.gz"),
           add_inputs(PARAMS["annotations_geneset"]),
           r"gene_profiles.dir/\1.exons.tsv.gz")
def calculateGeneExonProfiles(infiles, outfile):
    ''' Calculate metagene profiles over protein coding genes
    for each sample'''

    infile, reffile = infiles

    if PARAMS["reads_use_centre"]==1:
        use_centre = "--use-centre"
    else:
        use_centre = ""
        
    statement= '''python %(project_src)s/scripts/iCLIP_bam2geneprofile.py
                  -I %(reffile)s
                  -b 100
                  --flank-bins=50
                  --scale-flanks
                  %(use_centre)s
                  -S %(outfile)s
                  -L %(outfile)s.log
                  --bed=%(infile)s'''
    
    P.run()

    
###################################################################
###################################################################
@transform(get_indexed_bed,
           regex(".+/(.+).bed.gz"),
           add_inputs(PARAMS["annotations_geneset"]),
           r"gene_profiles.dir/\1.introns.tsv.gz")
def calculateGeneIntronProfiles(infiles, outfile):
    '''Get a metagene profile over concatenated introns'''

    bamfile, gtffile = infiles

    if PARAMS["reads_use_centre"]==1:
        use_centre = "--use-centre"
    else:
        use_centre = ""

        
    statement = '''cgat gtf2gtf -I %(gtffile)s
                                --method=exons2introns
                                -L %(outfile)s.log
                 | awk -F'\\t' 'BEGIN{OFS=FS} {$3="exon"; print}'
                 | python %(project_src)s/scripts/iCLIP_bam2geneprofile.py
                   -f 0
                   --bed=%(bamfile)s
                   %(use_centre)s
                   --exon-bins=100
                   -S %(outfile)s
                   -L %(outfile)s.log; '''

    P.run()

    
###################################################################
@merge((calculateGeneExonProfiles,
        calculateGeneIntronProfiles),
       "gene_profiles.load")
def loadGeneProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename='.+/(.+)\-(.+)\-(.+)\.(.+).tsv.gz',
                          cat = "factor,condition,rep,interval",
                          options = "-i factor -i condition -i rep -i interval")


###################################################################
@follows(loadGeneProfiles)
def profiles():
    pass


###################################################################
# Calling significant clusters
###################################################################
@follows(mkdir("clusters.dir"))
@transform(get_indexed_bed,
           regex(".+/(.+).bed.gz"),
           add_inputs(PARAMS["annotations_geneset"]),
           r"clusters.dir/\1.sig_bases.bed.gz")
def callSignificantBases(infiles, outfile):
    '''Call bases as significant based on mapping depth in window
    around base'''

    bam, gtf = infiles
    job_threads = PARAMS["clusters_threads"]
    job_memory="10G"
    if PARAMS["reads_use_centre"]==1:
        centre = "--centre"
    else:
        centre = ""
        
    statement = ''' python %(project_src)s/scripts/significant_bases_by_randomisation.py
                    -I %(gtf)s
                    --bed=%(bam)s
                    --spread=%(clusters_window_size)s
                    -p %(job_threads)s
                    %(centre)s
                    -L %(outfile)s.log
                    --threshold=%(clusters_pthresh)s
                | bgzip > %(outfile)s;

                checkpoint;

                tabix -p bed %(outfile)s
                '''

    P.run()

    
###################################################################
@transform(callSignificantBases,
           suffix("sig_bases.bed.gz"),
           "clusters.bed.gz")
def callSignificantClusters(infile, outfile):
    '''Join significant bases that are within in the cluster
    distance of each other, using the max score as the score for
    the cluster'''

    statement='''bedtools merge 
                  -i %(infile)s
                  -d %(clusters_window_size)s
                  -s
                  -c 4 
                  -o max
                | gzip > %(outfile)s'''

    P.run()
    
###################################################################
@collate(callSignificantClusters,
         regex("(.+/.+\-.+)\-(.+)\.bed.gz"),
         r"\1.reproducible.bed.gz")
def callReproducibleClusters(infiles, outfile):
    '''Find clusters that appear in more than one replicate'''

    PipelineiCLIP.callReproducibleClusters(infiles, outfile,
                                           PARAMS["clusters_min_reproducible"])


###################################################################
@transform(callSignificantBases, suffix(".bed.gz"), ".count_bases")
def countCrosslinkedBases(infile, outfile):
    ''' Count number of crosslinked bases within gene models'''

    statement = ''' zcat %(infile)s |
                    wc -l > %(outfile)s '''

    P.run()


###################################################################
@merge(countCrosslinkedBases, "cross_linked_bases.load")
def loadCrosslinkedBasesCount(infiles, outfile):
    
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).count_bases",
                         header="track,count",
                         cat="track",
                         has_titles=False)


###################################################################
@transform([callSignificantClusters, callReproducibleClusters],
           suffix(".clusters.bed.gz"),
           r".cluster_count")
def countClusters(infile, outfile):
    '''Count the number of significant clusters'''

    statement = '''zcat %(infile)s | wc -l > %(outfile)s'''
    P.run()


###################################################################
@merge(countClusters, "cluster_counts.load")
def loadClusterCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).(R[0-9]+|union|reproducible).cluster_count",
                         header="sample,replicate,count",
                         cat="sample,replicate",
                         has_titles=False)

    
###################################################################
@transform(callSignificantBases,
           suffix(".bed.gz"),
           add_inputs(generateContextBed),
           ".context_stats.tsv.gz")
def getSigBaseContextStats(infiles, outfile):
    '''Generate context stats for significant bases'''

    bases, context = infiles
    tmp = P.getTempFilename()

    statement =  '''  zcat %(bases)s | sort -k1,1 -k2,2n > %(tmp)s.bed;
                     checkpoint;
                     cgat bam_vs_bed
                           -a %(tmp)s.bed
                           -b %(context)s -S  %(outfile)s;
                     checkpoint;
                     rm %(tmp)s.bed'''

    P.run()

    
###################################################################
@transform([callSignificantClusters, callReproducibleClusters],
           suffix(".clusters.bed.gz"),
           add_inputs(generateContextBed),
           ".clusters.context_stats.tsv.gz")
def getClusterContextStats(infiles, outfile):
    '''Generate context stats for called clusters'''

    clusters, context = infiles
    tmp = P.getTempFilename()

    statement = '''  zcat %(clusters)s | sort -k1,1 -k2,2n > %(tmp)s.bed;
                     checkpoint;
                     cgat bam_vs_bed
                           -a %(tmp)s.bed
                           -b %(context)s -S  %(outfile)s;
                     checkpoint;
                     rm %(tmp)s.bed'''

    P.run()


###################################################################
@collate([getClusterContextStats, getSigBaseContextStats],
         regex(".+/(.+).(sig_bases|clusters).context_stats.tsv.gz"),
         r"clusters.dir/\2_context_stats.load")
def loadClusterContextStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=
                         "clusters.dir/(.+)(?:sig_bases|clusters).context_stats.tsv.gz")


###################################################################
@follows(callSignificantClusters,
         loadCrosslinkedBasesCount,
         loadClusterCounts,
         loadClusterContextStats)
def clusters():
    pass


###################################################################
# Motifs
###################################################################
@follows(mkdir("kmers.dir"))
@subdivide([dedup_alignments,
            get_union_bams,
            callSignificantBases],
           regex(".+/(.+)(.bam|.bed.gz)"),
           add_inputs(PARAMS["annotations_geneset"],
                      os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"])),
           [r"kmers.dir/\1.%smers.tsv.gz" % kmer
            for kmer in str(PARAMS["kmer_lengths"]).split(",")])
def get_kmer_enrichments(infiles, outfiles):
    ''' Look for enrichment of kmers in crosslinked bases. Will
    examine both significant bases and all bases and look for 
    kmers of a length up that specified in ini'''
    
    bamfile, geneset, genome = infiles

    if bamfile.endswith(".bam"):
        sites = "-b " + bamfile
    elif bamfile.endswith(".bed.gz"):
        sites = "--bed=" + bamfile

    if PARAMS["reads_use_centre"]==1:
        centre="--use-centre"
    else:
        centre=""

    job_threads=PARAMS["kmer_threads"]
        
    statement_template = '''python %%(project_src)s/scripts/iCLIP_kmer_enrichment.py
                                    %(sites)s
                                    %(centre)s
                                   -f %(genome)s
                                   -k %(kmer)s
                                   -s %%(clusters_window_size)s
                                   -p %%(kmer_threads)s
                                   -I %(geneset)s
                                   -S %(outfile)s
                                   -L %(outfile)s.log'''

    statements = []
    
    for outfile, kmer in zip(outfiles, str(PARAMS["kmer_lengths"]).split(",")):
        statements.append(statement_template % locals())

    P.run()


###################################################################
@collate(get_kmer_enrichments,
         regex("(.+-.+)-(.+)\.[0-9]+mers.tsv.gz"),
         "kmers.dir/all_bases.kmers.load")
def load_dedup_kmers(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="kmers.dir/(.+)-(.+)-(.+)\.([0-9]+mers).tsv.gz",
                         cat="factor,tag,replicate,k",
                         options="-i factor -i tag -i replicate -i kmer")


###################################################################   
@collate(get_kmer_enrichments,
         regex(".+/[^\.]+\.sig_bases.[0-9]+mers.tsv.gz"),
         "kmers.dir/sig_bases.kmers.load")
def load_sig_base_kmers(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="kmers.dir/(.+)-(.+)-(.+)\.sig_bases.([0-9]+mers).tsv.gz",
                         cat="factor,tag,replicate,k",
                         options="-i factor -i tag -i replicate -i kmer")

    
###################################################################
@follows(load_dedup_kmers,
         load_sig_base_kmers)
def kmers():
    pass


###################################################################
# Data Export
###################################################################
@follows(mkdir("bigWig"))
@subdivide([dedup_alignments,
            get_union_bams], 
           regex(".+/(.+).bam"),
          [r"bigWig/\1_plus.bw",
           r"bigWig/\1_minus.bw"])
def generateBigWigs(infile, outfiles):
    '''Generate plus and minus strand bigWigs from BAM files '''

    out_pattern = P.snip(outfiles[0], "_plus.bw")
    statement = '''python %(project_src)s/scripts/iCLIP2bigWig.py
                          -I %(infile)s
                          -L %(out_pattern)s.log
                          %(out_pattern)s '''

    P.run()
    
    
###################################################################
@follows(mkdir("export/hg19"))
@transform(generateBigWigs,
           regex("bigWig/(.+)"),
           r"export/hg19/\1")
def linkBigWig(infile, outfile):
    '''Link bigwig files to export directory'''
    
    try:
        os.symlink(os.path.abspath(infile), os.path.abspath(outfile))
    except OSError:
        os.unlink(outfile)
        os.symlink(os.path.abspath(infile), os.path.abspath(outfile))


###################################################################
@merge(linkBigWig, "export/hg19/tagwig_trackDb.txt")
def generateBigWigUCSCFile(infiles, outfile):
    '''Generate track configuration for exporting wig files '''


    track_template = '''
          track tagwig_%(track)s_%(strand)s
          parent tagwig_%(track)s
          bigDataUrl %(track_data_URL)s
          shortLabel %(short_label)s
          longLabel %(long_label)s
          color %(color)s
          type bigWig %(negate)s'''

    overlap_template = '''
       track tagwig_%(track)s
       parent clipwig
       shortLabel %(short_label)s
       longLabel %(long_label)s
       autoScale on
       visibility full
       container multiWig
       type bigWig
       aggregate solidOverlay
       maxHeightPixels 16:16:32
       alwaysZero on'''

    stanzas = {}
    for infile in infiles:
        track, strand = re.match(
            ".+/(.+-.+)_(plus|minus).bw", infile).groups()

        negate = ""

        if strand == "minus":
            negate = '''
          negateValues on'''
            color="255,0,0"
        else:
            color="0,0,255"
            
        track_data_URL = os.path.basename(infile)
        short_label = track + "iCLIP tags"
        long_label = "iCLIP tags from track %s" \
                     % track
        
        if track not in stanzas:
            stanzas[track] = overlap_template % locals()

        stanzas[track] += "\n" + track_template % locals()

    composite_stanaz = '''
    track clipwig
    shortLabel iCLIP tags
    longLabel Raw iCLIP tags
    superTrack on
    alwaysZero on
    maxHeightPixels 16:16:32'''

    output = "\n".join([composite_stanaz] + list(stanzas.values()))

    with IOTools.openFile(outfile, "w") as outf:
        outf.write(output)


###################################################################
@follows(mkdir("export/hg19"))
@transform([callSignificantClusters, callReproducibleClusters],
           regex("clusters.dir/(.+).bed.gz"),
           add_inputs(getContigSizes),
           r"export/hg19/\1.bigBed")
def exportClusters(infiles, outfile):
    ''' Add a track line to cluster files and export to export dir '''

    infile, genome = infiles
    PipelineiCLIP.clustersToBigBed(infile, genome, outfile)


###################################################################
@merge(exportClusters, "export/hg19/clusters_trackDb.txt")
def generateClustersUCSC(infiles, outfile):

    PipelineiCLIP.makeClustersUCSC(infiles, outfile, "pipelineClusters",
                                   "Clusters from iCLIP pipeline")


###################################################################
@merge([generateClustersUCSC, generateBigWigUCSCFile],
       "export/hg19/trackDb.txt")
def mergeTrackDbs(infiles, outfile):

    to_cluster = False
    infiles = " ".join(infiles)
    statement = "cat %(infiles)s > %(outfile)s"
    P.run()


###################################################################
@follows(mkdir("export"))
@originate(["export/hub.txt",
            "export/genomes.txt"])
def makeHubFiles(outfiles):

    hub_file = '''
    hub iCLIPPipeline%(version)s
    shortLabel CGAT iCLIP Pipelines
    longLabel All browser tracks CGAT iCLIP pipeline run
    genomesFile genomes.txt
    email i.sudbery@sheffield.ac.uk''' % PARAMS

    with IOTools.openFile("export/hub.txt", "w") as outf:
        outf.write(hub_file)

    genomes_file = '''
    genome hg19
    trackDb hg19/trackDb.txt'''

    with IOTools.openFile("export/genomes.txt", "w") as outf:
        outf.write(genomes_file)


###################################################################
@follows(mergeTrackDbs, makeHubFiles)
def export():
    pass


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(PrepareReads, mapping, MappingStats, reproducibility,
         profiles, clusters, kmers, export)
def full():
    pass


@follows( mkdir( "report" ))
def build_report():
    '''build report from scratch.'''

    try:
        os.symlink(os.path.abspath("conf.py"),
                   os.path.join(
                       os.path.abspath("mapping.dir"), "conf.py"))
    except OSError as e:
        E.warning(str(e))

    E.info("Running mapping report build from scratch")
#    statement = '''cd mapping.dir;
#                   python %(scripts_dir)s/CGATPipelines/pipeline_mapping.py
#                   -v5 -p1 make build_report '''
#    P.run()
    E.info("starting report build process from scratch")
    P.run_report( clean = True )


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("Updating Mapping Reported")
    statement = '''cd mapping.dir;
                   python %(pipelinedir)s/pipeline_mapping.py
                   -v5 -p1 make update_report '''
    E.info("updating report")
    P.run_report(clean=False)


@follows( update_report )
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
