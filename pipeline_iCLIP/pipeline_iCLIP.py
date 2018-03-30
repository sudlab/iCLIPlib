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
    ["%s.ini" % P.snip(__file__, ".py"),
     "pipeline.ini"])

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

PipelineiCLIP.PARAMS = PARAMS
PipelineiCLIP.PARAMS_ANNOTATIONS = PARAMS_ANNOTATIONS

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
           "umi_stats.load")
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
                          -i <( zcat %(infile)s | sed 's/ /_/g')
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
@follows(mkdir("fastqc"))
@transform(demux_fastq, regex(".+/(.+).fastq(.*)\.gz"),
           r"fastqc/\1\2.fastqc")
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
@follows(mkdir("mapping.dir"), demux_fastq)
@transform(demux_fastq,
           regex(".+/(.+).fastq.gz"),
           "mapping.dir/\1.bam")
def run_mapping(infile, outfiles):
    ''' Map reads with the specified read mapper '''


    if PARAMS["mapper"]=="star":
        job_threads=PARAMS["star_threads"]
        job_memory=PARAMS["star_memory"]
        star_mapping_genome = PARAMS["star_genome"] or PARAMS["genome"]
        m = PipelineMapping.STAR(
            executable=P.substituteParameters(**locals())["star_executable"],
            strip_sequence=0)
        
    elif PARAMS["mapper"]=="bowtie":
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

@transform(os.path.join(PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
           regex(".+/(.+).gtf.gz"),
           r"\1.context.bed.gz")
def generateContextBed(infile, outfile):
    ''' Generate full length primary transcript annotations to count
    mapping contexts '''

    genome = os.path.join(PARAMS["annotations_dir"],
                          PARAMS_ANNOTATIONS["interface_contigs"])
    statement = ''' zcat %(infile)s
                  | awk '$3=="exon"'
                  | python %(scripts_dir)s/gtf2gtf.py
                    --method=exons2introns
                    
                     -L %(outfile)s.log
                  | awk 'BEGIN{FS="\\t";OFS="\\t"} {$2="intron"; print}'
                  | gzip > %(outfile)s.tmp.gtf.gz;

                  checkpoint;

                  zcat %(infile)s %(outfile)s.tmp.gtf.gz
                  | awk '$3=="exon" || $3=="intron"'
                  | python %(scripts_dir)s/gff2bed.py
                    --set-name=source
                     -L %(outfile)s.log
                  | sort -k1,1 -k2,2n
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

    statement = ''' python %(scripts_dir)s/bed2stats.py
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
@follows(mapping_qc,loadContextIntervalStats )
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
def get_union_bams(infile, outfile):
    '''Merge replicates (as defined by thrid slot in name) to createView
    "union" tracks '''

    if len(infiles) == 1:
        infile = infiles[0]
        statement = ''' ln -sf %(infile)s %(outfile)s;
                        checkpoint;
                        ln -sf %(infile)s.bai %(outfile)s.bai; '''
    else:
        infiles = " ".join(infiles)
        statement = ''' samtools merge %(infiles)s > %(outfile)s;
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
    statement = '''python %(project_src)s/iCLIP2bigWig.py 
                     -I %(infile)s
                     %(outfile)s
                     --format=bed;

                     checkpoint;

                     bgzip %(outfile)s

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

    statement = ''' python %(project_src)s/length_stats.py
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

    statement = '''python %(scripts_dir)s/bam2stats.py
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

    statement = '''python %(project_src)s/umi_hist.py
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
    statement = ''' python %(scripts_dir)s/bam_vs_bed.py
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
         loadDedupedUMIStats,
         loadSplicingIndex)
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

    statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
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

    statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                 | gzip > %(outfile)s '''
    P.run()


###################################################################
@follows(mkdir("reproducibility.dir"))
@transform(dedup_alignments,
           regex(".+/(.+).bam"),
           add_inputs("deduped.dir/%s*.bam" % PARAMS["experiment_input"]),
           r"reproducibility.dir/\1_vs_control.reproducibility.tsv.gz")
def reproducibilityVsControl(infiles, outfile):
    '''Test what fraction of the locations in each bam also appear
    in the control files'''

    track = infiles[0]
    if track in infiles[1:]:
        P.touch(outfile)
    else:
        job_options = "-l mem_free=1G"

        infiles = " ".join(infiles)

        statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
                   %(infiles)s
                   -L %(outfile)s.log
                   -t %(track)s
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
@merge(reproducibilityVsControl,
       "reproducibility.dir/reproducibility_vs_control.load")
def loadReproducibilityVsControl(infiles, outfile):

    infiles = [infile for infile in infiles 
               if not PARAMS["experiment_input"] in infile]

    P.concatenateAndLoad(infiles, outfile, cat="Experiment",
                         regex_filename=".+/(.+\-.+)\-.+_vs_control.reproducibility.tsv.gz",
                         options = "-i File -i fold -i level")


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

    statement = '''python %(project_src)s/iCLIPlib/calculateiCLIPReproducibility.py
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
         loadReproducibilityVsControl,
         loadDistances)
def reproducibility():
    pass


###################################################################
@follows(mkdir("counts.dir"))
@transform(get_indexed_bed,
           regex(".+/(.+).bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"counts.dir/\1.tsv.gz")
def countReadsOverGenes(infiles, outfile):
    ''' use feature counts to quantify the number of tags on each gene'''

    bamfile, annotations = infiles
    if PARAMS["use_centre"]==1:
        use_centre = "--use-centre"
        
    statement = '''python %(SRCDIR)s/iCLIPlib/scripts/count_clip_sites.py
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
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"gene_profiles.dir/\1.exons.tsv")
def calculateGeneExonProfiles(infiles, outfile):
    ''' Calculate metagene profiles over protein coding genes
    for each sample'''

    infile, reffile = infiles

    if PARAMS["use_centre"]==1:
        use_centre = "--use-centre"
        
    statement= '''python %(SRCDIR)s/iCLIPlib/iCLIP_bam2geneprofile.py
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
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"gene_profiles.dir/\1.introns.tsv")
def calculateGeneIntronProfiles(infiles, outfile):
    '''Get a metagene profile over concatenated introns'''

    bamfile, gtffile = infiles

    if PARAMS["use_centre"]==1:
        use_centre = "--use-centre"
        
    statement = '''cgat gtf2gtf -I %(gtffile)s
                                --method=exons2introns
                                -L %(outfile)s.log
                 | awk -F'\\t' 'BEGIN{OFS=FS} {$3="exon"; print}'
                 | python %(project_src)s/iCLIPlib/scripts/iCLIP_bam2geneprofile.py
                   -f 0
                   --bed=%(bamfile)s
                   %(use_centre)s
                   --exon-bins=100
                   -S %(outfile)s
                   -L %(outfile)s; '''

    P.run()

    
###################################################################
@merge((calculateGeneExonProfiles,
        calculateGeneIntronProfiles),
       "gene_profiles.load")
def loadGeneProfiles(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                          regex_filename='.+/(.+)\-(.+)\-(.+)\.(.+).tsv.gz',
                          cat = "factor,condition,rep,",
                          options = "-i factor -i condition -i rep")


###################################################################
@follows(calculateExonTSSProfiles,
         calculateGeneProfiles,
         calculateExonProfiles,
         loadExonProfiles,
         loadGeneProfiles)
def profiles():
    pass


###################################################################
# Calling significant clusters
###################################################################
@follows(mkdir("clusters.dir"), mapping_qc)
@transform(get_indexed_bed,
           regex(".+/(.+).bed.gz"),
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"])),
           r"clusters.dir/\1.sig_bases.bed.gz")
def callSignificantBases(infiles, outfiles):
    '''Call bases as significant based on mapping depth in window
    around base'''

    bam, gtf = infiles
    job_threads = PARAMS["clusters_threads"]

    if PARAMS["use_centre"]==1:
        centre = "--centre"
    else:
        centre = ""
        
    statement = ''' python %(project_src)s/scripts/significant_bases_by_randomisation.py
                    -I %(gtf)s
                    -bed=%(bam)s
                    --spread=%(clusters_window_size)s
                    -p %(job_threads)s
                    %(centre)s
                    -L %(outfile)s.log
                    --threashold=%(clusters_pthresh)s
                | bgzip > %(outfile)s;

                checkpoint;

                tabix -p bed %(outfile)s
                '''

    P.run()

    
###################################################################
@transform(callSignificantBases,
           suffix("sig_bases.bed.gz"),
           "clusters.bed.gz")
def callSignificantClusters(infile, outifle):
    '''Join significant bases that are within in the cluster
    distance of each other, using the max score as the score for
    the cluster'''

    statement='''bedtools merge 
                  -i %(infile)s
                  -d %(cluster_window_size)s
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
           suffix(".bed.gz"),
           r".cluster_count")
def countClusters(infile, outfile):
    '''Count the number of significant clusters'''

    statement = '''zcat %(infile)s | wc -l > %(outfile)s'''
    P.run()


###################################################################
@merge(countClusters, "cluster_counts.load")
def loadClusterCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).(R[0-9]+|reproducible).cluster_count",
                         header="sample,replicate,count",
                         cat="sample,replicate",
                         has_titles=False)

    
###################################################################
@transform(callSignificantBases,
           suffix(".bg.gz"),
           add_inputs(generateContextBed),
           ".context_stats.tsv.gz")
def getSigBaseContextStats(infiles, outfile):
    '''Generate context stats for significant bases'''

    bases, context = infiles
    tmp = P.getTempFilename

    statement =  '''  zcat %(clusters)s | sort -k1,1 -k2,2n > %(tmp)s.bed;
                     checkpoint;
                     cgat bam_vs_bed
                           -a %(tmp)s.bed
                           -b %(context)s -S  %(outfile)s;
                     checkpoint;
                     rm %(tmp)s.bed'''

    P.run()

    
###################################################################
@transform([callSignificantClusters, callReproducibleClusters],
           suffix(".bed.gz"),
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
         regex(".+/(.+).(bases|clusteres).context_stats.tsv.gz"),
         r"clusters.dir/\2_context_stats.load")
def loadClusterContextStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=
                         "clusters.dir/(.+).context_stats.tsv.gz")


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
@follows(mkdir("kmers"))
@transform([dedup_alignments,
            get_union_bams,
            callSignificantBases],
           add_inputs(os.path.join(PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
                      os.path.join(PARAMS["genome_dir"],
                                   PARAMS["genome"])),
           regex(".+/(.+)(.bam|.bed.gz"),
           [r"kmers.dir/\1.%smers.tsv.gz" % kmer
            for kmer in PARAMS["kmer_lengths"].split(",")])
def get_kmer_enrichments(infiles, outfiles):
    ''' Look for enrichment of kmers in crosslinked bases. Will
    examine both significant bases and all bases and look for 
    kmers of a length up that specified in ini'''
    
    bamfile, geneset, genome = infiles

    if bamfile.endswith(".bam"):
        sites = "-b " + bamfile
    elif bamfile.endswith(".bg.gz"):
        sites = "--bed=" + bamfile

    if PARAMS["use_centre"]==1:
        centre="--use_centre"
    else:
        centre=""
        
    statement_template = '''python %(project_src)s/scripts/iCLIP_kmer_enerichment.py
                                    %(sites)s
                                    %(centre)s
                                   -f %(genome)s
                                   -k %(kmer)s
                                   -s %%(cluster_window_size)s
                                   -p %%(kmers_threads)s
                                   -I %(geneset)s
                                   -S %(outfile)s
                                   -L %(outfile)s.log'''

    statements = []
    
    for kmer in PARAMS["kmer_lengths"].split(","):
        statements.append(statement_template % locals())

    P.run()


###################################################################
@collate(get_kmer_enrichments,
         formatter("[^\.]+\.[0-9]+mer.tsv.gz"),
         "kmers.dir/all_bases.kmers.load")
def load_dedup_kmers(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.+)-(.+)[-\.](R.+|union)\.([0-9]+mer).tsv.gz",
                         cat="factor,tag,replicate,k",
                         options="-i factor -i tag -i replicate -i k -kmer")


###################################################################   
@collate(get_kmer_enrichments,
         formatter("[^\.]+\.sig_bases.[0-9]+mer.tsv.gz"),
         "kmers.dir/all_bases.kmers.load")
def load_sig_base_kmers(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename="(.+)-(.+)[-\.](R.+|union)\.([0-9]+mer).tsv.gz",
                         cat="factor,tag,replicate,k",
                         options="-i factor -i tag -i replicate -i k -kmer")

    
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
    statement = '''python %(project_src)s/iCLIP2bigWig.py
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
           r"export/hg19/\1.bigBed")
def exportClusters(infile, outfile):
    ''' Add a track line to cluster files and export to export dir '''
   
    PipelineiCLIP.clustersToBigBed(infile, outfile)


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


@follows( mkdir( "report" ), createViewMapping)
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


@follows(mkdir("report"), createViewMapping)
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
