################################################################################
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

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

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
|reaper              |                   |Used for demuxing and clipping reads            |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_iCLIP.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_iCLIP.tgz
   tar -xvzf pipeline_iCLIP.tgz
   cd pipeline_iCLIP
   python <srcdir>/pipeline_iCLIP.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *
from ruffus.combinatorics import *

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGATPipelines.PipelineUtilities as PUtils
import CGATPipelines.PipelineMapping as PipelineMapping
import PipelineiCLIP
###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# define some tracks if needed
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob("*.ini" ), "(\S+).ini" )


###################################################################
###################################################################
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
###################################################################
## worker tasks
###################################################################
@transform("*.fastq.gz", regex("(.+).fastq.gz"),
           add_inputs(os.path.join(PARAMS["bowtie_index_dir"],
                                   PARAMS["phix_genome"]+".fa")),
           r"\1.fastq.clean.gz")
def filterPhiX(infiles, outfile):
    ''' Use mapping to bowtie to remove any phiX mapping reads '''

    job_options = "-pe dedicated %i -R y -l mem_free=%s" % (
        PARAMS["phix_bowtie_threads"], PARAMS["phix_bowtie_memory"])
    infile, reffile = infiles
    outfile = P.snip(outfile,".gz")
    m = PipelineMapping.Bowtie(executable=PARAMS["phix_bowtie_exe"],
                               strip_sequence=False,
                               remove_non_unique=False)
    genome = PARAMS["phix_genome"]
    bowtie_options = PARAMS["phix_bowtie_options"] + \
                     " --un %s" % outfile
    bowtie_threads = PARAMS["phix_bowtie_threads"]
    bam_out = P.snip(infile,".fastq.gz") + ".phix.bam"
    statement = m.build((infile,),bam_out)
    statement += "checkpoint; gzip %(outfile)s"
    P.run()


@transform("sample_table.tsv", suffix(".tsv"), ".load")
def loadSampleInfo(infile, outfile):

    P.load(infile, outfile,
           options="--header=format,barcode,track -i barcode -i track")
###################################################################
@follows(mkdir("demux_fq"))
@transform(filterPhiX, regex("(.+).fastq.clean.gz"),
           r"demux_fq/\1.fastq.umi_trimmed.gz")
def extractUMI(infile, outfile):
    ''' Remove UMI from the start of each read and add to the read
    name to allow later deconvolving of PCR duplicates '''

    statement=''' zcat %(infile)s
                | python %(project_src)s/extract_umi.py
                        --bc-pattern=%(reads_bc_pattern)s
                        -L %(outfile)s.log
                | gzip > %(outfile)s '''

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
@transform("sample_table.tsv", regex("(.+)"),
           r"reaper_metadata.tsv")
def generateReaperMetaData(infile, outfile):
    '''Take the sample_table and use it to generate a metadata table
    for reaper in the correct format '''

    adaptor_5prime = PARAMS["reads_5prime_adapt"]
    adaptor_3prime = PARAMS["reads_3prime_adapt"]

    outlines = []
    for line in IOTools.openFile(infile):
        barcode = line.split("\t")[1]
        outlines.append([barcode, adaptor_3prime, adaptor_5prime, "-"])

    header = ["barcode", "3p-ad", "tabu", "5p-si"]
    PUtils.write(outfile, outlines, header)


###################################################################
@follows(loadUMIStats)
@split(extractUMI, regex(".+/(.+).fastq.umi_trimmed.gz"),
       add_inputs(generateReaperMetaData, "sample_table.tsv"),
       r"demux_fq/*_\1.fastq.1.gz")
def demux_fastq(infiles, outfiles):
    '''Demultiplex each fastq file into a seperate file for each
    barcode/UMI combination'''

    infile, meta, samples = infiles
    track = re.match(".+/(.+).fastq.umi_trimmed.gz", infile).groups()[0]

    statement = '''reaper -geom 5p-bc
                          -meta %(meta)s
                          -i <( zcat %(infile)s | sed 's/ /_/g')
                          --noqc
                          -basename demux_fq/%(track)s_
                          -clean-length 15 > %(track)s_reapear.log;
                   checkpoint;
                   rename _. _ demux_fq/*clean.gz;
                 '''

    for line in IOTools.openFile(samples):
        line = line.split("\t")
        bc, name = line[1:]
        name = name.strip()
        statement += '''checkpoint;
                        mv demux_fq/%(track)s_%(bc)s.clean.gz
                           demux_fq/%(name)s_%(track)s.fastq.1.gz; ''' % locals()

    P.run()


###################################################################
@active_if(PARAMS["reads_paired"]==1)
@transform("*.fastq.2.gz", suffix(".fastq.2.gz"),
           ".fastq.reaped.2.gz")
def reapRead2(infile,outfile):

    track = P.snip(outfile,".fastq.reaped.2.gz")
    statement = ''' reaper -geom no-bc
                           -3pa %(reads_3prime_adapt)s
                           -i %(infile)s
                           -basename %(track)s
                           --noqc
                           -clean-length 15 > %(track)s_pair_reaper.log;
                    checkpoint;
                    mv %(track)s.lane.clean.gz %(outfile)s '''

    P.run()


###################################################################
@active_if(PARAMS["reads_paired"]==1)
@follows(mkdir("reconciled.dir"), reapRead2)
@transform(demux_fastq,
           regex(".+/(.+)_(.+).fastq.1.gz"),
           add_inputs(r"\2.fastq.2.gz"),
           [r"reconciled.dir/\1_\2.fastq.2.gz",
            r"reconciled.dir/\1_\2.fastq.1.gz"])
def reconsilePairs(infiles, outfiles):
    ''' Pull reads read 2 file that are present in each read 1 file '''

    track = P.snip(os.path.basename(infiles[0]), ".fastq.1.gz")
    infiles = " ".join(infiles)
    job_options = "-l mem_free=1G"
    statement = '''python %(scripts_dir)s/fastqs2fastqs.py
                          --method=reconcile
                          --id-pattern-1='(.+)_.+_[ATGC]+'
                          --output-pattern=reconciled.dir/%(track)s.fastq.%%s.gz
                           %(infiles)s > reconciled.dir/%(track)s.log '''

    P.run()
      
###################################################################

###################################################################
@follows(mkdir("fastqc"))
@transform(demux_fastq, regex(".+/(.+).fastq.gz"),
           r"fastqc/\1.fastqc")
def qcDemuxedReads(infile, outfile):
    ''' Run fastqc on the post demuxing and trimmed reads'''

    m = PipelineMapping.FastQc(nogroup=False)
    statement = m.build((infile,),outfile)
    exportdir = "."
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

    infiles = glob.glob("%s/*.fastq.1.gz" % PARAMS["input"])

    outfiles = ["mapping.dir/%(mapper)s.dir/%(track)s.%(mapper)s.bam"
                % {"mapper": PARAMS["mappers"],
                   "track": re.match("%s/(.+).fastq.1.gz" % PARAMS["input"], infile).groups()[0]}
                for infile in infiles]

    yield (infiles, outfiles)


###################################################################
@follows(mkdir("mapping.dir"), demux_fastq)
@files(mapping_files)
def run_mapping(infiles, outfiles):
    ''' run the mapping target of the mapping pipeline '''

    to_cluster = False
    statement = ''' ln -f pipeline.ini mapping.dir/pipeline.ini;
                    checkpoint;
                    cd mapping.dir;
                    nice python %(scripts_dir)s/../CGATPipelines/pipeline_mapping.py
                    make mapping
                    -v5 -p%(pipeline_mapping_jobs)s  '''

    P.run()


###################################################################
@merge(run_mapping, "mapping.sentinal")
def mapping_qc(infiles, outfile):
    ''' run mapping pipeline qc targets '''

    to_cluster = False

    statement = '''cd mapping.dir;
                   nice python %(scripts_dir)s/../CGATPipelines/pipeline_mapping.py
                   make qc -v5 -p%(pipeline_mapping_jobs)s '''
    P.run()

    P.touch("mapping.sentinal")


###################################################################
@follows(mapping_qc)
@transform("mapping.dir/geneset.dir/reference.gtf.gz",
           regex(".+/(.+).gtf.gz"),
           r"\1.context.bed.gz")
def generateContextBed(infile, outfile):
    ''' Generate full length primary transcript annotations to count
    mapping contexts '''

    genome = os.path.join(PARAMS["annotations_dir"], "contigs.tsv")
    statement = ''' zcat %(infile)s
                  | python %(scripts_dir)s/gtf2gtf.py
                    --exons2introns
                    --with-utr
                     -L %(outfile)s.log
                  | awk 'BEGIN{FS="\\t";OFS="\\t"} {$2="intron"; print}'
                  | gzip > %(outfile)s.tmp.gtf.gz;

                  checkpoint;

                  zcat %(infile)s %(outfile)s.tmp.gtf.gz           
                  | python %(scripts_dir)s/gff2bed.py
                    --name=source
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
                            --per-name
                            -I %(infile)s
                    | gzip > %(outfile)s '''

    P.run()


###################################################################
@transform(getContextIntervalStats, suffix(".tsv.gz"),
           ".load")
def loadContextIntervalStats(infile, outfile):

    P.load(infile, outfile)


###################################################################
@transform(mapping_qc, regex("(.+)"),
           "mapping.dir/view_mapping.load")
def createViewMapping(infile, outfile):
    ''' Create tables neccessary for mapping report '''

    to_cluster = False
    statement = '''cd mapping.dir;
                   nice python %(scripts_dir)s/../CGATPipelines/pipeline_mapping.py
                   make createViewMapping -v5 -p1 '''
    P.run()


###################################################################
@collate("mapping.dir/*.dir/*.bam",
         regex("(.+)/(.+\-.+\-[^_]+)_(.+)\.([^.]+)\.bam"),
         r"\1/\2.\4.bam")
def mergeBAMFiles(infiles, outfile):
    '''Merge reads from the same library, run on different lanes '''

    if len(infiles) == 1:
        P.clone(infiles[0], outfile)
        P.clone(infiles[0] + ".bai", outfile + ".bai")
        return

    infiles = " ".join(infiles)
    statement = ''' samtools merge %(outfile)s %(infiles)s >& %(outfile)s.log;

                    checkpoint;

                    samtools index %(outfile)s '''

    P.run()


###################################################################
@follows(mapping_qc,loadContextIntervalStats )
def mapping():
    pass


###################################################################
# Deduping, Counting, etc
###################################################################
@follows(mkdir("deduped.dir"), run_mapping)
@transform(mergeBAMFiles, regex(".+/(.+)\.[^\.]+\.bam"),
           r"deduped.dir/\1.bam")
def dedup_alignments(infile, outfile):
    ''' Deduplicate reads, taking UMIs into account'''

    outfile = P.snip(outfile, ".bam")

    statement = ''' python %(project_src)s/dedup_umi.py
                    -I %(infile)s
                    -S %(outfile)s.tmp.bam
                    -L %(outfile)s.log;
 
                    checkpoint;

                    samtools sort %(outfile)s.tmp.bam %(outfile)s;
                   
                    checkpoint;

                    samtools index %(outfile)s.bam ;
 
                    checkpoint;

                    rm %(outfile)s.tmp.bam'''

    P.run()


###################################################################
@transform([dedup_alignments,mergeBAMFiles], suffix(".bam"),
           ".frag_length.tsv")
def getFragLengths(infile, outfile):
    ''' estimate fragment length distribution from read lengths'''

    statement = ''' python %(project_src)s/length_stats.py
                           -I %(infile)s
                           -S %(outfile)s
                           -L %(outfile)s.log
               '''
    if PARAMS["reads_paired"] == 1:
        statement += "--paired=%s" % (
            int(PARAMS["reads_length"]) - len(PARAMS["reads_bc_pattern"]))
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
                         --force
                          < %(infile)s > %(outfile)s '''

    P.run()


###################################################################
@merge(dedupedBamStats, "deduped_bam_stats.load")
def loadDedupedBamStats(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".+/(.+).bam_stats.tsv")


###################################################################
@transform(dedup_alignments,
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
@follows(mkdir("saturation.dir"), run_mapping)
@split(mergeBAMFiles, regex(".+/(.+)\.[^\.]+\.bam"),
       [r"saturation.dir/\1.%.3f.bam" % (1.0/(2 ** x))
        for x in range(1, 6)] +
       [r"saturation.dir/\1.%.3f.bam" % (x/10.0)
        for x in range(6, 11, 1)])
def subsetForSaturationAnalysis(infile, outfiles):
    '''Perform subsetting of the original BAM files, dedup and
    return the context stats. Test for resturn on investment for
    further sequencing of the same libraries '''

    track = re.match(".+/(.+)\.[^\.]+\.bam", infile).groups()[0]

    statements = []
    statement_template = '''
                            python %%(project_src)s/dedup_umi.py
                              -I %(infile)s
                              -L %(outfile)s.log
                              -S %(outfile)s.tmp.bam
                              --subset=%(subset).3f;
                             checkpoint;
                             samtools sort %(outfile)s.tmp.bam %(outfile)s;
                             checkpoint;
                             samtools index %(outfile)s.bam '''
    for x in range(1, 6):
        subset = 1.0/(2 ** x)
        outfile = "saturation.dir/%%(track)s.%.3f" % subset
        statements.append(statement_template % locals())

    for x in range(6, 11):
        subset = x/10.0
        outfile = "saturation.dir/%%(track)s.%.3f" % subset
        statements.append(statement_template % locals())

    P.run()


###################################################################
@transform(subsetForSaturationAnalysis, suffix(".bam"), ".bamstats.tsv")
def subsetBamStats(infile, outfile):
    ''' Stats on the subset BAMs '''

    job_options = "-l mem_free=500M"
    statement = ''' python %(scripts_dir)s/bam2stats.py 
                    --force < %(infile)s > %(outfile)s '''
    P.run()


###################################################################
@merge(subsetBamStats, "subset_bam_stats.load")
def loadSubsetBamStats(infiles, outfile):
    P.concatenateAndLoad(infiles, outfile,
                         regex_filename= ".+/(.+-.+-.+)\.([0-9]+\.[0-9]+).bamstats.tsv",
                         cat="track,subset")


###################################################################
@transform([mergeBAMFiles, dedup_alignments, subsetForSaturationAnalysis],
           suffix(".bam"),
           add_inputs(generateContextBed),
           ".reference_context.tsv")
def buildContextStats(infiles, outfile):
    ''' Find context from reads '''

    infile, reffile = infiles

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
         loadSubsetBamStats,
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

    statement = '''python %(project_src)s/calculateiCLIPReproducibility.py
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
# Analysis
###################################################################
@follows(mkdir("gene_profiles.dir"))
@transform(dedup_alignments, regex(".+/(.+).bam"),
           add_inputs("mapping.dir/geneset.dir/refcoding.gtf.gz"),
           r"gene_profiles.dir/\1.tsv")
def calculateGeneProfiles(infiles, outfile):
    ''' Calculate metagene profiles over protein coding genes
    for each sample'''

    infile, reffile = infiles
    statement = '''python %(scripts_dir)s/bam2geneprofile.py
                           --method=geneprofilewithintrons
                           --bamfile=%(infile)s
                           --gtffile=%(reffile)s
                           --normalization=total-sum
                           --normalize-profile=area
                           --log=%(outfile)s.log
                           --output-filename-pattern=%(outfile)s.%%s
                           > %(outfile)s '''

    P.run()


###################################################################
@transform("mapping.dir/geneset.dir/refcoding.gtf.gz",
           suffix(".gtf.gz"),
           ".exons.gtf.gz")
def transcripts2Exons(infile, outfile):
    ''' Make each exon a seperate gene to allow quantitaion over exons
    only keeps those exons longer than 100 base pairs'''
    
    tmp_outfile = P.snip(outfile, ".gtf.gz") + ".tmp.gtf.gz"
    PipelineiCLIP.removeFirstAndLastExon(infile, tmp_outfile)

    statement = '''python %(scripts_dir)s/gff2bed.py 
                          --is-gtf 
                           -I %(tmp_outfile)s 
                           -L %(outfile)s.log
                  | sort -k1,1 -k2,2n
                  | mergeBed -s -d 100 -nms
                  | awk '($3-$2) > 100 {print}'
                  | awk 'BEGIN{OFS="\\t"} {$4=NR; print}'
                  | python %(scripts_dir)s/bed2gff.py --as-gtf -L %(outfile)s.log
                  | gzip -c > %(outfile)s '''
    P.run()

    os.unlink(tmp_outfile)


###################################################################
@transform("mapping.dir/geneset.dir/refcoding.gtf.gz",
           suffix(".gtf.gz"),
           ".introns.gtf.gz")
def transcripts2Introns(infile, outfile):
    ''' Make each exon a seperate gene to allow quantitaion over exons'''

    tmp_outfile = P.snip(outfile, ".gtf.gz") + ".tmp.gtf.gz"
    PipelineiCLIP.removeFirstAndLastExon(infile, tmp_outfile)

    statement = '''python %(scripts_dir)s/gtf2gtf.py
                           -I %(tmp_outfile)s
                          --log=%(outfile)s.log
                           --exons2introns
                  | python %(scripts_dir)s/gff2bed.py
                          --is-gtf
                           -L %(outfile)s.log
                  | sort -k1,1 -k2,2n
                  | mergeBed -s -d 100 -nms
                  | awk 'BEGIN{OFS="\\t"} {$4=NR; print}'
                  | python %(scripts_dir)s/bed2gff.py
                          --as-gtf -L %(outfile)s.log
                  | gzip -c > %(outfile)s '''
    P.run()

    os.unlink(tmp_outfile)

###################################################################
@product(dedup_alignments,
         formatter(".+/(?P<TRACK>.+).bam"),
         [transcripts2Exons, transcripts2Introns],
         formatter(".+/refcoding.(?P<INTERVALTYPE>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK[0][0]}.{INTERVALTYPE[1][0]}.log")
def calculateExonProfiles(infiles, outfile):

    infile, reffile = infiles

    outfile = P.snip(outfile, ".log")
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=intervalprofile
                          --bamfile=%(infile)s
                          --gtffile=%(reffile)s
                          --normalization=total-sum
                          --base-accuracy
                          --normalize-profile=area
                          --resolution-upstream=50
                          --resolution-downstream=50
                          --extension-upstream=50
                          --extension-downstream=50
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.log '''

    P.run()


###################################################################
@product(dedup_alignments,
         formatter(".+/(?P<TRACK>.+).bam"),
         [transcripts2Exons, transcripts2Introns],
         formatter(".+/refcoding.(?P<INTERVALTYPE>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK[0][0]}.{INTERVALTYPE[1][0]}.tssprofile.log")
def calculateExonTSSProfiles(infiles, outfile):

    infile, reffile = infiles

    outfile = P.snip(outfile, ".tssprofile.log")
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=tssprofile
                          --bamfile=%(infile)s
                          --gtffile=%(reffile)s
                          --normalization=total-sum
                          --base-accuracy
                          --normalize-profile=area
                          --resolution-upstream=100
                          --resolution-downstream=100
                          --extension-inward=100
                          --extension-outward=100
                          --output-filename-pattern=%(outfile)s.%%s
                          --log=%(outfile)s.tssprofile.log '''

    P.run()


###################################################################
@follows(calculateExonTSSProfiles,
         calculateGeneProfiles,
         calculateExonProfiles)
def profiles():
    pass


###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(PrepareReads, mapping, MappingStats, reproducibility,
         profiles)
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
    statement = '''cd mapping.dir;
                   python %(scripts_dir)s/../CGATPipelines/pipeline_mapping.py
                   -v5 -p1 make build_report '''
    P.run()
    E.info("starting report build process from scratch")
    P.run_report( clean = True )


@follows(mkdir("report"), createViewMapping)
def update_report():
    '''update report.'''

    E.info("Updating Mapping Reported")
    statement = '''cd mapping.dir;
                   python %(scripts_dir)s/../CGATPipelines/pipeline_mapping.py
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
