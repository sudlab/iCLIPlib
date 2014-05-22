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

:Author: Andreas Heger
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

import sys, glob, gzip, os, itertools, re, math, types, collections, time
import optparse, shutil
import sqlite3

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGATPipelines.PipelineUtilities as PUtils
import CGATPipelines.PipelineMapping as PipelineMapping

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
    statement += "; checkpoint; gzip %(outfile)s"
    P.run()


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
                           demux_fq/%(name)s_%(track)s.fastq.gz; ''' % locals()

    P.run()


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
@follows(demux_fastq,qcDemuxedReads, loadUMIStats)
def PrepareReads():
    pass


###################################################################
# Mapping
###################################################################
def mapping_files():

    infiles = glob.glob("demux_fq/*.fastq.gz")

    outfiles = ["mapping.dir/%(mapper)s.dir/%(track)s.%(mapper)s.bam"
                % {"mapper": PARAMS["mappers"],
                   "track": re.match("demux_fq/(.+).fastq.gz", infile).groups()[0]}
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
                    --merge-transcripts
                    --with-utr
                     -L %(outfile)s.log
                  | python %(scripts_dir)s/gff2bed.py
                    --name=source
                     -L %(outfile)s.log
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
                  
                    rm %(outfile)s.tmp.bed %(outfile)s.tmp2.bed '''

    P.run()


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
@follows(mapping_qc)
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

    statement = ''' python %(project_src)s/dedup_umi.py
                    -I %(infile)s
                    -S %(outfile)s
                    -L %(outfile)s.log;
 
                    checkpoint;

                    samtools index %(outfile)s '''

    P.run()


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
                              -S %(outfile)s
                              --subset=%(subset).3f;
                             checkpoint;
                             samtools index %(outfile)s '''
    for x in range(1, 6):
        subset = 1.0/(2 ** x)
        outfile = "saturation.dir/%%(track)s.%.3f.bam" % subset
        statements.append(statement_template % locals())

    for x in range(6, 11):
        subset = x/10.0
        outfile = "saturation.dir/%%(track)s.%.3f.bam" % subset
        statements.append(statement_template % locals())

    P.run()


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
                           > %(outfile)s '''

    P.run()


###################################################################
@transform("mapping.dir/geneset.dir/refcoding.gtf.gz",
           suffix(".gtf.gz"),
           ".exons.gtf.gz")
def transcripts2Exons(infile, outfile):
    ''' Make each exon a seperate gene to allow quantitaion over exons'''


    statement = '''python %(scripts_dir)s/gff2bed.py 
                          --is-gtf 
                           -I %(infile)s 
                           -L %(outfile)s.log
                  | awk 'BEGIN{OFS='\\t'} {$4=NR; print}'
                  | python %(scripts_dir)s/bed2gff.py --is-gtf -L %(outfile)s.log
                  | gzip -c > %(outfile)s '''
    P.run()


###################################################################
@transform("mapping.dir/geneset.dir/refcoding.gtf.gz",
           suffix(".gtf.gz"),
           ".introns.gtf.gz")
def transcripts2Introns(infile, outfile):
    ''' Make each exon a seperate gene to allow quantitaion over exons'''


    statement = '''python %(scripts_dir)s/gff2bed.py 
                          --is-gtf 
                           -I %(infile)s 
                           -L %(outfile)s.log
                  | awk 'BEGIN{OFS='\\t'} {$4=NR; print}'
                  | python %(scripts_dir)s/bed2gff.py --is-gtf -L %(outfile)s.log
                  | gzip -c > %(outfile)s '''
    P.run()


###################################################################
@product(dedup_alignments, 
         formatter(".+/(?P<TRACK>).bam"),
         [transcripts2Exons, transcripts2Introns], 
         formatter(".+/refcoding.(?<INTERVALTYPE>.+).gtf.gz"),
         "gene_profiles.dir/{TRACK}[0][0].INTERVALTYPE[1][0].log")
def calculateExonProfiles(infiles, outfile):

    infile, reffile = infiles

    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                          --method=intervalprofile
                          --bamfile=%(infile)s
                          --gtffile=%(reffile)s
                          --normalization=total-sum
                          --normalize-profile=area
                          --log=%(outfile)s.log '''

    P.run()

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows(PrepareReads, mapping)
def full(): 
    pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( mkdir( "report" ), createViewMapping )
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
