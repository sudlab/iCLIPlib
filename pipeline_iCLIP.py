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
@files("sample_table.tsv", "barcodes_table.tsv")
def generateBarcodesTable(infile, outfile):
    ''' Take specification of barcodes for each sample, which may
    contain Ns in UMI positions and generate all possible barcodes
    by subistuting Ns for all possibilities, and outputs repear
    metadata file '''

    adaptor_5prime = PARAMS["reads_5prime_adapt"]
    adaptor_3prime = PARAMS["reads_3prime_adapt"]

    DNA_bases = 'ATCG'
    outlines = []
    for line in IOTools.openFile(infile):
        barcode = line.split("\t")[0]
        base_list = []
        for base in barcode:
            if base == "N":
                base_list.append(DNA_bases)
            else:
                base_list.append(base)
        possible_barcodes = ["".join(x) for x in itertools.product(*base_list)]
        for bc in possible_barcodes:
            outlines.append([bc, adaptor_5prime, adaptor_3prime, "-"])

    header = ["barcode", "3p-ad", "tabu", "5p-si"]
    PUtils.write(outfile, outlines, header)
###################################################################


@follows(mkdir("demuxed_fq"))
@split("*.fastq.gz", regex("(.+).fastq.gz"),
       add_inputs(generateBarcodesTable),
       r"demuxed_fq/\1_*.fastq.gz")
def demux_fastq(infiles, outfiles):
    '''Demultiplex each fastq file into a seperate file for each
    barcode/UMI combination'''

    infile, meta = infiles
    track = re.match("(.+).fastq.gz", infile).groups()[0]

    statement = '''reaper -geom 5p-bc
                          -meta %(meta)s
                          -i %(infile)s
                          --noqc
                          -basename demuxed_fq/%(track)s_
                          -clean-length 5;
                   checkpoint;
                   rename _. _ demuxed_fq/*clean.gz;
                   checkpoint;
                   rename .clean .fastq demuxed_fq/*'''

    P.run()


###################################################################

@follows(mkdir("mapping.dir"))
@merge(demux_fastq, "mapping.sentinal")
def run_mapping(infiles, outfiles):
    ''' run the mapping target of the mapping pipeline '''

    to_cluster = False
    statement = ''' ln -f pipeline.ini mapping.dir/pipeline.ini;
                    checkpoint;
                    cd mapping.dir; 
                    nice python %(scripts_dir)s/../CGATPipelines/pipeline_mapping.py
                    make mapping
                    -v5 -p%(pipeline_mapping_jobs)s > ../mapping.sentinal '''

    P.run()

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows( )
def full(): pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish():
    '''publish report and data.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )
