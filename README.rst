Tools for dealing with iCLIP data
==================================

This package contains a selection of programming tools for dealing with iCLIP data
and a bunch of scripts based on them. 

There are several pipelines available for dealing with iCLIP data - turning reads into
significant bases, but often we any iCLIP project will require a number of non-standard 
analysis. Here we have a set of tools for doing ad-hoc analysis.

In general profiles are represented as pandas Series where the index represents 
genomic bases and the values represents tag counts. 

Documentation at for the classes and modules at readthedocs `here <http://icliplib.readthedocs.io/en/latest/index.html>`_. At the moment doucmenation for the scripts is in the `--help`. Proper documentation at readthedocs coming soon.

Proper installation and distribution coming very soon. 

The principle functions that produce these are:
    * count_intervals
    * count_transcript

Useful tools for analysis are:
    * randomizeSites - randomise the sites in a profile
    * rand_apply - randomize a profile a number of times and apply a function to it
    * spread - take a profile and extend each tag some bases in each direction
    
Also useful is 
   * TranscriptCoordInterconverter - as class for converting genomic coordinates to transcript ones

In addition to this are implementations for a number of published analyses:
   * kmer_enrichment - for looking for enriched kmers compared to randomised profiles
   * meta_gene - calculate a meta gene profile 
   * get_crosslink_fdr_by_randomisation - find significantly crosslinked bases by comparison to randomised profiles
     
and a number of scripts that implement example analyses. 

The requirements for the module are cgat-apps, pysam, numpy, bxpython and pandas, all installable from conda/bioconda
The pipeline has numerous other requirements. 

The documents are currently in the docs folder and in the code. Look out for them on readthedocs very soon.



  

