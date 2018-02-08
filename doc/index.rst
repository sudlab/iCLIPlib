.. iCLIPlib documentation master file, created by
   sphinx-quickstart on Thu Jan  4 15:21:51 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. currentmodule:: iCLIP
Welcome to iCLIPlib's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


This library presents a series of tools for manipulating data from RNA
Cross-Linking and ImmunoPrecipitation (CLIP) like experiments. Also includes
example scripts performing common analyses.  Most of the
code was developed with iCLIP in mind. This documentation is currently 
mostly a stub while I get my act together to produce more complete narrative
documentation. 

Installation
------------

Currently installation relies mostly on the user to deal with the dependencies.
These are currently:

   * CGAT (www.github.com/CGATOxford/cgat)
     
   * numpy
     
   * pysam
     
   * pandas

   * bx-python

After installation of the dependencies, run `python setup.py install`, and put
the contents of the scripts directory somewhere your shell can find them.

Look out very soon for proper dependency management and installers from
conda and PyPI.

  
Getter functions
----------------

The key functions are the *_getter functions. These functions are used to 
return Series of signal across a specified genomic region. The usual way to 
use these functions is to pass the file(s) containing the data to the
factory function which will return a function that implements expected 
interface and can then be passed to various other functions without those
functions being aware of the source of the data.  

.. autosummary::
    :toctree: generated/
    
    getters.make_getter
    getters.getter
    
the following should not need to be called directly but contain details of
how bamfiles are converted into cross-link profiles:

.. autosummary::
    :toctree: generated/
    
    counting.getCrosslink
    counting.find_first_deletion


Counter functions
-----------------
 
The getter functions themselves return conuts over regions. These functions
return `pandas.Series` over various objects, such as lists of genomic
intervals or transcripts. Transcripts are genearally handled by the
CGAT.GTF interface.
 
.. autosummary::
   :toctree: generated/
  
   counting.count_intervals
   counting.count_transcript
    
Metagene profiles
-----------------

Metagenes are representations of average profiles. Counts can be binned
into a set number of bins across a transcript or using a set bin width. 
Various normalisation approaches can be applied.

.. autosummary::
    :toctree: generated/
    
    meta.meta_gene
    meta.get_binding_matrix
    meta.processing_index
    meta.get_window
    
users with more complex requirements might also want to look at the
underlying binning function.

.. autosummary::
    :toctree: generated/
    
    meta.bin_counts
    
normalisation functions:

.. autosummary::
    :toctree: generated/
    
    meta.quantile_row_norm
    meta.sum_row_norm
    meta.compress_matrix
    
Randomisation
-------------

Since we know very little about the statistical distributions involved with
many aspects of this sort of data, the signficance of patterns is often
measured by comparison to randomised patterns. These tools help create that
type of analysis:

.. autosummary::
    :toctree: generated/
    
    random.bootstrap
    random.boot_ci
    random.count_ratio
    random.ratio_and_ci
    utils.randomiseSites
    utils.rand_apply

Coordinates Converter
---------------------

There are many good reasons to align reads to the genome. But many of our
analyses take place on spliced transcripts. The 
:class:`TranscriptCoordInterconverter` is a tool for moving between these
two worlds, converting genome coords to transcript coords and vice-versa

.. autosummary::
    :toctree: generated/
    
    iCLIP.utils.TranscriptCoordInterconverter
    

Clusters
--------

Although there are now many more advanced ways to calculate clusters, the
method described by Wang et al, PLoS Biol. 2010:e1000530 :PMID:`2104891`
remains the most popular. These tools help to produce analyses of this
type:

.. autosummary::
    :toctree: generated/
    
    clusters.Ph
    clusters.fdr
    clusters.get_crosslink_fdr_by_randomisation

    
Measures of similarity and difference between profiles
-------------------------------------------------------

There is no agreed best way to measure distance between profiles, but these
techniques should give you some ideas. Warning, they can be slow...

.. autosummary::
    :toctree: generated/
    
    distance.calcAverageDistance
    distance.findMinDistance
    distance.corr_profile
    
 These methods are implemented as a full functioning script that compares
profiles across whole transcript sets/genomes in 
`reproducibility_by_exon.py`.

Analysis of kmers
-----------------

Tools for examining the distribution and enrichment of kmers and other
landmarks around clipped bases

.. autosummary::
    :toctree: generated/
    
    kmers.find_all_matches
    kmers.pentamer_frequency
    kmers.pentamer_enrichment
    
Other
------

.. autosummary::
    :toctree: generated/
    
    utils.spread
    
 
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
