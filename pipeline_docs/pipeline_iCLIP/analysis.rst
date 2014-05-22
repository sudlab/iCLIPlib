

Interpretive Analysis
=====================

Introduction
-------------


This report comes in two parts. The first part is the pipeline section. This serves as a sort of data appendix/lab book which contains a large amount of data and analysis results with very little explaination/interpretation. The second part of the report is the interpretive section which will contain disucission and interpretation of the results as well as summaries.

I will endevour to update this report as new data and analysis become availible. 

Currently this report is protected by obfuscation - that it only those we choose know the address, and it is not indexed by google. However, this does not amount to proper security as it would be quite possible to discover the report if a concerted effort was made. We do offer a password protection facility if it become neccessary.

Contents
---------

.. toctree::
   :maxdepth: 2

   analysis/mapping.rst
   analysis/saturation.rst
   analysis/context.rst
   analysis/reproducibility.rst


Executive Summaries
-------------------

14/5/2014 - Pilot Experiment Assessment
+++++++++++++++++++++++++++++++++++++++

Read processing and mapping
............................

*  The vast majority of reads were sucessfull procesed.
*  There is a a strong inbalance in the number of reads from each of the 12 libraries, with replicate 1, and Nxf1 libaries having substaniallly fewer reads.
*  The usage of UMIs is not random, this may limited their usefulness. 
*  Dispite this, estimates of library exhaustion are greater for libraries from replicate 1, Nxf1 and the negative cotrol samples. 
    * Library exhaustion may be over estimated due to the non-random nature of UMI usage. 
* For more info see :ref:`mapping` and :ref:`saturation`

Mapping contexts
.................

*  For Alyref and Chtop there is a strong enrichment of reads mapping to exons as opposed to introns, and a proportion of reads map over splice junctions, suggesting these factors bind spliced transcripts rather than primary transcripts. 
*  There is also a strong enrichment for binding to none lincRNA-like non-coding RNA (such as rRNA, snoRNA etc) in all samples. 
*  For more info see :ref:`context`

Reproducibility
.................

*  The levels of reproducibility between replicates for all samples increases as the number of reads mapping to each crosslink site increases.
*  However, links are also more likely to reproduce in the control samples as depth increases. 
*  The reproducibility in replicates is higher than the reproduciblity in controls for Alyref and Chtop, but not so much for Nxf1
*  This is confirmed by clustering on reproduciblity, where replicates 2 and 3 of Alyref and Chtop cluster together, as do replicate 1 of these factors plus replicates 2 and 3 of Nxf1 and one of the negative controls. 
*  For full details see :ref:`reproducibility`

Gene Profiles
...............

*  The distribution of reads over transcripts suggests that Alyref may bind at the beginning of transcripts and Chtop at the end. 
*  Full details see :ref:`metagene`

Recommendation
................

*  I recommend we ditch replicate 1 and go ahead with further sequencing for replicates 2 and 3, plus the un-tested replicate 4.


