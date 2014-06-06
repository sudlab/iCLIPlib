Distribution of UMIs
---------------------


.. report:: Sample_QC.UMI_stats
   :render: r-ggplot
   :statement: aes(x=track,y=freq) + geom_violin() + scale_y_log10() + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + geom_hline(yintercept=1/(4^5), lty=2)

   Distribution of UMIs


Reads Per Sample
-----------------


.. report:: Sample_QC.ReadsPerSample
   :render: r-ggplot
   :statement: aes(y=total, x=track) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x,...) format(x,...,big.mark=",", scientific= F, trim = T)) + ylab("Reads")

   Number of reads for each barcode


.. report:: Sample_QC.PercentDemuxed
   :render: r-ggplot
   :statement: aes(y=demuxed * 100, x=track) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x) sprintf("%.0f%%",x)) + ylab("Percent Passed Filter")

   Reads succesfully demuxed and filtered


.. report:: Sample_QC.PercentMapped
   :render: r-ggplot
   :statement: aes(y=mapped * 100, x=track) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x) sprintf("%.0f%%",x), limits = c(0,100)) + ylab("Percent reads mapped")

   Percent Reads mapped



.. report:: Sample_QC.PercentDeDuped
   :render: r-ggplot
   :statement: aes(y=p_unique * 100, x=track) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x) sprintf("%.0f%%",x)) + ylab("Percent reads unique")

   Percent of Reads Unique



.. report:: Sample_QC.FinalReads
   :render: r-ggplot
   :statement: aes(y=reads_mapped, x=track) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle=90)) + scale_y_continuous(labels = function(x,...) format(x,...,big.mark=",", scientific= F, trim = T)) + ylab("Total unique mapped reads")

   Final total mapped reads

.. report:: Sample_QC.PercentSpliced
   :render: r-ggplot
   :statement: aes(Track, pspliced) + geom_bar(stat="identity") +  scale_y_continuous(labels = function(x) sprintf("%.0f%%",x*100)) + ylab("Percent reads spliced") + theme_bw() + theme(axis.text.x=element_text(angle=90))

   Percent of deduped reads spliced

.. report:: Sample_QC.ReadLengths
   :render: r-ggplot
   :transform: melt
   :statement: aes(x=Track, y=Data/sum(Data), fill=Slice) + geom_bar(position="fill", stat="identity") + ylab("Fraction of reads") + scale_fill_discrete(name="Length bin (bp)") + coord_flip() + theme_bw()

   Read length distribution of unmapped reads




Saturation Analysis
--------------------

.. report:: Sample_QC.AlignmentSaturation
   :render: r-ggplot
   :transform: label-paths
   :statement: aes(x=subset, y=counts, color = factor, shape = factor) + geom_point() + geom_line() + facet_wrap(~replicate) + theme_bw() + theme(aspect.ratio = 1)

   Subsampling of alignments


.. report:: Sample_QC.AlignmentSaturation
   :render: r-ggplot
   :transform: label-paths
   :statement: aes(x=counts/subset, y=counts, color = factor, shape = factor) + geom_point() + geom_line() + facet_wrap(~replicate) + theme_bw() + theme(aspect.ratio = 1)

   Tests for model assumptions


.. report:: Sample_QC.LibrarySize_Binom
   :render: r-ggplot
   :statement: aes(x=subset, y=alignments) + geom_point() + geom_line(aes(y=expected_unique)) + geom_hline(yintercept=rframe$library_size[1]) + theme_bw()
   :width: 200
   :layout: column-4
   

   curve fits for saturation using Binomal distribution



.. report:: Sample_QC.LibrarySize_mm
   :render: r-ggplot
   :statement: aes(x=subset, y=alignments) + geom_point() + geom_line(aes(y=expected_unique)) + geom_hline(yintercept=rframe$library_size[1]) + theme_bw()
   :width: 200
   :layout: column-4

   curve fits for saturation using reciprical fit


.. report:: Sample_QC.mm_fit_stats
   :render: r-ggplot
   :transform: label-paths
   :statement: aes(x=Slice, y=Library.Size) + geom_bar(stat="identity") + geom_bar(aes(y=Library.Size*Percent.Saturation/100), stat="identity", fill = "orange") + theme_bw() + theme(axis.text.x = element_text(angle=90))

   Library size estimates

Context Stats
---------------

.. report:: Sample_QC.ContextStats
   :render: pie-plot
   :layout: column-4

   Mapping Contexts for deduped reads


.. report:: Sample_QC.ContextRepresentation
   :render: r-ggplot
   :statement: aes(category, log2(precent_alignments/percent_bases)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=90,hjust=1)) + ylab("log2 enrichment")
   :layout: column-3
   :groupby: track

   Enrichments of contexted over expectation


Reproducibility
----------------

Reproducilbity measures the number of sites with at least n reads mapping to them in one replicate that have reads mapping to them in 1 or 2 of the other replicates as a fraction of the total number of sites with that depth in that replicate. 

.. report:: Sample_QC.Reproducibility
   :render: r-ggplot
   :transform: label-paths
   :statement: aes(level, reproducibility, color=Replicate) + geom_line() + geom_point() + facet_grid(Slice ~ Track) + coord_cartesian(xlim=c(0,5)) + theme_bw()
   :tf-label-level: 3

   Reproduciblity


The problem with the measure above (which is the one outlined in Sutomui et al) is that the largest rep will always have a lower reproducibility because all those extra locations can't possibly be replicated. Below I normalise the reproduciblility by the maximum possible level of reproduciblitity.

.. report:: Sample_QC.NormReproducibility
   :render: r-ggplot
   :transform: label-paths
   :statement: aes(level, reproducibility, color=Replicate) + geom_line() + geom_point() + facet_grid(Slice ~ Track) + coord_cartesian(xlim=c(0,5), ylim=c(0,1)) + theme_bw()
   :tf-label-level: 3

   Normalised Reproduciblity

The next plot shows how reproducible cross-linked bases are in the control samples rather than in other replicates of the sample cell line. 

.. report:: Sample_QC.ReproducibilityVsControl
   :render: r-ggplot
   :transform: label-paths
   :slices: 1,3
   :statement: aes(level,reproducibility,color=Replicate) + geom_line() + geom_point() + facet_grid(Slice~Track) + theme_bw() + coord_cartesian(xlim = c(0,25))
   :tf-label-level: 3

   Reproducibility vs. Controls


Given that there is some reproducibility between one replicate and others pulling down the same factors and also some between that same replicate and the negative controls, how much infomation is there in the sample that is due to the correct pull down. Assuming that infomation shared between a sample and a control will also be shared by another replicate of the same pull down, the ratio of replicating bases between A) a replicate and a the control and B) one replicate and another should be above one, and the excess should speak to how much extra, factor specific information there is. 


.. report:: Sample_QC.ReproducibilityReplicateVsControl
   :render: r-ggplot
   :transform: label-paths
   :statement: aes(depth,ratio,color=Replicate) + geom_line() + geom_point() + facet_grid(Slice~Track) + scale_y_log10() + coord_cartesian(xlim=c(0,10)) + theme_bw()
   :tf-label-level: 3
   :slices: 1

   Ratio of reproducibility in replicates of same factor to that in other factors.



The reproducibility can also be used to calculate a distance metric between samples. The jaccard index is the interection of two sets divided by the union. By applying this accross each pair of samples at the 1 level we can build a clustering of samples.

