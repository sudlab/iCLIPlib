library(rtracklayer)
library(DEXSeq)
library(BiocParallel)
library(experimentr)
library(GenomicRanges)
args = Experiment.start()
gffs <- import(args$infiles[1])
counts <- read.delim(args$infiles[2])
counts <- counts[,-1]
counts <- as.matrix(counts)
col_data = read.delim(args$infile[3], header=T, row.names = 1)
col_data$condition <- relevel(col_data$condition, "Control")
cat("# design is:\n")
print(col_data)

counts <- counts[, rownames(col_data)]
cat('# got ', dim(counts)[2], ' samples\n')
cat('# building data set\n')
dxd <- DEXSeqDataSet(counts, col_data, design= ~ sample + exon + condition:exon,
                           featureID = mcols(gffs)$exon_id, 
                           groupID = mcols(gffs)$gene_id,
                           featureRanges=gffs)
BPPARAM = MulticoreParam(workers=6)
#cat('#subsetting for only retained introns\n')
#dxd <- dxd[grepl("I",featureIDs(dxd)),]
cat('# computing results\n')
dxd_results <- DEXSeq(dxd, BPPARAM=BPPARAM, quiet=F)
cat('#Subset results to only introns')
dxd_results <- dxd_results[grepl("I[0-9]+", dxd_results$featureID),]
dxd_results$padj <- p.adjust(dxd_results$pvalue, method="BH")
sig_exons <- rownames(subset(dxd_results, padj < 0.05))
cat('# got', length(sig_exons), ' significant exons\n')
sig_granges <- rowRanges(dxd[sig_exons,])
dxd_df <-as.data.frame(dxd_results)
cat('# outputting results')
write.table(dxd_df[,1:15],args$outfiles[1],
            quote=F, sep = "\t", row.names=FALSE)
export(sig_granges, args$outfiles[2])
save.image(args$outfiles[3])
Experiment.stop()
