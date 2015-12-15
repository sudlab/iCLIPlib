library(edgeR)
library(experimentr)

args <- Experiment.start()

counts <- read.delim(args$infiles[1])

sample_data <- strsplit(colnames(counts)[-1], "_")
sample_data <- do.call(rbind, sample_data)
sample_data <- as.data.frame(sample_data)
colnames(sample_data) <- c("Condition", "Fraction", "Replicate")
sample_data$Group = factor(paste(sample_data$Condition, sample_data$Fraction, sep="."))
design <- model.matrix(~0+Group, data=sample_data) 
colnames(design) <- levels(sample_data$Group)

rownames(counts) = counts$Geneid
counts <- as.matrix(counts[,-1])

dge <- DGEList(counts)

dge <- calcNormFactors(dge)

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
my.contrasts = makeContrasts(
  Fraction.Change=(ALYREF.Nuclear - ALYREF.Cytoplasmic) - (Control.Nuclear - Control.Cytoplasmic), levels=design)

lrt <- glmLRT(fit, contrast=my.contrasts)

lrt$table$ALYREF.Fraction <- (lrt$coefficients[,"ALYREF.Nuclear"] - lrt$coefficients[,"ALYREF.Cytoplasmic"])/log(2)
lrt$table$Control.Fraction <- (lrt$coefficients[,"Control.Nuclear"] - lrt$coefficients[,"Control.Cytoplasmic"])/log(2)
lrt$table$padj = p.adjust(lrt$table$PValue, method = "BH")

lrt$table$gene_id = rownames(lrt$table)

write.table(lrt$table, args$outfiles[1], sep="\t", row.names = FALSE)
Experiment.stop()