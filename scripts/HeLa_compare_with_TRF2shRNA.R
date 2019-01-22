
library("DESeq2")
library("tximport")
library("EnsDb.Hsapiens.v86")
library("tidyverse")
library("ggrepel")
library("gridExtra")
library("genefilter")
library("geneplotter")
library("VennDiagram")
library("viridis")
library("biomaRt")

pdf("HeLa_overexpression_TRF2_shRNA_overlaps.pdf")

samples <- read_tsv("20190103_received_samples.txt")
samples$cell <- sapply(samples$LongSampleName, function(x) str_split(x, pattern="-")[[1]][[1]])
samples$sample_replicate_id <- paste0(samples$Sample, "_", samples$Replicate)
samples_HeLa <- filter(samples, cell == "HeLa")
files <- file.path("salmon_on_hg38_output_cdna_ncrna", samples_HeLa$sample_replicate_id, "quant.sf")
names(files) <- samples_HeLa$sample_replicate_id
tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples_HeLa,
                                   design = ~ Sample)
ddsTxi <- DESeq(ddsTxi)
results_TERRA <- results(ddsTxi, contrast = c("Sample", "TERRA", "EV"))
results_ARIA <- results(ddsTxi, contrast = c("Sample", "ARIA", "EV"))
results_HP <- results(ddsTxi, contrast = c("Sample", "HP", "EV"))
genes_TERRA <- rownames(results_TERRA[which(results_TERRA$padj < 0.01), ])
genes_ARIA <- rownames(results_ARIA[which(results_ARIA$padj < 0.01), ])
genes_HP <- rownames(results_HP[which(results_HP$padj < 0.01), ])

setwd("~/TRF2_siRNA")
samples <- read_tsv("sample_info.txt")
samples$condition <- paste(samples$siRNA, samples$TimePoint, sep = "_")
samples$condition <- relevel(factor(samples$condition), ref = "TRF2_0")
files <- file.path("salmon_on_hg19_output", samples$Filename, "quant.sf")
names(files) <- samples$Filename
tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
ddsTxi <- DESeq(ddsTxi)
# Test for genes differentially expressed between timePoints in the TRF2 shRNA data
res48h <- results(ddsTxi, contrast = c("condition", "TRF2_48", "TRF2_0"))
res96h <- results(ddsTxi, contrast = c("condition", "TRF2_96", "TRF2_0"))
# Test for genes differentially expressed betweent timePoints in the control shRNA data
res48hControl <- results(ddsTxi, contrast = c("condition", "control_48", "control_0"))
res96hControl <- results(ddsTxi, contrast = c("condition", "control_96", "control_0"))
# Explicitly test for genes with a low LFC in the control timecourse
res48hControl_lowLFC <- results(ddsTxi,
                                lfcThreshold = log2(1.25),
                                altHypothesis = "lessAbs",
                                contrast = c("condition", "control_0", "control_48"))
res96hControl_lowLFC <- results(ddsTxi,
                                lfcThreshold = log2(1.25),
                                altHypothesis = "lessAbs",
                                contrast = c("condition", "control_0", "control_96"))
selected_genes_48h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res48h[which(res48h$padj < 0.1),])))
selected_genes_96h <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res96h[which(res96h$padj < 0.1),])))
selected_genes_48h_control <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res48h[which(res48hControl$padj < 0.1),])))
selected_genes_96h_control <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                             filter = GeneIdFilter(rownames(
                                               res96h[which(res96hControl$padj < 0.1),])))
selected_genes_48h_control_low_LFC <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                                     filter = GeneIdFilter(rownames(
                                                       res48h[which(res48hControl_lowLFC$padj < 0.1),])))
selected_genes_96h_control_low_LFC <- GenomicFeatures::genes(EnsDb.Hsapiens.v86,
                                                     filter = GeneIdFilter(rownames(
                                                       res96h[which(res96hControl_lowLFC$padj < 0.1),])))
# Include only the genes that also have a low LFC in the control samples
selected_genes_48h <- selected_genes_48h[(names(selected_genes_48h) %in% names(selected_genes_48h_control_low_LFC))]
selected_genes_96h <- selected_genes_96h[(names(selected_genes_96h) %in% names(selected_genes_96h_control_low_LFC))]
selected_genes_48h
genes_TRF2_48h <- names(selected_genes_48h)
genes_TRF2_96h <- names(selected_genes_96h)

setwd("~/scratch/TALES_RNAseq_analysis/")


venn <- draw.quad.venn(area1 = length(genes_TERRA),
                       area2 = length(genes_ARIA),
                       area3 = length(genes_HP),
                       area4 = length(genes_TRF2_96h),
                       n12 = length(intersect(genes_TERRA,
                                              genes_ARIA)),
                       n13 = length(intersect(genes_TERRA,
                                              genes_HP)),
                       n14 = length(intersect(genes_TERRA,
                                              genes_TRF2_96h)),
                       n23 = length(intersect(genes_ARIA,
                                              genes_HP)),
                       n24 = length(intersect(genes_ARIA,
                                              genes_TRF2_96h)),
                       n34 = length(intersect(genes_HP,
                                              genes_TRF2_96h)),
                       n123 = length(Reduce(intersect,
                                            list(genes_TERRA,
                                                 genes_ARIA,
                                                 genes_HP))),
                       n124 = length(Reduce(intersect,
                                            list(genes_TERRA,
                                                 genes_ARIA,
                                                 genes_TRF2_96h))),
                       n134 = length(Reduce(intersect,
                                            list(genes_TERRA,
                                                 genes_HP,
                                                 genes_TRF2_96h))),
                       n234 = length(Reduce(intersect,
                                            list(genes_ARIA,
                                                 genes_HP,
                                                 genes_TRF2_96h))),
                       n1234 = length(Reduce(intersect,
                                             list(genes_TERRA,
                                                  genes_ARIA,
                                                  genes_HP,
                                                  genes_TRF2_96h))),
                       category = c("TERRA",
                                    "ARIA",
                                    "HP",
                                    "TRF2"),
                       col = c("red", "blue", "green", "black"))

print(venn)

dev.off()


