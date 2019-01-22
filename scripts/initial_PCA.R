
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

# Function to get more than the first two PCs
# Adapted from the DESeq2::plotPCA function
plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, pcs = c(1, 2))
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup,
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, pcs[1]], PC2 = pca$x[, pcs[2]], group = group,
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[pcs[1]:pcs[2]]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
        geom_point(size = 3) + xlab(paste0("PC", pcs[1], ": ", round(percentVar[pcs[1]] *
        100), "% variance")) + ylab(paste0("PC", pcs[2], ": ", round(percentVar[pcs[2]] *
        100), "% variance")) + coord_fixed() + scale_colour_viridis_d()
}

pdf("initial_PCA.pdf", width = 10)


samples <- read_tsv("20190103_received_samples.txt")
samples$cell <- sapply(samples$LongSampleName, function(x) str_split(x, pattern="-")[[1]][[1]])
samples$sample_replicate_id <- paste0(samples$Sample, "_", samples$Replicate)


files <- file.path("salmon_on_hg38_output", samples$sample_replicate_id, "quant.sf")
names(files) <- samples$sample_replicate_id

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Sample)

vsd <- vst(ddsTxi, blind=TRUE)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell, vsd$Sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- rev(viridis(255))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("cell", "Sample"))


samples_HeLa <- filter(samples, cell == "HeLa")

files <- file.path("salmon_on_hg38_output", samples_HeLa$sample_replicate_id, "quant.sf")
names(files) <- samples_HeLa$sample_replicate_id

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples_HeLa,
                                   design = ~ Sample)

vsd <- vst(ddsTxi, blind=TRUE)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell, vsd$Sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- rev(viridis(255))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("cell", "Sample"))

plotPCA(vsd, intgroup=c("cell", "Sample"), pcs = c(2,3))

plotPCA(vsd, intgroup=c("cell", "Sample"), pcs = c(3,4))


samples_U2OS <- filter(samples, cell == "U2OS")

files <- file.path("salmon_on_hg38_output", samples_U2OS$sample_replicate_id, "quant.sf")
names(files) <- samples_U2OS$sample_replicate_id

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples_U2OS,
                                   design = ~ Sample)

vsd <- vst(ddsTxi, blind=TRUE)

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell, vsd$Sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- rev(viridis(255))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("cell", "Sample"))

plotPCA(vsd, intgroup=c("cell", "Sample"), pcs = c(2,3))

plotPCA(vsd, intgroup=c("cell", "Sample"), pcs = c(3,4))



dev.off()

