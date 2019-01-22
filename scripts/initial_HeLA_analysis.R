
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

pdf("HeLa_rna_seq_analysis.pdf")

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

plotMA(results_TERRA, main = "TERRA vs EV")
plotMA(results_ARIA, main = "ARIA vs EV")
plotMA(results_HP, main = "HP vs EV")


result_ranks <- data_frame(gene = rownames(results_TERRA),
                           padj_TERRA = results_TERRA$padj,
                           ranks_TERRA = rank(results_TERRA$padj, na.last = TRUE, ties.method = "max"),
                           padj_ARIA = results_ARIA$padj,
                           ranks_ARIA = rank(results_ARIA$padj, na.last = TRUE, ties.method = "max"),
                           padj_HP = results_HP$padj,
			   ranks_HP = rank(results_HP$padj, na.last = TRUE, ties.method = "max"))

ggplot(result_ranks, aes(x = ranks_TERRA, y = ranks_ARIA)) + geom_point()
ggplot(result_ranks, aes(x = ranks_TERRA, y = ranks_HP)) + geom_point()

ggplot(result_ranks, aes(x = -log10(padj_TERRA), y = -log10(padj_ARIA))) + geom_point()
ggplot(result_ranks, aes(x = -log10(padj_TERRA), y = -log10(padj_HP))) + geom_point()

genes_TERRA <- rownames(results_TERRA[which(results_TERRA$padj < 0.01), ])
genes_ARIA <- rownames(results_ARIA[which(results_ARIA$padj < 0.01), ])
genes_HP <- rownames(results_HP[which(results_HP$padj < 0.01), ])

venn.plot <- draw.triple.venn(area1 = length(genes_TERRA),
                              area2 = length(genes_ARIA),
                              area3 = length(genes_HP),
                              n12 = length(intersect(genes_TERRA,
                                                     genes_ARIA)),
                              n23 = length(intersect(genes_ARIA,
                                                     genes_HP)),
                              n13 = length(intersect(genes_TERRA,
                                                     genes_HP)),
                              n123 = length(Reduce(intersect,
                                                   list(genes_TERRA,
                                                        genes_ARIA,
                                                        genes_HP))),
                              category = c("TERRA",
                                           "ARIA",
                                           "HP"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

grid.newpage()
grid.draw(venn.plot)


genes_up_TERRA <- rownames(results_TERRA[which(results_TERRA$padj < 0.01 & results_TERRA$log2FoldChange > 0), ])
genes_down_TERRA <- rownames(results_TERRA[which(results_TERRA$padj < 0.01 & results_TERRA$log2FoldChange < 0), ])

genes_up_ARIA <- rownames(results_ARIA[which(results_ARIA$padj < 0.01 & results_ARIA$log2FoldChange > 0), ])
genes_down_ARIA <- rownames(results_ARIA[which(results_ARIA$padj < 0.01 & results_ARIA$log2FoldChange < 0), ])

genes_up_HP <- rownames(results_HP[which(results_HP$padj < 0.01 & results_HP$log2FoldChange > 0), ])
genes_down_HP <- rownames(results_HP[which(results_HP$padj < 0.01 & results_HP$log2FoldChange < 0), ])


venn.plot <- draw.triple.venn(area1 = length(genes_up_TERRA),
                              area2 = length(genes_up_ARIA),
                              area3 = length(genes_up_HP),
                              n12 = length(intersect(genes_up_TERRA,
                                                     genes_up_ARIA)),
                              n23 = length(intersect(genes_up_ARIA,
                                                     genes_up_HP)),
                              n13 = length(intersect(genes_up_TERRA,
                                                     genes_up_HP)),
                              n123 = length(Reduce(intersect,
                                                   list(genes_up_TERRA,
                                                        genes_up_ARIA,
                                                        genes_up_HP))),
                              category = c("TERRA up",
                                           "ARIA up",
                                           "HP up"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

grid.newpage()
grid.draw(venn.plot)

venn.plot <- draw.triple.venn(area1 = length(genes_down_TERRA),
                              area2 = length(genes_down_ARIA),
                              area3 = length(genes_down_HP),
                              n12 = length(intersect(genes_down_TERRA,
                                                     genes_down_ARIA)),
                              n23 = length(intersect(genes_down_ARIA,
                                                     genes_down_HP)),
                              n13 = length(intersect(genes_down_TERRA,
                                                     genes_down_HP)),
                              n123 = length(Reduce(intersect,
                                                   list(genes_down_TERRA,
                                                        genes_down_ARIA,
                                                        genes_down_HP))),
                              category = c("TERRA down",
                                           "ARIA down",
                                           "HP down"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

grid.newpage()
grid.draw(venn.plot)

dev.off()

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

library("topGO")

library("gplots")

for (sample in c("TERRA", "ARIA", "HP")) {

genes_names <- getBM(attributes = c('external_gene_name',
                                'description',
				'ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand'
                                ),
                 filters = 'ensembl_gene_id',
                 values = get(paste0("genes_", sample)),
                 mart = ensembl)

genes_names <- cbind(genes_names,
                           as.data.frame(get(paste0("results_", sample))[genes_names$ensembl_gene_id,]))

genes_names <- genes_names %>% arrange(padj)

write_tsv(genes_names,
          paste0("HeLa_results/genes_", sample, ".tsv"))

# GO analysis

# Trying to more explicitly compare genes with similar distribution of expression
overallBaseMean <- as.matrix(get(paste0("results_", sample))[, "baseMean", drop = F])

sig_idx <- match(get(paste0("genes_", sample)), rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)

}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

multidensity( list(
  all= log2(get(paste0("results_", sample))[is.na(get(paste0("results_", sample))$padj) == FALSE ,"baseMean"]) ,
  foreground =log2(get(paste0("results_", sample))[get(paste0("genes_", sample)), "baseMean"]),
  background =log2(get(paste0("results_", sample))[backG, "baseMean"])),
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(get(paste0("genes_", sample)),  backG)
inSelection =  geneIDs %in% get(paste0("genes_", sample))
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tgd.BP <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
            annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.BP <- runTest(tgd.BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.BP <- runTest(tgd.BP, algorithm = "classic", statistic = "Fisher" )

allRes.BP <- GenTable(tgd.BP, Fisher.elim = resultTopGO.elim.BP,
                   Fisher.classic = resultTopGO.classic.BP,
                   orderBy = "Fisher.classic", topNodes = 75)

pdf(paste0("HeLa_results/", sample, "_top_BP_GO_results.pdf"), width = 12)
showSigOfNodes(tgd.BP,
               score(resultTopGO.classic.BP),
               firstSigNodes = sum(allRes.BP$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()


write.table(allRes.BP %>% filter(Fisher.classic < 0.01),
            file = paste0("HeLa_results/", sample, "top_BP_GO_results.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

tgd.MF <- new( "topGOdata", ontology="MF", allGenes = alg, nodeSize=5,
            annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.MF <- runTest(tgd.MF, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.MF <- runTest(tgd.MF, algorithm = "classic", statistic = "Fisher" )

allRes.MF <- GenTable(tgd.MF, Fisher.elim = resultTopGO.elim.MF,
                   Fisher.classic = resultTopGO.classic.MF,
                   orderBy = "Fisher.classic", topNodes = 75)

pdf(paste0("HeLa_results/", sample, "_top_MF_GO_results.pdf"), width = 12)
showSigOfNodes(tgd.MF,
               score(resultTopGO.classic.MF),
               firstSigNodes = sum(allRes.MF$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()


write.table(allRes.MF %>% filter(Fisher.classic < 0.01),
            file = paste0("HeLa_results/", sample, "top_MF_GO_results.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


}


