
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

pdf("U2OS_rna_seq_analysis.pdf")

samples <- read_tsv("20190103_received_samples.txt")
samples$cell <- sapply(samples$LongSampleName, function(x) str_split(x, pattern="-")[[1]][[1]])
samples$sample_replicate_id <- paste0(samples$Sample, "_", samples$Replicate)

samples_U2OS <- filter(samples, cell == "U2OS")

files <- file.path("salmon_on_hg38_output_cdna_ncrna", samples_U2OS$sample_replicate_id, "quant.sf")
names(files) <- samples_U2OS$sample_replicate_id

tx2gene <- values(transcripts(EnsDb.Hsapiens.v86))[, c("tx_id", "gene_id")]

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples_U2OS,
                                   design = ~ Sample)

ddsTxi <- DESeq(ddsTxi)


results_AD6 <- results(ddsTxi, contrast = c("Sample", "AD6", "AD6_NT"))
results_NLS3 <- results(ddsTxi, contrast = c("Sample", "N3", "N3_NT"))
results_SID4 <- results(ddsTxi, contrast = c("Sample", "S4", "S4_NT"))

plotMA(results_AD6, main = "AD6")
plotMA(results_NLS3, main = "NLS3")
plotMA(results_SID4, main = "SID4")

genes_all_AD6 <- rownames(results_AD6[which(results_AD6$padj < 0.01), ])
genes_all_NLS3 <- rownames(results_NLS3[which(results_NLS3$padj < 0.01), ])
genes_all_SID4 <- rownames(results_SID4[which(results_SID4$padj < 0.01), ])

grid.newpage()
venn.plot <- draw.triple.venn(area1 = length(genes_all_AD6),
                              area2 = length(genes_all_NLS3),
                              area3 = length(genes_all_SID4),
                              n12 = length(intersect(genes_all_AD6,
                                                     genes_all_NLS3)),
                              n23 = length(intersect(genes_all_NLS3,
                                                     genes_all_SID4)),
                              n13 = length(intersect(genes_all_AD6,
                                                     genes_all_SID4)),
                              n123 = length(Reduce(intersect,
                                                   list(genes_all_AD6,
                                                        genes_all_NLS3,
                                                        genes_all_SID4))),
                              category = c("AD6",
                                           "NLS3",
                                           "SID4"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

genes_up_AD6 <- rownames(results_AD6[which(results_AD6$padj < 0.01 & results_AD6$log2FoldChange > 0), ])
genes_down_AD6 <- rownames(results_AD6[which(results_AD6$padj < 0.01 & results_AD6$log2FoldChange < 0), ])

genes_up_NLS3 <- rownames(results_NLS3[which(results_NLS3$padj < 0.01 & results_NLS3$log2FoldChange > 0), ])
genes_down_NLS3 <- rownames(results_NLS3[which(results_NLS3$padj < 0.01 & results_NLS3$log2FoldChange < 0), ])

genes_up_SID4 <- rownames(results_SID4[which(results_SID4$padj < 0.01 & results_SID4$log2FoldChange > 0), ])
genes_down_SID4 <- rownames(results_SID4[which(results_SID4$padj < 0.01 & results_SID4$log2FoldChange < 0), ])

grid.newpage()
venn.plot <- draw.triple.venn(area1 = length(genes_up_AD6),
                              area2 = length(genes_up_NLS3),
                              area3 = length(genes_up_SID4),
                              n12 = length(intersect(genes_up_AD6,
                                                     genes_up_NLS3)),
                              n23 = length(intersect(genes_up_NLS3,
                                                     genes_up_SID4)),
                              n13 = length(intersect(genes_up_AD6,
                                                     genes_up_SID4)),
                              n123 = length(Reduce(intersect,
                                                   list(genes_up_AD6,
                                                        genes_up_NLS3,
                                                        genes_up_SID4))),
                              category = c("AD6 up",
                                           "NLS3 up",
                                           "SID4 up"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

grid.newpage()
venn.plot <- draw.triple.venn(area1 = length(genes_down_AD6),
                              area2 = length(genes_down_NLS3),
                              area3 = length(genes_down_SID4),
                              n12 = length(intersect(genes_down_AD6,
                                                     genes_down_NLS3)),
                              n23 = length(intersect(genes_down_NLS3,
                                                     genes_down_SID4)),
                              n13 = length(intersect(genes_down_AD6,
                                                     genes_down_SID4)),
                              n123 = length(Reduce(intersect,
                                                   list(genes_down_AD6,
                                                        genes_down_NLS3,
                                                        genes_down_SID4))),
                              category = c("AD6 down",
                                           "NLS3 down",
                                           "SID4 down"),
                              euler.d = TRUE,
                              overrideTriple = 1,
                              col = c("red", "blue", "green"))

#grid.newpage()
grid.draw(venn.plot)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes_names_AD6 <- getBM(attributes = c('external_gene_name',
                                'ensembl_gene_id'),
                 filters = 'ensembl_gene_id',
                 values = genes_all_AD6,
                 mart = ensembl)

genes_names_SID4 <- getBM(attributes = c('external_gene_name',
                                'ensembl_gene_id'),
                 filters = 'ensembl_gene_id',
                 values = genes_all_SID4,
                 mart = ensembl)

genes_names_AD6 <- as.data.frame(cbind(genes_names_AD6, results_AD6[genes_names_AD6$ensembl_gene_id,]))
genes_names_SID4 <- as.data.frame(cbind(genes_names_SID4, results_SID4[genes_names_SID4$ensembl_gene_id,]))


# make volcano plots
p <- ggplot(results_AD6 %>% as.data.frame(),
       aes(x = log2FoldChange, y = -log10(padj), colour = padj < 0.01)) +
       geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
       geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
       geom_point() +
       geom_label_repel(data = rbind(genes_names_AD6 %>% arrange(padj) %>% head(10),
                                     genes_names_AD6 %>% arrange(desc(abs(log2FoldChange))) %>% head(10)) %>% unique(),
                        mapping = aes(label = external_gene_name)) +
       ggtitle("AD6") +
       scale_x_continuous(limits = c(-10, 10)) +
       theme_bw()
print(p)


p <- ggplot(results_SID4 %>% as.data.frame(),
       aes(x = log2FoldChange, y = -log10(padj), colour = padj < 0.01)) +
       geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
       geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
       geom_point() +
       geom_label_repel(data = rbind(genes_names_SID4 %>% arrange(padj) %>% head(10),
                                     genes_names_SID4 %>% arrange(desc(abs(log2FoldChange))) %>% head(10)) %>% unique(),
                        mapping = aes(label = external_gene_name)) +
       ggtitle("SID4") +
       scale_x_continuous(limits = c(-10, 10)) +
       theme_bw()
print(p)



# Plot AD6 vs SID4 expressions using VST
vsd <- vst(ddsTxi, blind = FALSE) %>%
         assay() %>%
         data.frame(gene_ID = rownames(.),
                    .) %>%
         gather(key = "sample",
                value = "vst",
                AD6_1:S4_NT_3) %>%
         separate(sample,
                  into = c("sample", "treated"),
                  extra = "merge") %>%
         mutate(treated = ifelse(grepl("NT", treated), "Non-treated", "Treated"))

vsd_by_gene <- vsd %>%
                 group_by(treated, sample, gene_ID) %>%
                 summarise(mean_vst = mean(vst))


ggplot(vsd_by_gene %>% filter(treated == "Treated", sample %in% c("AD6", "S4")) %>% spread(sample, mean_vst),
       aes(x = AD6, y = S4)) +
  geom_point()

dev.off()

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

library("topGO")

library("gplots")

for (sample in c("AD6", "SID4")) {

for (direction in c("all", "up", "down")) {

genes_names <- getBM(attributes = c('external_gene_name',
                                'description',
                                'ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand'
                                ),
                 filters = 'ensembl_gene_id',
                 values = get(paste0("genes_", direction, "_", sample)),
                 mart = ensembl)

genes_names <- cbind(genes_names,
                           as.data.frame(get(paste0("results_", sample))[genes_names$ensembl_gene_id,]))

genes_names <- genes_names %>% arrange(padj)

write_tsv(genes_names,
          paste0("U2OS_results/genes_", sample, ".tsv"))

# GO analysis

# Trying to more explicitly compare genes with similar distribution of expression
overallBaseMean <- as.matrix(get(paste0("results_", sample))[, "baseMean", drop = F])

sig_idx <- match(get(paste0("genes_", direction, "_", sample)), rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)

}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

multidensity( list(
  all= log2(get(paste0("results_", sample))[is.na(get(paste0("results_", sample))$padj) == FALSE ,"baseMean"]) ,
  foreground =log2(get(paste0("results_", sample))[get(paste0("genes_", direction, "_", sample)), "baseMean"]),
  background =log2(get(paste0("results_", sample))[backG, "baseMean"])),
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(get(paste0("genes_", direction, "_", sample)),  backG)
inSelection =  geneIDs %in% get(paste0("genes_", direction, "_", sample))
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

pdf(paste0("U2OS_results/", sample, "_", direction, "_top_BP_GO_results.pdf"), width = 12)
showSigOfNodes(tgd.BP,
               score(resultTopGO.classic.BP),
               firstSigNodes = sum(allRes.BP$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()


write.table(allRes.BP %>% filter(Fisher.classic < 0.01),
            file = paste0("U2OS_results/", sample, "_", direction, "_top_BP_GO_results.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

tgd.MF <- new( "topGOdata", ontology="MF", allGenes = alg, nodeSize=5,
            annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )

## run tests
resultTopGO.elim.MF <- runTest(tgd.MF, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic.MF <- runTest(tgd.MF, algorithm = "classic", statistic = "Fisher" )

allRes.MF <- GenTable(tgd.MF, Fisher.elim = resultTopGO.elim.MF,
                   Fisher.classic = resultTopGO.classic.MF,
                   orderBy = "Fisher.classic", topNodes = 75)

pdf(paste0("U2OS_results/", sample, "_", direction, "_top_MF_GO_results.pdf"), width = 12)
showSigOfNodes(tgd.MF,
               score(resultTopGO.classic.MF),
               firstSigNodes = sum(allRes.MF$Fisher.classic < 0.01),
               useInfo = 'all')
dev.off()



write.table(allRes.MF %>% filter(Fisher.classic < 0.01),
            file = paste0("U2OS_results/", sample, "_", direction, "_top_MF_GO_results.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


}

}
