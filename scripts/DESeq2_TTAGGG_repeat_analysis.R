
library("DESeq2")

library("tximport")
library("EnsDb.Hsapiens.v86")

library("tidyverse")
library("ggrepel")
library("gridExtra")

library("genefilter")
library("geneplotter")
library("VennDiagram")

library("Rsamtools")

library("biomaRt")

pdf("HeLa_TTAGGG_repeat_analysis.pdf")


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


Dna <- getGenomeFaFile(EnsDb.Hsapiens.v86)

genes_in_data <- GenomicFeatures::genes(EnsDb.Hsapiens.v86, filter = GeneIdFilter(rownames(results_TERRA)))

genes_in_data <- genes_in_data[seqnames(genes_in_data) %in% seqnames(seqinfo(Dna))]

gene_telreps <- data_frame(gene_ID = names(genes_in_data),
                           gene_width = width(genes_in_data))

gene_telreps$genes_TERRA <- gene_telreps$gene_ID %in% genes_TERRA
gene_telreps$genes_ARIA <- gene_telreps$gene_ID %in% genes_ARIA
gene_telreps$genes_HP <- gene_telreps$gene_ID %in% genes_HP


geneSeqs <- getSeq(Dna, genes_in_data)

for (number_reps in 1:5) {

  gene_telrep_matches <- vmatchPattern(paste(rep("TTAGGG", number_reps), collapse = ""), geneSeqs)

  gene_telrep_number_matches <- unlist(lapply(gene_telrep_matches, length))

  gene_telreps[[paste0("telrep_matches_", number_reps)]] <- gene_telrep_number_matches

}


#p <- ggplot(mapping = aes(x = gene_width, y = telrep_matches, colour = genes_TERRA)) +
#  geom_point(data = gene_telreps %>%
#                     gather(key = "number_reps", value = "telrep_matches",
#                            telrep_matches_1:telrep_matches_5) %>%
#                     filter(genes_TERRA == FALSE)) +
#  geom_point(data = gene_telreps %>%
#                     gather(key = "number_reps", value = "telrep_matches",
#                            telrep_matches_1:telrep_matches_5) %>%
#                     filter(genes_TERRA == TRUE)) +
#  facet_wrap(~ number_reps, scales = "free_y")
#
#print(p)
#
#print(p + scale_x_log10())
#
#print(p + scale_x_log10() + scale_y_log10())
#
#p <- ggplot(gene_telreps %>%
#         gather(key = "number_reps", value = "telrep_matches",
#                telrep_matches_1:telrep_matches_5),
#       aes(x = gene_width, y = telrep_matches / gene_width)) +
#  geom_point() +
#  facet_wrap(~ number_reps, scales = "free_y")
#
#print(p)
#
#p <- ggplot(gene_telreps %>%
#         gather(key = "number_reps", value = "telrep_matches",
#                telrep_matches_1:telrep_matches_5),
#       aes(x = gene_width, y = telrep_matches / gene_width)) +
#  geom_point() +
#  scale_y_log10() +
#  facet_wrap(~ number_reps, scales = "free_y")
#
#print(p)
#
#p <- ggplot(mapping = aes(x = gene_width, y = telrep_matches / gene_width, colour = genes_TERRA)) +
#  geom_point(data = gene_telreps %>%
#               gather(key = "number_reps", value = "telrep_matches",
#                      telrep_matches_1:telrep_matches_5) %>% filter(genes_TERRA == FALSE)) +
#  geom_point(data = gene_telreps %>%
#               gather(key = "number_reps", value = "telrep_matches",
#                      telrep_matches_1:telrep_matches_5) %>% filter(genes_TERRA == TRUE)) +
#  scale_y_log10() +
#  facet_wrap(~ number_reps, scales = "free_y")
#
#print(p)
#
#p <- ggplot(mapping = aes(x = gene_width, y = telrep_matches / gene_width, colour = genes_TERRA)) +
#  geom_point(data = gene_telreps %>%
#               gather(key = "number_reps", value = "telrep_matches",
#                      telrep_matches_1:telrep_matches_5) %>% filter(genes_TERRA == FALSE)) +
#  geom_point(data = gene_telreps %>%
#               gather(key = "number_reps", value = "telrep_matches",
#                      telrep_matches_1:telrep_matches_5) %>% filter(genes_TERRA == TRUE)) +
#  scale_x_log10() +
#  scale_y_log10() +
#  facet_wrap(~ number_reps, scales = "free_y")
#
#print(p)
#

# Let's try to use linear modelling for the one repeat (the other don't show such a linear correlation)

mod <- lm(telrep_matches_1 ~ gene_width, data = gene_telreps)

gene_telreps <- modelr::add_predictions(gene_telreps, mod) %>%
                modelr::add_residuals(mod)


#p <- ggplot(gene_telreps, aes(x = gene_width)) +
#  geom_point(aes(y = telrep_matches_1)) +
#  geom_line(aes(y=pred), colour = "red", size = 1)
#
#print(p)
#
#p <- ggplot(gene_telreps, aes(x = gene_width)) +
#  geom_point(aes(y = resid))
#
#print(p)
#
#p <- ggplot(mapping = aes(x = gene_width)) +
#  geom_point(data = filter(gene_telreps, genes_TERRA == FALSE), mapping = aes(y = telrep_matches_1, colour = genes_TERRA)) +
#  geom_point(data = filter(gene_telreps, genes_TERRA == TRUE), mapping = aes(y = telrep_matches_1, colour = genes_TERRA)) +
#  geom_line(data = gene_telreps, mapping = aes(y=pred), colour = "red", size = 1)
#
#print(p)
#
#p <- ggplot(mapping = aes(x = gene_width)) +
#  geom_point(data = filter(gene_telreps, genes_TERRA == FALSE), mapping = aes(y = resid, colour = genes_TERRA)) +
#  geom_point(data = filter(gene_telreps, genes_TERRA == TRUE), mapping = aes(y = resid, colour = genes_TERRA))
#
#print(p)
#
#p <- ggplot(mapping = aes(x = gene_width)) +
#  geom_point(data = filter(gene_telreps, genes_TERRA == FALSE), mapping = aes(y = resid, colour = genes_TERRA)) +
#  geom_point(data = filter(gene_telreps, genes_TERRA == TRUE), mapping = aes(y = resid, colour = genes_TERRA)) +
#  scale_x_log10()
#
#print(p)
#

gathered_gene_telreps <- gene_telreps %>%
      gather(key = "ectopic_gene", value = "diff_expressed",
             starts_with("genes_")) %>%
      gather(key = "number_reps", value = "gene_telrep_matches",
      starts_with("telrep_matches"))


p <- ggplot(data = gathered_gene_telreps,
            mapping = aes(x = ectopic_gene,
                     y = gene_telrep_matches,
                     colour = diff_expressed)) +
     geom_boxplot() +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)


p <- ggplot(mapping = aes(x = ectopic_gene,
                     y = gene_telrep_matches,
                     colour = diff_expressed)) +
     geom_jitter(data = gathered_gene_telreps %>% filter(diff_expressed == FALSE)) +
     geom_jitter(data = gathered_gene_telreps %>% filter(diff_expressed == TRUE)) +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)


p <- ggplot(data = gathered_gene_telreps,
            mapping = aes(x = number_reps,
                     y = resid,
                     colour = diff_expressed)) +
     geom_boxplot() +
     facet_wrap(~ ectopic_gene) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

p <- ggplot(mapping = aes(x = number_reps,
                     y = resid,
                     colour = diff_expressed)) +
     geom_jitter(data = gathered_gene_telreps %>% filter(diff_expressed == FALSE)) +
     geom_jitter(data = gathered_gene_telreps %>% filter(diff_expressed == TRUE)) +
     facet_wrap(~ ectopic_gene) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)



# Now do this just for promoter sequences.

promoterSeqs <- getSeq(Dna, promoters(genes_in_data[seqnames(genes_in_data) %in% c(1:22, "X", "Y")],
                                      downstream = 1000))

promoter_telreps <- data_frame(gene_ID = names(genes_in_data[seqnames(genes_in_data) %in% c(1:22, "X", "Y")]),
                           gene_width = width(genes_in_data[seqnames(genes_in_data) %in% c(1:22, "X", "Y")]))

promoter_telreps$genes_TERRA <- promoter_telreps$gene_ID %in% genes_TERRA
promoter_telreps$genes_ARIA <- promoter_telreps$gene_ID %in% genes_ARIA
promoter_telreps$genes_HP <- promoter_telreps$gene_ID %in% genes_HP

for (number_reps in 1:5) {

  promoter_telrep_matches <- vmatchPattern(paste(rep("TTAGGG", number_reps), collapse = ""), promoterSeqs)

  promoter_telrep_number_matches <- unlist(lapply(promoter_telrep_matches, length))

  promoter_telreps[[paste0("promoter_telrep_matches_", number_reps)]] <- promoter_telrep_number_matches

}


#p <- ggplot(mapping = aes(x = gene_width, y = promoter_telrep_matches, colour = genes_TERRA)) +
#  geom_point(data = promoter_telreps %>%
#               gather(key = "number_reps", value = "promoter_telrep_matches",
#                      promoter_telrep_matches_1:promoter_telrep_matches_5) %>% filter(genes_TERRA == FALSE)) +
#  geom_point(data = promoter_telreps %>%
#               gather(key = "number_reps", value = "promoter_telrep_matches",
#                      promoter_telrep_matches_1:promoter_telrep_matches_5) %>% filter(genes_TERRA == TRUE)) +
#  facet_wrap(~ number_reps, scales = "free_y")
#
#print(p)
#
#print(p + scale_x_log10())
#

gathered_promoter_telreps <- promoter_telreps %>%
      gather(key = "ectopic_gene", value = "diff_expressed",
             starts_with("genes_")) %>%
      gather(key = "number_reps", value = "promoter_telrep_matches",
      starts_with("promoter_telrep"))

p <- ggplot(data = gathered_promoter_telreps,
            mapping = aes(x = ectopic_gene,
                     y = promoter_telrep_matches,
                     colour = diff_expressed)) +
     geom_boxplot() +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

p <- ggplot(mapping = aes(x = ectopic_gene,
                     y = promoter_telrep_matches,
                     colour = diff_expressed)) +
     geom_jitter(data = gathered_promoter_telreps %>% filter(diff_expressed == FALSE)) +
     geom_jitter(data = gathered_promoter_telreps %>% filter(diff_expressed == TRUE)) +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

# Let's see if there is any enrichment in exons.

exons_in_data <- GenomicFeatures::exons(EnsDb.Hsapiens.v86, filter = GeneIdFilter(rownames(results_TERRA)))

exons_in_data <- exons_in_data[seqnames(exons_in_data) %in% seqnames(seqinfo(Dna))]

exons_in_data <- lapply(unique(exons_in_data$gene_id), function(x) {

  gene_ranges <- GenomicRanges::reduce(exons_in_data[exons_in_data$gene_id == x])

  gene_ranges$gene_id <- x

  return(gene_ranges)

  })

names(exons_in_data) <- unique(exons_in_data$gene_id)

exon_sequences <- lapply(exons_in_data, function(x) getSeq(Dna, x))

exon_matches <- t(sapply(exon_sequences, function(x) {

     sapply(1:5, function(e) {
              telrep_matches <- vmatchPattern(paste(rep("TTAGGG", e), collapse = ""), x)
              sum(unlist(lapply(telrep_matches, length)))})

}))


rownames(exon_matches) <- unlist(lapply(exons_in_data, function(x) unique(x$gene_id)))

colnames(exon_matches) <- paste0("exon_telrep_matches_", 1:5)

reduced_exon_width <- unlist(lapply(exons_in_data, function(x) sum(width(x))))

exon_match_frame <- cbind(data_frame(gene_ID = rownames(exon_matches),
                                     exon_width = reduced_exon_width,
                                     genes_TERRA = rownames(exon_matches) %in% genes_TERRA,
                                     genes_ARIA = rownames(exon_matches) %in% genes_ARIA,
                                     genes_HP = rownames(exon_matches) %in% genes_HP),
                          exon_matches)


gathered_exon_telreps <- exon_match_frame %>%
      gather(key = "ectopic_gene", value = "diff_expressed",
             starts_with("genes_")) %>%
      gather(key = "number_reps", value = "exon_telrep_matches",
      starts_with("exon_telrep"))


p <- ggplot(data = gathered_exon_telreps,
            mapping = aes(x = ectopic_gene,
                     y = exon_telrep_matches,
                     colour = diff_expressed)) +
     geom_boxplot() +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)


# 5' and 3'UTR

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

results <- getBM(attributes = c('ensembl_gene_id',
                                'chromosome_name',
                                '5_utr_start', '5_utr_end',
                                '3_utr_start', '3_utr_end'),
                 filters = 'ensembl_gene_id',
                 values = rownames(results_TERRA),
                 mart = ensembl)


utr5_ranges <- with(results %>% filter(is.na(`5_utr_start`) == FALSE, is.na(`5_utr_end`) == FALSE),
                    GRanges(seqnames = chromosome_name,
                            ranges = IRanges(start = `5_utr_start`,
                                             end = `5_utr_end`),
                            strand = "*",
                            gene_id = ensembl_gene_id))

utr5_ranges <- utr5_ranges[seqnames(utr5_ranges) %in% seqnames(seqinfo(Dna))]

utr5_ranges_reduced <- lapply(unique(utr5_ranges$gene_id), function(x) {

  gene_ranges <- GenomicRanges::reduce(utr5_ranges[utr5_ranges$gene_id == x])

  gene_ranges$gene_id <- x

  return(gene_ranges)

})


names(utr5_ranges_reduced) <- unique(utr5_ranges$gene_id)

utr5_sequences <- lapply(utr5_ranges_reduced, function(x) getSeq(Dna, x))

utr5_matches <- t(sapply(utr5_sequences, function(x) {
     sapply(1:5, function(e) {
              telrep_matches <- vmatchPattern(paste(rep("TTAGGG", e), collapse = ""), x)
              sum(unlist(lapply(telrep_matches, length)))})
}))


rownames(utr5_matches) <- names(utr5_ranges_reduced) 

colnames(utr5_matches) <- paste0("utr5_telrep_matches_", 1:5)

reduced_utr5_width <- unlist(lapply(utr5_ranges_reduced, function(x) sum(width(x))))

utr5_match_frame <- cbind(data_frame(gene_ID = rownames(utr5_matches),
                                     utr5_width = reduced_utr5_width,
                                     genes_TERRA = rownames(utr5_matches) %in% genes_TERRA,
                                     genes_ARIA = rownames(utr5_matches) %in% genes_ARIA,
                                     genes_HP = rownames(utr5_matches) %in% genes_HP),
                          utr5_matches)


gathered_utr5_telreps <- utr5_match_frame %>%
      gather(key = "ectopic_gene", value = "diff_expressed",
             starts_with("genes_")) %>%
      gather(key = "number_reps", value = "utr5_telrep_matches",
      starts_with("utr5_telrep"))


p <- ggplot(data = gathered_utr5_telreps,
            mapping = aes(x = ectopic_gene,
                     y = utr5_telrep_matches,
                     colour = diff_expressed)) +
     geom_boxplot() +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)


utr3_ranges <- with(results %>% filter(is.na(`3_utr_start`) == FALSE, is.na(`3_utr_end`) == FALSE),
                    GRanges(seqnames = chromosome_name,
                            ranges = IRanges(start = `3_utr_start`,
                                             end = `3_utr_end`),
                            strand = "*",
                            gene_id = ensembl_gene_id))

utr3_ranges <- utr3_ranges[seqnames(utr3_ranges) %in% seqnames(seqinfo(Dna))]

utr3_ranges_reduced <- lapply(unique(utr3_ranges$gene_id), function(x) {

  gene_ranges <- GenomicRanges::reduce(utr3_ranges[utr3_ranges$gene_id == x])

  gene_ranges$gene_id <- x

  return(gene_ranges)

})


names(utr3_ranges_reduced) <- unique(utr3_ranges$gene_id)

utr3_sequences <- lapply(utr3_ranges_reduced, function(x) getSeq(Dna, x))

utr3_matches <- t(sapply(utr3_sequences, function(x) {
     sapply(1:5, function(e) {
              telrep_matches <- vmatchPattern(paste(rep("TTAGGG", e), collapse = ""), x)
              sum(unlist(lapply(telrep_matches, length)))})
}))


rownames(utr3_matches) <- names(utr3_ranges_reduced) 

colnames(utr3_matches) <- paste0("utr3_telrep_matches_", 1:5)

reduced_utr3_width <- unlist(lapply(utr3_ranges_reduced, function(x) sum(width(x))))

utr3_match_frame <- cbind(data_frame(gene_ID = rownames(utr3_matches),
                                     utr3_width = reduced_utr3_width,
                                     genes_TERRA = rownames(utr3_matches) %in% genes_TERRA,
                                     genes_ARIA = rownames(utr3_matches) %in% genes_ARIA,
                                     genes_HP = rownames(utr3_matches) %in% genes_HP),
                          utr3_matches)


gathered_utr3_telreps <- utr3_match_frame %>%
      gather(key = "ectopic_gene", value = "diff_expressed",
             starts_with("genes_")) %>%
      gather(key = "number_reps", value = "utr3_telrep_matches",
      starts_with("utr3_telrep"))


p <- ggplot(data = gathered_utr3_telreps,
            mapping = aes(x = ectopic_gene,
                     y = utr3_telrep_matches,
                     colour = diff_expressed)) +
     geom_boxplot() +
     facet_wrap(~ number_reps, scales = "free_y") +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)




dev.off()


