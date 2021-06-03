# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: annotate.R
# time: 2021-04-20
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(dplyr))

# parsing args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop(paste(
        "Usage: Rscript", "annotate.R",
        "genome.chrom_size", "txdb.sqlite", "tss.bed"
    ))
}

genome <- read.table(args[1])
txdb <- loadDb(args[2])

filter <- list(tx_chrom = genome$V1)
cols <- c("gene_id")
p <- promoters(
    transcripts(txdb, columns=cols, filter=filter, use.names=F), upstream=50, downstream=50
)
p <- as.data.frame(p)
names(p$gene_id) <- row.names(p)
genes <- unlist(p$gene_id)
p <- p[names(genes),]
p <- distinct(p, gene_id, .keep_all = T)
p <- p[, c("seqnames", "start", "end")]
write.table(
    p, args[3], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
