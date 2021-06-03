# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: annotate.R
# time: 2021-04-20
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(GenomicFeatures))

# parsing args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop(paste(
        "Usage: Rscript", "annotate.R",
        "peaks.bed", "txdb.sqlite", "annotated.bed"
    ))
}

peaks <- readPeakFile(args[1])
txdb <- loadDb(args[2])

annotation <- annotatePeak(peaks, TxDb = txdb, annoDb = "org.Hs.eg.db")
annotation <- as.data.frame(annotation)
columns <- c(
    "seqnames", "start", "end", "V4", "V5", "V6",
    "SYMBOL", "annotation", "distanceToTSS",
    "geneId", "ENSEMBL", "transcriptId",
    "geneStart", "geneEnd", "geneLength", "geneStrand"
)
annotation <- annotation[, columns]
annotation$geneStrand <- ifelse(annotation$geneStrand == 1, "+", "-")
write.table(
    annotation, args[3],
    quote = FALSE, sep = "\t", na = "", row.names = FALSE,
    col.names = c(
        "#chrom", "start", "end", "name", "score", "strand",
        "symbol", "annotation", "distance2tss",
        "entrez", "ensembl", "transcript",
        "gene_start", "gene_end", "gene_length", "gene_strand"
    )
)