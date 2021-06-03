# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: compare.R
# time: 2021-04-22
suppressPackageStartupMessages(library(DESeq2))

# parse args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop(paste(
        "Usage: Rscript", "compare.R", "plan.txt",
        "counts.txt", "result.txt"
    ))
}

# read data
metadata <- read.table(args[1])
counts <- read.table(args[2])

# check sample
if (any(rownames(metadata) != colnames(counts))) {
    stop("Samples in counts.txt did not match samples in metadata.txt")
}

# process data
col_data <- cbind(rownames(metadata), metadata)
col_data[] <- lapply(col_data, as.factor)
counts_data <- as.matrix(counts)
dds <- DESeqDataSetFromMatrix(
    countData = counts_data, colData = col_data, design = ~condition
)
dds <- DESeq(dds)
res <- results(dds, c("condition", "numerator", "denominator"))
res <- as.data.frame(res)
write.table(
    res, args[3],
    sep = "\t", quote = FALSE, col.names = NA
)