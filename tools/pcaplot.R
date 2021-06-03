# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: pcaplot.R
# time: 2021-04-20
suppressPackageStartupMessages(library(MatrixGenerics))
suppressPackageStartupMessages(library(ggplot2))

# parse args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop(paste(
        "Usage: Rscript", "pcaplot.R", "metadata.txt",
        "counts.txt", "plot.pdf"
    ))
}

# read data
metadata <- read.table(args[1])
counts <- read.table(args[2])

# select upper quantile variance peaks
rv <- rowVars(as.matrix(counts))
select <- rv > quantile(rv, 0.75)

# plot
pca <- prcomp(t(counts[select, ]), scale. = TRUE)
percent <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)
groups <- apply(metadata, 1, paste, collapse = ":")
d <- data.frame(
    PC1 = pca$x[, 1], PC2 = pca$x[, 2], Groups = groups
)
g <- ggplot(
    data = d, aes_string(x = "PC1", y = "PC2", color = "Groups")
) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percent[1], "% variance")) +
    ylab(paste0("PC2: ", percent[2], "% variance")) +
    coord_fixed()
ggsave(args[3], g)