#!/usr/bin/env Rscript
library(dplyr)
library(SingleCellExperiment)

# Functions involved in plotting
library(ggplot2)
library(ggpointdensity)
library(ggdist)
library(ggpubr)
library(viridis)
library(scales)

args <- commandArgs(trailingOnly=TRUE)
print(args)

load(args[1])

# Filter count matrices to retain cells and features with at least 100 counts
cat("\t..filtering matrices to retain cells and genes with >99 counts\n")
mats_filtered <- lapply(mats[c(1,length(mats))], function(dat) {
  dat <- dat$mat
  gene_totals <- Matrix::rowSums(counts(dat))
  cell_totals <- Matrix::colSums(counts(dat))
  return(dat[gene_totals > 99, cell_totals > 99])
})

# Common cells and genes
cat("\t..identifying common barcodes between matrices\n")
barcodes <- colnames(mats_filtered[[1]])[colnames(mats_filtered[[1]]) %in% colnames(mats_filtered[[2]])]
genes <- rownames(mats_filtered[[1]])[rownames(mats_filtered[[1]]) %in% rownames(mats_filtered[[2]])]
mats_filtered[[1]] <- mats_filtered[[1]][genes, barcodes]
mats_filtered[[2]] <- mats_filtered[[2]][genes, barcodes]

knn_results <- lapply(c(2:10,20,30,50), function(k) {
	cat("\t..running k =", k, "kNN clustering\n")
	lapply(mats_filtered, function(m) FNN::get.knn(as.matrix(t(counts(m))), k = k)) 
	})

save(knn_results, file=gsub(".RData", "_knn.RData", args[1]))
cat("RData saved to", gsub(".RData", "_knn.RData", args[1]), "\n")
