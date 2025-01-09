#!/usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(ggpubr)



# Merged Cell barcodes
cell.counts.file <- "~/Documents/scarecrow_test/merged.csv"
data <- read.csv(cell.counts.file)
names(data) <- c("CellBarcode", "Count")
ggplot(data[which(data$Count>1),], aes(x=Count, y=after_stat(density))) + geom_histogram(bins=100) + geom_density() + 
  labs(subtitle=paste0("Unique cell barcodes: ", dim(data)[1], "\nmedian observation count: ", median(data$Count))) +
  xlab("Barcode observation count") + ylab("Density") + theme_classic()
data <- data[-grep("null", data$CellBarcode),]
ggplot(data, aes(x=Count, y=after_stat(density))) + geom_histogram(bins=100) + geom_density() + 
  labs(subtitle=paste0("Unique cell barcodes (excluding any nulls): ", dim(data)[1], "\nmedian observation count: ", median(data$Count))) +
  xlab("Barcode observation count") + ylab("Density") + theme_classic()

# Barcode combination cumulative distribution
ggplot(data, aes(x=Count)) + geom_step(stat="ecdf") + theme_classic()



# Merged individual barcodes
cell.counts.file <- "~/Documents/scarecrow_test/barcodes_merged.csv"
data <- read.csv(cell.counts.file, header=F)
names(data) <- c("Index", "Barcode", "Count")
data$Index <- factor(data$Index)
ggplot(data[which(data$Barcode!="null"),], aes(x=Count, y=after_stat(density), group=Index, fill=Index, color=Index)) + 
  facet_grid(.~Index, scales="free") + geom_density() + scale_fill_jco() + scale_color_jco() +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), n.breaks=3) +
  xlab("Barcode observation count") + ylab("Density") + theme_classic() + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5), legend.position = "none")



# Cell barcodes
cell.counts.file <- "~/Documents/scarecrow_test/cDNA.fq.cell.counts.csv"
data <- read.csv(cell.counts.file)
ggplot(data[which(data$Count>1),], aes(x=Count, y=after_stat(density))) + geom_histogram(bins=100) + geom_density() + 
  labs(subtitle=paste0("Unique cell barcodes: ", dim(data)[1], "\nmedian observation count: ", median(data$Count))) +
  xlab("Barcode observation count") + ylab("Density") + theme_classic()
data <- data[-grep("null", data$CellBarcode),]
ggplot(data[which(data$Count>1),], aes(x=Count, y=after_stat(density))) + geom_histogram(bins=100) + geom_density() + 
  labs(subtitle=paste0("Unique cell barcodes (excluding any nulls): ", dim(data)[1], "\nmedian observation count: ", median(data$Count))) +
  xlab("Barcode observation count") + ylab("Density") + theme_classic()


# Individual barcodes
barcode.counts.file <- "~/Documents/scarecrow_test/cDNA.fq.barcode.counts.csv"
data <- read.csv(barcode.counts.file)
data$Index <- factor(data$Index)
ggplot(data[which(data$Barcode!="null"),], aes(x=Count, y=after_stat(density), group=Index, fill=Index, color=Index)) + geom_density() +
  scale_fill_jco() + scale_color_jco() +
  xlab("Barcode observation count") + ylab("Density") + theme_classic()
