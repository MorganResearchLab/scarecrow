#!/usr/bin/env Rscript
library(ggplot2)
library(ggsci)

args = commandArgs(trailingOnly=TRUE)
indir <- args[1]
outfile <- args[2]

# Read barcode results files
data <- lapply(
  list.files(indir, full=T, pattern="*.csv"),
  function(x) {
    read.csv(x, sep="\t")
    }
  )
data <- do.call("rbind", data)

# Generate plot
data$barcode_whitelist <- factor(data$barcode_whitelist, 
                                 levels = unique(data$barcode_whitelist)[order(nchar(unique(data$barcode_whitelist)))])
p <- ggplot(data, aes(x=start, fill=read)) + geom_histogram(binwidth=1) + 
  facet_grid(barcode_whitelist~read, scales="free") +
  xlab("barcode start position") + ylab("count") + scale_fill_jco(guide = "none") +
  scale_x_continuous(n.breaks = 10)+
  theme_classic()
ggsave(p, file=outfile)

