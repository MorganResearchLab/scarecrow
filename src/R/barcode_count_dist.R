#!/usr/bin/env Rscript
library(ggplot2)
library(ggsci)
library(ggpubr)

#args = commandArgs(trailingOnly=TRUE)
#indir <- args[1]
#outfile <- args[2]

args <- NULL
indir <- "~/Documents/scarecrow_test/results"

# Read barcode results files
data <- lapply(
  list.files(indir, full=T, pattern="*.csv"),
  function(x) {
    read.csv(x, sep="\t")
    }
  )
data <- do.call("rbind", data)

# Generate plot
#data[grep("BC", data$barcode_whitelist),]$barcode_whitelist <- "['BC']"
#data$barcode_whitelist <- factor(data$barcode_whitelist, 
#                                 levels = unique(data$barcode_whitelist)[order(nchar(unique(data$barcode_whitelist)))])

p1 <- ggplot(data[which(data$orientation=="forward"),], aes(x=start, fill=read)) + geom_histogram(binwidth=1) + 
  facet_grid(barcode_whitelist~read) +
  xlab("barcode start position") + ylab("count") + scale_fill_jco(guide = "none") +
  scale_x_continuous(n.breaks = 10)+
  theme_classic() + theme(strip.text.y = element_text(angle=0))
p2 <- ggplot(data[which(data$orientation=="reverse"),], aes(x=start, fill=read)) + geom_histogram(binwidth=1) + 
  facet_grid(barcode_whitelist~read) +
  xlab("barcode start position") + ylab("count") + scale_fill_jco(guide = "none") +
  scale_x_continuous(n.breaks = 10)+
  theme_classic() + theme(strip.text.y = element_text(angle=0))
ggarrange(p1,p2,ncol=2,labels=c("forward", "reverse")) %>% 
  ggsave(file="~/Documents/scarecrow_test/plot.png",
         height=7, width=11, units="in")



counts <- data.frame(do.call("rbind", lapply(1:100, function(i) 
  table(data[which(data$read == "read2" & data$orientation == "forward" & data$start == i),]$barcode_whitelist))))
counts$max <- apply(counts,1,max)

  


# Notes
# BC1 and BC2 are on read 1 (forward) - counts ~ 7500 on 'R3_v3' and 'v1'
# BC3 is also on read 1 but in reverse - counts ~ 7500 on '198' and '255'
# Evercode WT v2 indicates BC1, BC2 and BC3 are on read 2, and that read 1 contains cDNA:
# https://support.parsebiosciences.com/hc/en-us/articles/14846676930452-What-are-the-run-configuration-and-sequencing-requirements-for-WT-libraries
# P5 | BC4 | R1 | .... cDNA .... | BC1 | . | BC2 | . | BC3 | R2 | BC4 | P7
# BC4 is the UDI-WT which is used for demultiplexing libraries/sublibraries
# R1 and R2 are the TruSeq primers (should be consistent across reads)


