library(ggplot2)
library(ggsci)
data <- read.csv("/Users/s14dw4/Documents/Repos/barcode_counts.csv", sep="\t")

ggplot(data, aes(x=start, fill=read)) + geom_histogram(binwidth=1) + facet_grid(read~.) +
  xlab("barcode start position") + ylab("count") + scale_fill_jco(guide = "none") +
  theme_classic()
