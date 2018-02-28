#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE) #input file and output directory
}

ggplot_plot_recall_prec_correctrate <- function()
{


yy <- read.table(args[1],h=T)

  ggplot(data=yy, aes(x=factor(yy$metric), y=yy$score, fill=yy$metric))  + geom_boxplot()  + 
  labs(title = "Impact of correction method on base-wise quality of sequences\n", x = "", y = "Distribution of recalls, precisions and correct base rate values in all corrected reads") +
  theme_bw() +
  guides(fill=guide_legend(title="Metrics"))

}

outName <- paste(args[2], "/plot_recall_precision.png", sep="")
ggsave(
  outName,
  ggplot_plot_recall_prec_correctrate(),
  width = 8,
  height = 8,
  dpi = 250
)




