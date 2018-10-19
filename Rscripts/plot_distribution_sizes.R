#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)< 2) {
  stop("At least 2 arguments must be supplied", call.=FALSE) #input file and output directory
}

ggplot_plot_distr_sizes <- function()
{


  yy <- read.table(args[1],h=T)
 ggplot(data=yy, aes(x=factor(yy$type), y=yy$size, fill=yy$type))  + geom_boxplot() +
#~   ggplot(data=yy, aes(x=yy$size, fill=yy$type)) + geom_histogram(alpha=.5, position="identity") + 
#~   ggplot(data=yy, aes(x=yy$size, fill=yy$type)) + geom_histogram(alpha=.5, position="identity") + 
  labs(title = "Corrected reads/sequence sizes\n", x = "Reads or sequences", y = "Sizes distribution") +
  guides(fill=guide_legend(title="Reads type")) + theme_bw()

}

outName <- paste(args[2], "/plot_size_distribution.png", sep="")
ggsave(
  outName,
  ggplot_plot_distr_sizes(),
  width = 8,
  height = 8,
  dpi = 250
)




