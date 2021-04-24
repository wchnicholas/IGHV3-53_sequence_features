#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape2)
library(stringr)
library(dplyr)
require(cowplot)

plot_CDRH3_len_dist <- function(CDRH3_table, graphname){
  textsize <- 7
  p <- ggplot(CDRH3_table,aes(x=CDRH3_length,y=count)) +
         geom_bar(stat="identity", width=1) +
         theme_cowplot(12) +
         theme(axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'), 
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
     scale_fill_manual(values=c('black'),drop=FALSE) +
     labs(y=expression(bold('count')),x=expression()) +
     xlim(3,24)
  ggsave(graphname, p, height=0.9, width=1.8)
  }

CDRH3_len_table <- read_tsv('result/CDRH3_length_dist.tsv')
plot_CDRH3_len_dist(filter(CDRH3_len_table,class=='enriched'), 'graph/CDRH3_length_dist_enriched.png')
plot_CDRH3_len_dist(filter(CDRH3_len_table,class=='nonenriched'), 'graph/CDRH3_length_dist_nonenriched.png')
