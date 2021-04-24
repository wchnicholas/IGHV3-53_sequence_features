#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(data.table)
library(GGally)
library(e1071)
library(ggforce)
library(ggbeeswarm)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_F58Y_dist <- function(dist_table, graphname){
  print (dist_table)
  data_summary <- ddply(dist_table,c("resi58"), summarise, mean = mean(dist), sd = sd(dist))
  print (data_summary)
  textsize <- 9
  dodge_value <- 1
  p <-  ggplot() +
          geom_beeswarm(data=dist_table, aes(x=resi58, y=dist),
                        dodge.width=dodge_value,cex=8,size=0.5,color="gray40") + 
          geom_errorbar(data=data_summary, mapping=aes(x=resi58, ymin=mean-sd, ymax=mean+sd),
                        size=0.3, color="black", width=0.3, position=position_dodge(dodge_value)) +
          geom_point(data=data_summary, mapping = aes(x=resi58, y=mean),
                     size=4, color="black", shape="-", position=position_dodge(dodge_value)) +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(angle=0,hjust=1,vjust=0.5,colour = 'black'),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.title.x=element_blank(),
                legend.key.size=unit(2.2,"mm"),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize-1,face="bold"),
                legend.position='right') +
          scale_color_manual(values=c('black','black'),drop=FALSE) +
          labs(x=expression(bold(V['H']~residue~58)),y=expression(bold("Distance (A)")))
  ggsave(graphname,p,width=2,height=2,dpi=600)
  }

dist_table <- read_tsv('result/Y58F_Tstacking_dist.tsv')
plot_F58Y_dist(dist_table, 'graph/Dist_resi58.png')
