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
require(cowplot)

plot_CDRH3_vs_LC_usage <- function(Ab_table, graphname){
  textsize <- 7
  p <-  ggplot(data=Ab_table, aes(x=LC_germline, y=CDRH3_length, color=Reference)) +
          geom_point(aes(fill=count,size=count),color='black',pch=21) +
          scale_fill_gradientn(colours=c("white", 'green'),
                limits=c(0,max(Ab_table$count)),
                values=rescale(c(0,1)),
                guide="colorbar",
                na.value="black") +
          scale_size_continuous(range = c(0,3.5)) +
          theme_cowplot(12) +
          theme(text = element_text(size=textsize,face="bold"),
                legend.key.size = unit(0.5, 'lines'),
                axis.text.y=element_text(size=textsize,face="bold",color='black'),
                axis.text.x=element_text(size=textsize,face="bold",color='black',angle=90,hjust=0.5,vjust=0.5),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold",color='black'),
                axis.text.x.top=element_text(angle=90,hjust=1,vjust=0.5,size=textsize,face="bold",color='black')) +
          ylab("CDR H3 length") +
          xlab("")
  ggsave(graphname,p,width=4,height=2.5,dpi=600)
  }

Ab_table  <- read_tsv("result/LC_germline_vs_CDRH3_length.tsv")
LC_germline_levels  <- sort(unique(Ab_table$LC_germline))
Ab_table <- Ab_table %>%
              mutate(LC_germline=factor(LC_germline,level=LC_germline_levels))
plot_CDRH3_vs_LC_usage(Ab_table, "graph/LC_usage_vs_CDRH3.png")
