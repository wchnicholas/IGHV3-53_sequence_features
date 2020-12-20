#R code
library(ggplot2)
library(reshape)
library(scales)
library(cowplot)
library(gridExtra)
library(dplyr)
library(readr)
library(grid)
pdf(NULL)

plot_affinity_compare <- function(data_affinity, Fab_name){
  textsize <- 7
  p  <- ggplot(data_affinity, aes(x=Variant, y=log10(Affinity), group=1)) +
	  geom_point(stat='summary', fun.y=sum) +
          stat_summary(fun.y=sum, geom="line") +
	  theme_cowplot(12) +
	  theme(axis.title=element_text(size=textsize,face="bold"),
		    axis.text=element_text(size=textsize,face="bold"),
            plot.title = element_text(size=textsize,face="bold",hjust=0.5),
	  	    legend.title=element_blank(),
		    legend.key.size = unit(0.5, 'lines'),
		    legend.text=element_text(size=textsize,face="bold"),
		    legend.position='none') +
          geom_hline(yintercept=c(3), linetype="dotted") +
          scale_y_continuous(breaks=c(0,1,2,3),labels=c('1','10','100','1000'),limit=c(-0.4,3.4)) +
          ggtitle(Fab_name) +
          ylab(bquote(bold("K"[D]~"(nM)"))) +
          xlab('')
  ggsave(filename=paste('graph/Affinity_compare_',Fab_name,'.png',sep=''),p,height=1.5,width=1.2,dpi=600)
  }

data_affinity <- read_tsv('result/Fab_affinity.tsv') %>%
                   mutate(Variant=factor(Variant,levels=c('Y58','F58')))
Fabs <- unique(data_affinity$Fab)
for (Fab_name in Fabs){
  data_sample <- data_affinity %>%
                   filter(Fab==Fab_name)
  plot_affinity_compare(data_sample,Fab_name)
  }
