#R code
library(ggplot2)
library(dplyr)
library(reshape)
library(scales)
library(cowplot)
library(gridExtra)
library(grid)
pdf(NULL)

plotsensorgram <- function(data_sign,data_fit,sample,ymax, ymin){
  textsize <- 9
  #textsize <- 7
  p <- ggplot() +
        #geom_line(data=data_sign, aes(x=Time, y=Signal,group=SampleID),size=0.5,color='blue') +
        #geom_line(data=data_fit, aes(x=Time, y=Signal,group=SampleID),size=0.3,color='#CC6666') +
        geom_line(data=data_sign, aes(x=Time, y=Signal,group=SampleID),size=0.8,color='blue') +
        geom_line(data=data_fit, aes(x=Time, y=Signal,group=SampleID),size=0.5,color='#CC6666') +
        theme_cowplot(12) +
        theme(plot.title=element_text(size=textsize,face="bold", hjust=0.5),
              axis.title=element_text(size=textsize,face="bold"),
              axis.text=element_text(size=textsize,face="bold"),
              axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
              legend.title=element_blank(),
              legend.key.size = unit(0.5, 'lines'),
              legend.text=element_text(size=textsize,face="bold"),
              legend.position='none') + 
        ylim(ymin,ymax) +
        xlab(bquote(bold(Time~'(s)'))) +
        ylab(bquote(bold(Response~'(nm)'))) +
        ggtitle(sample)
  return (p)
  }

adjusttime <- function(t){
  t$Time <- t$Time-min(t$Time)
  return (t)
  }

adjustsignal <- function(t, factor){
  t$Signal <- t$Signal*factor
  return (t)
  }

plot_wrapper <- function(sample){
  folder <- "result"
  fit    <- adjusttime(read.table(paste(folder,'/',sample,'_fit_All.compile',sep=''),header=1))
  sign   <- adjusttime(read.table(paste(folder,'/',sample,'_sign_All.compile',sep=''),header=1))
  ymax <- 0.6
  ymin <- -0.01
  if (grepl('COV107-23-F58Y',sample)){ymax <- 0.1}
  if (grepl('COVD21-C8-F58Y',sample)){ymax <- 0.1}
  if (grepl('CC12.3-WT',sample)){ymax <- 0.8; ymin <- -0.1}
  if (grepl('CC12.3-F58Y',sample)){ymax <- 0.8; ymin <- -0.1}
  if (grepl('COVA2-20-WT',sample)){ymax <- 0.1}
  if (grepl('COVA2-20-Y58F',sample)){ymax <- 0.1}
  p <- plotsensorgram(sign,fit,sample,ymax, ymin)
  return (p)
  }

targets <- c('RBD')
#binders <- c('COV107-23', 'COVD21-C8', 'COV107-23-swap', 'COVD21-C8-swap')
binders <- c('COV107-23-F58Y','COVD21-C8-F58Y','CC12.3-WT','CC12.3-F58Y','COVA2-20-WT','COVA2-20-Y58F')
samples <- c()
for (target in targets){
  for (binder in binders){
    sample <- paste(target,'_',binder,sep='')
    samples <- c(samples, sample)
    }
  }
plots <- lapply(samples,plot_wrapper)
p <- grid.arrange(grobs=plots,ncol=length(binders)/3)
ggsave(filename='graph/Sensorgram_All.png',p,height=6,width=5,dpi=600)
