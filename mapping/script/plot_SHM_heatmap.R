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

plot_SHM_heatmap <- function(SHM_data, graphname){
  textsize    <- 6
  p <-  ggplot(data=SHM_data, aes(x=variable, y=Ab, fill=value)) +
          geom_tile() +
          theme_cowplot(12) +
          scale_fill_gradientn(colours=c("white", "red"),
                         values=rescale(c(0, 1)),
                         guide="colorbar",
                         na.value="grey") +
          theme(plot.title=element_blank(),
                axis.text.x=element_text(size=textsize,angle=90,hjust=0.5,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(size=textsize,angle=0,hjust=1,vjust=0.5,colour = 'black'),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.title.x=element_text(size=textsize,face="bold"),
                strip.text = element_text(margin = margin(5)),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.key.size=unit(2.2,"mm"),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize-1,face="bold"),
                legend.position='none') +
          #scale_x_discrete(position = "top") +
          xlab("SHM") +
          ylab("Antibody")
  ggsave(graphname,p,width=8,height=10,dpi=600)
  }

plot_SHM_freq <- function(SHM_data_freq, graphname){
  textsize    <- 6
  colorscale  <- brewer.pal(3,"Accent")
  p <-  ggplot(data=SHM_data_freq, aes(x=SHM,y=freq,group=CDRH3_length)) +
          geom_bar(stat="identity", aes(fill=CDRH3_length),position=position_dodge()) +
          theme_cowplot(12) +
          scale_fill_manual(values=colorscale,drop=FALSE) +
          theme(axis.text.x=element_text(size=textsize-1, angle=90,hjust=0.5,vjust=0.5,colour = 'black', face="bold"),
                axis.text.y=element_text(size=textsize, angle=0,hjust=1,vjust=0.5,colour = 'black', face="bold"),
                axis.title=element_text(size=textsize,face="bold"),
                legend.key.size=unit(2.2,"mm"),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification = "center",
                legend.position='top') + 
          ylab(bquote(bold('Frequency (%)'))) +
          xlab("SHM")
  ggsave(graphname,p,width=9,height=2,dpi=600)
  }

classifying_CDRH3 <- function(CDRH3_length){
  if (CDRH3_length > 15){return ('long')}
  else {return ('short')}
  }

plot_SHM_freq_by_CDRH3_length <- function(SHM_mut_freq_by_CDRH3_length, graphname){
  textsize    <- 6
  p <-  ggplot(data=SHM_mut_freq_by_CDRH3_length, aes(x=variable, y=CDRH3_length, fill=value)) +
          geom_tile() +
          theme_cowplot(12) +
          scale_fill_gradientn(colours=c("white", "red"),
                         values=rescale(c(0, 1)),
                         guide="colorbar",
                         na.value="grey") +
          theme(plot.title=element_blank(),
                axis.text.x=element_text(size=textsize,angle=90,hjust=0.5,vjust=0.5,colour='black',face="bold"),
                axis.text.y=element_text(size=textsize,angle=0,hjust=1,vjust=0.5,colour='black',face="bold"),
                axis.title.y=element_text(size=textsize+1,face="bold"),
                axis.title.x=element_text(size=textsize+1,face="bold"),
                #strip.text = element_text(margin = margin(5)),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.key.size=unit(2.2,"mm"),
                legend.title=element_blank(),
                legend.text=element_text(size=textsize-1,face="bold"),
                legend.position='none') +
          #scale_x_discrete(position = "top") +
          xlab("Somatic mutation") +
          ylab("CDR H3 length")
  ggsave(graphname,p,width=6.5,height=1.8,dpi=600)
  }

SHM_data <- read_tsv("result/AntibodyDataSM_heatmap.tsv") %>%
              arrange(CDRH3_length, LC) %>%
              mutate(length_class=mapply(classifying_CDRH3, CDRH3_length))
Ab_level <- rev(SHM_data$Ab)
SHM_data_melt <- SHM_data  %>%
                   select(-length_class) %>%
                   melt(id.vars=c('Ab','LC','CDRH3_length')) %>%
                   mutate(Ab=factor(Ab,levels=Ab_level))
plot_SHM_heatmap(SHM_data_melt, 'graph/SHM_heatmap.png')

SHM_data_short <- SHM_data %>% filter(length_class == 'short') %>%
                    select(-Ab, -CDRH3_length, -LC, -length_class) %>%
                    colSums() %>%
                    melt() %>%
                    mutate(SHM=rownames(.)) %>%
                    select(SHM, value) %>%
                    data.table() %>%
                    rename(freq=value) %>%
                    mutate(freq=freq/length(filter(SHM_data, CDRH3_length=='short'))*100) %>%
                    mutate(CDRH3_length = 'HCDR3 length < 15') %>%
                    arrange(freq)
print (filter(SHM_data_short,SHM=='S53P'))
SHM_data_long  <- SHM_data %>% filter(length_class == 'long') %>%
                    select(-Ab, -CDRH3_length, -LC, -length_class) %>%
                    colSums() %>%
                    melt() %>%
                    mutate(SHM=rownames(.)) %>%
                    select(SHM, value) %>%
                    data.table() %>%
                    rename(freq=value) %>%
                    mutate(freq=freq/length(filter(SHM_data, CDRH3_length=='long'))*100) %>%
                    mutate(CDRH3_length = 'HCDR3 length ≥ 15')
SHM_levels <- rev(SHM_data_short$SHM)
SHM_data_freq  <- rbind(SHM_data_short, SHM_data_long) %>%
                    mutate(SHM=factor(SHM, levels=(SHM_levels)))
plot_SHM_freq(SHM_data_freq, 'graph/SHM_freq.png')

SHM_mut_freq_by_CDRH3_length <- SHM_data %>% 
                                  select(-Ab, -LC, -length_class) %>% 
                                  aggregate(list(.$CDRH3_length), mean) %>%
                                  melt(id="CDRH3_length") %>%
                                  filter(variable != 'Group.1') %>%
                                  mutate(CDRH3_length=factor(CDRH3_length,levels=sort(unique(CDRH3_length))))
SHM_mut_freq_dcast <- SHM_mut_freq_by_CDRH3_length %>%
                        dcast(variable~CDRH3_length) 
SHM_mut_freq_dcast[, "max"] <- apply(SHM_mut_freq_dcast[, 2:14], 1, max)
SHM_high_freq <- SHM_mut_freq_dcast %>%
                   filter(max > 0.05) %>%
                   .$variable
SHM_mut_freq_by_CDRH3_length <- SHM_mut_freq_by_CDRH3_length %>%
                                  filter(variable %in% SHM_high_freq)
plot_SHM_freq_by_CDRH3_length(SHM_mut_freq_by_CDRH3_length, 'graph/SHM_freq_by_CDRH3_length.png')

