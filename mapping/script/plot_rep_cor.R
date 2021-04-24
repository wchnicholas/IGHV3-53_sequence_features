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

Plot_RepvsRep <- function(p,xlab,ylab,title){
  textsize <- 9
  p <- p +
         geom_point(size=0.5,color='black') +
         xlab(xlab) +
		 ylab(ylab) +
         theme_cowplot(12) +
		 theme(plot.title=element_text(size=11,face="bold",hjust=0.5),
               axis.title=element_text(size=textsize,face="bold"),
			   axis.text=element_text(size=textsize,face="bold"),
			   legend.title=element_blank(),
			   legend.text=element_text(size=textsize,face="bold"),
			   legend.position='none') +
         ggtitle(title)
  return (p)
  }

plot_bind_vs_exp <- function(enrich_table, title, graphname){
  textsize <- 7
  p <- ggplot(enrich_table, aes(x=bind,y=exp,color=LC)) +
         geom_rect(data=NULL,aes(xmin=0,xmax=1.1,ymin=-0.36,ymax=0.3), color=NA, fill=alpha('grey85', 0.5)) +
         geom_point(size=0.5) +
         xlab(expression(bold(log['10']~enrichment~'(binding)'))) +
         ylab(expression(bold(log['10']~enrichment~'(expression)'))) +
         ylim(-0.36,0.3) +
         theme_cowplot(12) +
		 theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title=element_text(size=textsize,face="bold"),
			   axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
			   legend.title=element_blank(),
			   legend.text=element_text(size=textsize,face="bold"),
			   legend.position='bottom') +
         #scale_color_manual(values=c('blue','black'),drop=FALSE) +
         scale_color_manual(values=c('black','black'),drop=FALSE) +
         guides(color=guide_legend(nrow=2, override.aes = list(size=0.5))) +
         ggtitle(title)
  #ggsave(graphname,p,height=2,width=2)
  ggsave(graphname,p,height=2,width=1.7)
  }

plot_bind_compare_class <- function(p, data_summary, ylab, graphname){
  textsize <- 7
  p <- p +
         geom_errorbar(data=data_summary, mapping=aes(x=LC, ymin=mean-sd, ymax=mean+sd),
                        size=0.3, color="black", width=0.3, position=position_dodge(1)) +
         geom_point(data=data_summary, mapping = aes(x=LC, y=mean),
                    size=4, color="black", shape="-", position=position_dodge(1)) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=textsize,face="bold",hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold"),
               axis.title.x=element_blank(),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='bottom') +
         ylab(ylab) 
  ggsave(graphname,p,height=1.5,width=1.3)
  }

enrich_table <- read_csv("result/CDRH3_Lib_RFindex.csv") %>%
                  select(-X2) %>%
                  mutate(bind=rowMeans(cbind(`IGGK1-9_bind_rep1_RF`,`IGGK1-9_bind_rep2_RF`))) %>%
                  mutate(exp=rowMeans(cbind(`IGKV1-9_exp_rep1_RF`,`IGKV1-9_exp_rep2_RF`))) %>%
                  mutate(AA=mapply(function(seq){return (str_sub(seq,3,-1))}, AA))
IGKV19_CDRH3  <- read_tsv('data/IGKV1-9_CDRH3.txt')$CDRH3

LC_list <- c()
for (AA in enrich_table$AA){
  if (AA %in% IGKV19_CDRH3){
    LC_list <- c(LC_list, 'from IGKV1-9 antibody')
    }
  else {
    LC_list <- c(LC_list, 'from non-IGKV1-9 antibody')
    }
  }


enrich_table <- enrich_table %>%
                  mutate(LC=LC_list)

p_KV19_bind_cor  <- Plot_RepvsRep(ggplot(enrich_table,aes(x=`IGGK1-9_bind_rep1_RF`,y=`IGGK1-9_bind_rep2_RF`)),
                                  expression(bold(log['10']~enrichment~'(replicate 1)')), 
                                  expression(bold(log['10']~enrichment~'(replicate 2)')),
                                  'Binding')
p_KV19_exp_cor   <- Plot_RepvsRep(ggplot(enrich_table,aes(x=`IGKV1-9_exp_rep1_RF`,y=`IGKV1-9_exp_rep2_RF`)),
                                  expression(bold(log['10']~enrichment~'(replicate 1)')), 
                                  expression(bold(log['10']~enrichment~'(replicate 2)')),
                                  'Expression')
print ("Pearson correlation between binding replicates:")
print (cor(enrich_table$`IGGK1-9_bind_rep1_RF`,enrich_table$`IGGK1-9_bind_rep2_RF`))
print ("Pearson correlation between expression replicates:")
print (cor(enrich_table$`IGKV1-9_exp_rep1_RF`,enrich_table$`IGKV1-9_exp_rep2_RF`))
p <- grid.arrange(p_KV19_bind_cor,p_KV19_exp_cor,ncol=2)
ggsave('graph/rep_compare.png',p,height=2.5,width=5)

print ("Pearson correlation between expression and binding:")
print (cor(enrich_table$bind,enrich_table$exp))
#plot_bind_vs_exp(enrich_table, 'CDR H3 variants in B38', 'graph/CDRH3lib_exp_vs_bind.png')
plot_bind_vs_exp(enrich_table, 'HCDR3 variants in B38', 'graph/CDRH3lib_exp_vs_bind_v2.png')

CDRH3_nonKV19_enriched <- enrich_table %>%
                            filter(LC=='from non-IGKV1-9 antibody') %>%
                            filter(bind>log10(2))
print ("enriched CDRH3s that are from non-IGKV1-9")
cat(as.character(CDRH3_nonKV19_enriched$AA),sep="\n")

print ("Difference in binding between IGKV1-9 and non-IGKV1-9:")
print (t.test(filter(enrich_table,LC=="from IGKV1-9 antibody")$bind, filter(enrich_table,LC=="from non-IGKV1-9 antibody")$bind))
plot_bind_compare_class(ggplot()+geom_beeswarm(data=enrich_table, aes(x=LC,y=bind), dodge.width=1,cex=3,size=0.1,color="grey"),
                        ddply(enrich_table,c("LC"), summarise, mean = mean(bind), sd = sd(bind)),
                        expression(bold(log['10']~enrichment~'(binding)')),
                        'graph/compare_class_binding.png')
print ("Difference in expression between IGKV1-9 and non-IGKV1-9:")
print (t.test(filter(enrich_table,LC=="from IGKV1-9 antibody")$exp, filter(enrich_table,LC=="from non-IGKV1-9 antibody")$exp))
plot_bind_compare_class(ggplot()+geom_beeswarm(data=enrich_table, aes(x=LC,y=exp), dodge.width=1,cex=3,size=0.1,color="grey"),
                        ddply(enrich_table,c("LC"), summarise, mean = mean(exp), sd = sd(exp)),
                        expression(bold(log['10']~enrichment~'(expression)')),
                        'graph/compare_class_expression.png')
