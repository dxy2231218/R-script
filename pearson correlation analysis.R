library(tidyverse)
library(magrittr)
lncRNA_exp <- read.csv('差异分析/lncRNA_logcpm.csv',header = T,check.names = F,row.names = 1)
load('clean_rna_exp.Rdata')
lncRNA_exp %<>% select_if(str_detect(colnames(lncRNA_exp),'01A$'))
genes <- read.table('gene.txt',sep='\t',header = T)[,1]
exp <- exp[genes,]
load('差异分析/DiffGenes.Rdata')
lncRNA_exp <- lncRNA_exp[dge_list,]
common <- intersect(colnames(exp),colnames(lncRNA_exp))
exp <- exp[,common] %>% t()
exp<-exp[,apply(exp,2,function(m)!all(is.na(m)))] 
lncRNA_exp <- lncRNA_exp[,common] %>% t()
cor_list<-list()
for(i in 1:ncol(lncRNA_exp)){
  l=c()
  for(j in 1:ncol(exp)){
    l=append(l,cor(lncRNA_exp[,i],exp[,j]))
  }
  cor_list[[i]]<-l
}
cor_dataframe <- as.data.frame(cor_list)
colnames(cor_dataframe)<-colnames(lncRNA_exp)
rownames(cor_dataframe)<- colnames(exp)
cor_test_list<-list()
for(i in 1:ncol(lncRNA_exp)){
  l=c()
  for(j in 1:ncol(exp)){
    l=append(l,cor.test(lncRNA_exp[,i],exp[,j])$p.value)
  }
  cor_test_list[[i]]<-l
}

cor_test_dataframe <- as.data.frame(cor_test_list)
colnames(cor_test_dataframe)<-colnames(lncRNA_exp)
rownames(cor_test_dataframe)<- colnames(exp)
cor_dataframe %<>% mutate(id=rownames(cor_dataframe))
cor_dataframe <- reshape2::melt(cor_dataframe,id='id')
cor_test_dataframe %<>% mutate(id=rownames(cor_test_dataframe))
cor_test_dataframe <- reshape2::melt(cor_test_dataframe)

data <- cbind(cor_dataframe,p.value=cor_test_dataframe$value)
data %<>% rename(c('gene'='variable','cor_value'='value'))
data_filter <- data %>% filter(cor_value>0.3,p.value<0.05)
data_filter %>% count(gene)
write.csv(data_filter,file = '差异分析/基因_lncRNA相关性.csv',quote=F,row.names=F)
library(ggalluvial)
data_filter %>% ggplot(aes(id,gene))+geom_flow()+geom_stratum()
