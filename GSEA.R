##基因富集分析
dir.create('GSEA')
rm(list = ls())
library(pacman)
setwd('R/project2/')
options(connectionObserver = NULL) ##防止加载org.Hs.eg.db出现错误
devtools::install_github('GuangchuangYu/clusterProfiler') ##重装clusterProfiler包
BiocManager::install('org.Hs.eg.db')
BiocManager::install('clusterProfiler')
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
p_load(tidyverse)

##读入mRNA表达数据
load('clean_rna_exp.Rdata')##加载的是FPKM数据，不是count值数据
table(str_sub(colnames(exp),1,12)) %>% table() ##判断列名截短后是否有重复
exp %<>% rename_all(~str_sub(.,1,12))
##读入风险数据
risk <- read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)

##高低风险比较得到logFC
common <- intersect(rownames(risk),colnames(exp))
risk <- risk[common,]
exp_data<-exp[,common]
high_risk_sample <- risk %>%filter(Risk=='high') %>% rownames()
low_risk_sample <- risk %>%filter(Risk=='low') %>% rownames()
exp_data_high <- exp_data[,high_risk_sample]
exp_data_low <- exp_data[,low_risk_sample]
mean_h <- rowMeans(exp_data_high)
mean_h[mean_h<0.00001]<-0.00001
mean_l <- rowMeans(exp_data_low)
mean_l[mean_l<0.00001]<-0.00001
logFC <- log2(mean_h)-log2(mean_l)
logFC<-sort(logFC,decreasing = T)
genes<- names(logFC)

##读入基因集文件件
gmt <- read.gmt('c2.cp.kegg.v7.4.symbols.gmt')

##基因集富集分析(GESA)，重要！！！
kk <- GSEA(logFC,TERM2GENE = gmt,pvalueCutoff = 1,eps = 0)
kktab <- as.data.frame(kk) %>% filter(pvalue<0.05)
write.table(kktab,file="GSEA/GSEA.result.txt",sep="\t",quote=F,row.names = F)
save(kk,file = 'GSEA/GSEA.result.Rdata')
load('GSEA.result.Rdata')
##输出高风险富集图形，只展示前5通路
kkup <- kktab %>% filter(NES>0)
show_term <-rownames(kkup)[1:5]
gseaplot2(kk, show_term, base_size=8, title="Enriched in high risk group")
ggsave('GSEA.highRisk.pdf', width=7, height=5.5)

##输出低风险富集图形，只展示前5通路
kkdown <- kktab %>% filter(NES<0)
show_term <-rownames(kkdown)[1:5]
gseaplot2(kk, show_term, base_size=8, title="Enriched in low risk group")
ggsave('GSEA.lowRisk.pdf', width=7, height=5.5)

kktab_sort<-kktab %>% dplyr::arrange(desc(abs(NES)))
