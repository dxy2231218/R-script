library(magrittr)
##免疫浸润与riskscore相关分析
dir.create('免疫与Risk相关')
rm(list = ls())
risk<-read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)
load('clean_rna_exp.Rdata')
exp %<>% rename_all(~str_sub(.,1,12))
exp <- exp[,rownames(risk)]
exp<-exp[rowMeans(exp)>1,]
write.csv(cbind(gene=rownames(exp),exp),'免疫与Risk相关/mrna_exp.csv',row.names = F)
##mrna_exp.csv载入网页工具TIMER2分析免疫细胞含量，分析结果是estimation_matrix.csv
##读入免疫数据
immue <- read.csv('estimation_matrix.csv',header = T,row.names = 1,check.names = F)


##筛选低风险的病例
risk_l <- risk %>% filter(Risk=='low')
immue_l <- immue %>% select(rownames(risk_l))
immue_T <- t(immue_l)
platform <- str_split(colnames(immue_T),'_',simplify = T)[,2]
cell_type <-str_split(colnames(immue_T),'_',simplify = T)[,1] 
corr<-c()
for(i in 1:ncol(immue_T)){
  corr<-append(corr,cor(immue_T[,i],risk_l$risk_score))
}
corr_data <- data.frame(platform=platform,celltype=cell_type,corr=corr)
corr_test<-list()
for(i in 1:ncol(immue_T)){
  corr_test[[i]]<-cor.test(immue_T[,i],risk_l$risk_score)
}
corr_test_p.value<-sapply(corr_test,function(m)m$p.value)
corr_data <- cbind(corr_data,p.value=corr_test_p.value)
corr_data2 <- corr_data %>% filter(p.value<0.05)
corr_data2 %>% ggplot(aes(celltype,corr))+geom_point(aes(color=platform))+coord_flip()
corr_data2<-corr_data2[,-5]
corr_data2 %<>% mutate(celltype2=str_c(celltype,platform,sep='_'))
corr_data2 %>% ggplot(aes(celltype2,corr))+geom_point(aes(color=platform),size=3)+
  labs(y='Correlation coefficient',x='')+
  scale_x_discrete(limits=corr_data2$celltype2)+scale_color_brewer(palette = 'Dark2',name='software')+theme_bw()+
  theme(legend.title = element_text(size=16),legend.text = element_text(size=14),axis.title.x = element_text(size=14),axis.text.x = element_text(size=14),axis.text.y = element_text(color=rep(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")[c(6,1,2,5,7,3,4)],c(2,2,5,4,10,4,4))))+coord_flip()
ggsave('免疫与Risk相关/相关性点图.pdf',width = 10,height = 10)
write.csv(corr_data,file = '风险指数与免疫细胞相关性.csv',row.names = F,quote = F)
##绘制相关系数超过-0.3免疫细胞相关性图
corr_data3 <- corr_data2 %>% filter(corr< -0.2)
data_for_plot <- immue_T %>%as.data.frame()%>% dplyr::select(corr_data3$celltype2) %>% cbind(.,risk_score=risk_l$risk_score)
ggplot(aes(data_for_plot[,10],data_for_plot[,1]),data = data_for_plot)+geom_point()+geom_smooth(method = 'lm')+
  annotate(x=1,y=0.25,geom = 'text',label='R=-0.318   p.value=1.19e-07')+theme_bw()+
  labs(x='risk score',y=colnames(data_for_plot)[1])+theme(axis.title = element_text(size=14),axis.text = element_text(size=14))
ggsave(paste0(colnames(data_for_plot)[1],'_risk.pdf'),height = 6,width = 6)

for(i in 1:9){
  R=as.character(signif(corr_data3[i,3],3))
  p.value=formatC(corr_data3[i,4],digits = 3)
  print(list(R,p.value))
  y=max(data_for_plot[,i])
  ggplot(aes(data_for_plot[,10],data_for_plot[,i]),data = data_for_plot)+geom_point()+geom_smooth(method = 'lm')+
    annotate(x=0.25,y=y,geom = 'text',label=paste0('R=',R,'\t','p.value=',p.value))+theme_bw()+
    labs(x='Risk score',y=colnames(data_for_plot)[i])+theme(axis.title = element_text(size=14),axis.text = element_text(size=14))
  ggsave(paste0('免疫与Risk相关/',colnames(data_for_plot)[i],'_risk.pdf'),height = 6,width = 6)
}
paste0('R=',as.character(signif(corr_data3[1,3],3)),'\t','p.value=',formatC(corr_data3[1,4],digits = 3))

##
library(estimate)
exp %<>% mutate(GeneSymbol=rownames(exp)) 
exp %<>% relocate(GeneSymbol)
write.table(exp,file = '免疫与Risk相关/estimate_input.txt',quote = F,sep = '\t')
filterCommonGenes(input.f = '免疫与Risk相关/estimate_input.txt',output.f = 'estimate_group.gct',id = 'GeneSymbol')
estimateScore(input.ds = 'estimate_group.gct',output.ds = 'estimateScore.gct',platform = 'illumina')
scores<-read.table('estimateScore.gct',sep='\t',header = T,skip = 2)
rownames(scores)<-scores$NAME
scores <- t(scores[,4:ncol(scores)])
row_names <- rownames(scores) %>% str_replace_all('\\.','\\-')
rownames(scores)<- row_names
scores<- cbind(scores[,-4],risk_score=risk$Risk)
scores[,1:3] %<>% sapply(function(m)as.double(m))
scores %<>% as.data.frame() 
scores[,1:3]<-scores[,1:3] %>% sapply(function(m)as.double(m))
scores %>% sapply(typeof)
scores %>% ggplot(aes(risk_score,ImmuneScore))+geom_boxplot(aes(fill=risk_score))+
  ggsignif::geom_signif(comparisons = list(c('high','low')))+theme_bw()+
  labs(x='Risk score')+scale_fill_discrete(name='')+theme(axis.title = element_text(size=14),
                                                          axis.text=element_text(size=12),legend.position = 'none')
ggsave('免疫与Risk相关/ImmueScore_risk.pdf',height = 5,width = 5)

scores %>% ggplot(aes(risk_score,ESTIMATEScore))+geom_boxplot(aes(fill=risk_score))+
  ggsignif::geom_signif(comparisons = list(c('high','low')))+theme_bw()+
  labs(x='Risk score')+scale_fill_discrete(name='')+theme(axis.title = element_text(size=14),
                                                          axis.text=element_text(size=12),legend.position = 'none')
ggsave('免疫与Risk相关/ESTIMATEScore_risk.pdf',height = 5,width = 5)

scores %>% ggplot(aes(risk_score,StromalScore))+geom_boxplot(aes(fill=risk_score))+
  ggsignif::geom_signif(comparisons = list(c('high','low')))+theme_bw()+
  labs(x='Risk score')+scale_fill_discrete(name='')+theme(axis.title = element_text(size=14),
                                                          axis.text=element_text(size=12),legend.position = 'none')
ggsave('免疫与Risk相关/StromalScore_risk.pdf',height = 5,width = 5)

##得分与riskScore相关性
scores %<>% dplyr::rename('Risk'='risk_score')
scores %<>% mutate(risk_score=risk$risk_score)
corr<-c()
corr_test<-list()
for(i in 1:3){
  corr<-append(corr,cor(scores[,i],scores[,5]))
  corr_test[[i]]<-cor.test(scores[,i],scores[,5])
}
corr_test_p.value <- sapply(corr_test,function(m)m$p.value)
scores %>% ggplot(aes(risk_score,ImmuneScore))+geom_point()+geom_smooth(method = 'lm')

for(i in 1:3){
  R=as.character(signif(corr[i],3))
  p.value=formatC(corr_test_p.value[i],digits = 2)
  y=max(scores[,i])-500
  ggplot(aes(scores[,5],scores[,i]),data=scores)+geom_point()+geom_smooth(method = 'lm')+
    labs(x='Risk score',y=colnames(scores)[i])+theme_bw()+
    annotate(geom = 'text',x=6,y=y,label=paste0('R=',R,'\t','p.value=',p.value))+
    theme(axis.text = element_text(size=12),axis.title = element_text(size=14))
  ggsave(paste0('免疫与Risk相关/corr_',colnames(scores)[i],'_riskScore.pdf'),height = 5,width = 5)  
}
write.csv(scores,'ESTIMATE得分.csv',row.names = T)
