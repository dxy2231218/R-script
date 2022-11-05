rm(list = ls())
dir.create('ROC')
p_load(survival,survminer,timeROC)
##读入数据
tcga_train_risk <- read.table('模型/risk.TCGAtrain.txt',sep='\t',header = T,check.names = F,row.names = 1)
tcga_test_risk <- read.table('模型/risk.TCGAtest.txt',sep='\t',header = T,check.names = F,row.names = 1)
tcga_test_risk %<>% rename(c('risk_score'='risk_score_test'))
tcga_all_risk <- read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)
#geo_risk <- read.table('model/risk.GEO.txt',sep='\t',header = T,check.names = F)
do_ROC<-function(dat){
  data<-get(dat)
  ROC<- timeROC(T=data$futime,delta=data$fustat,
                marker=data$risk_score,cause=1,
                weighting='aalen',
                times=c(1,3,5),ROC=TRUE)
  #return(ROC)
  #pdf('ROC/ROC_tcga_train_risk.pdf',height = 5,width = 5)
  pdf(paste0('ROC/ROC_',dat,'.pdf'),height = 5,width = 5)
  plot(ROC,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC$AUC[3]))),
         col=c("green",'blue','red'),lwd=2,bty = 'n')
  dev.off()
}
do_ROC('tcga_train_risk')
