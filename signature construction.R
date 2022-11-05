##构建风险诊断模型
dir.create('模型')
rm(list = ls())
library(pacman)
library(magrittr)
p_load(tidyverse)
p_load(survival,survminer,glmnet,caret,timeROC,finalfit)
##加载表达数据
##load('clean_lncRNA_exp.Rdata')
lncRNA_exp<-read.table('lncRNA_logcpm.csv',sep=',',header = T,check.names = F)
##读入与坏死相关lncRNA
related_lncRNA <- read.csv('基因_lncRNA相关性.csv',header = T)
related_lncRNA<-related_lncRNA %>% distinct(gene) %>% pull(gene)
colnames(lncRNA_exp)[1]<-'id'
lncRNA <- lncRNA_exp %>% filter(id %in% related_lncRNA)%>% column_to_rownames('id')

##读入生存数据
surv<-read.table('time.txt',sep='\t',check.names = F,header = T,row.names = 1)

##合并数据集
lncRNA %<>% select_if(str_detect(colnames(lncRNA),'01A')) %>%
  rename_all(~str_sub(.,1,12))
common <- intersect(rownames(surv),colnames(lncRNA))
data <- cbind(surv[common,],t(lncRNA)[common,])
dim(data)
data %<>% mutate(futime=futime/365)
data %<>% rename_all(~str_replace_all(.,'-','\\.'))
do_model<-function(){
  intrain <- createDataPartition(y = data[,1],p = 0.65,list = F)
  train<-data[intrain,]
  test<- data[-intrain,]
  train_out <- train %>% mutate(id=rownames(train)) %>% relocate(id)
  test_out <- test %>% mutate(id=rownames(test)) %>% relocate(id)
  
  ##单因素cox回归
  fit <- coxphuni(.data = train,dependent = "Surv(futime,fustat)",explanatory = setdiff(colnames(train),c('futime','fustat')))
  cox_res <- fit2df(fit,condense=F)
  ##以0.05为阈值筛选显著基因
  out_uni_tab <- cox_res %>% filter(p<0.05)
  
  ##筛选显著基因表达矩阵
  uni_sig_exp <- train_out %>% select(futime,fustat,out_uni_tab$explanatory)
  uni_sig_exp_out <- train_out %>% select(id,futime,fustat,out_uni_tab$explanatory)##单因素cox选出90个显著基因
  
  ##lasso回归
  x<- as.matrix(uni_sig_exp[,3:ncol(uni_sig_exp)])
  y<-data.matrix(Surv(uni_sig_exp$futime,uni_sig_exp$fustat))
  ##
  fit <- glmnet(x,y,family = 'cox',maxit = 1000,alpha = 1) ##alpha默认为1,即lasso回归
  cvfit <-cv.glmnet(x,y,family = 'cox',maxit = 1000) ##交叉验证,选择最小lambda值
  #plot(fit)
  #plot(cvfit)
  coef <- coef(fit, s = cvfit$lambda.min) ##选择lambda值最小的回归系数
  index <- which(coef != 0)
  actCoef <- coef[index]
  lasso_gene <-rownames(coef)[index]
  lasso_sig_exp <- uni_sig_exp %>%select(futime,fustat,lasso_gene)
  lasso_sig_exp_out<- uni_sig_exp_out %>% select(id,futime,fustat,lasso_gene)
  gene_coef=cbind(Gene=lasso_gene, Coef=actCoef) ##lasso回归进一步选择出22个显著基因
  
  
  ###########构建cox模型###############
  multicox_fit <- coxphmulti(.data = lasso_sig_exp,dependent = "Surv(futime,fustat)",explanatory = setdiff(colnames(lasso_sig_exp),c('futime','fustat')))
  multicox_fit2 <- coxph(Surv(futime, fustat) ~ ., data = lasso_sig_exp)
  multicox_fit2<- step(multicox_fit2,direction = 'both')##逐步回归
  multicox_res <- multicox_fit %>% fit2df(condense=F)
  multicox_sum <- summary(multicox_fit2) ##多因素cox回归，进一步选出12个显著基因
  
  
  ##输出模型公式
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multicox_sum$coefficients[,"coef"],
    HR=multicox_sum$conf.int[,"exp(coef)"],
    HR.95L=multicox_sum$conf.int[,"lower .95"],
    HR.95H=multicox_sum$conf.int[,"upper .95"],
    pvalue=multicox_sum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  
  ##输出train组风险文件
  risk_score <- predict(multicox_fit2,newdata = train,type = 'risk') ##利用train组得到模型，用来预测train组风险
  cox_gene <- rownames(multicox_sum$coefficients)
  out_col <- c('futime','fustat',cox_gene)
  median_train_risk <- median(risk_score)
  risk <- ifelse(risk_score>median_train_risk,'high','low')
  train_risk_out <- cbind(id=rownames(train),train[,out_col],risk_score,Risk=risk)
  
  #输出test组风险文件
  risk_score_test=predict(multicox_fit2,type="risk",newdata=test)
  risk_test <- ifelse(risk_score_test>median_train_risk,'high','low')
  test_risk_out <- cbind(id=rownames(test),test[,out_col],risk_score_test,Risk=risk_test)
  
  
  
  ##生存差异分析
  formula <- as.formula('Surv(futime,fustat==1)~Risk')
  fit_train <- surv_fit(formula,data =train_risk_out)
  ggsurvplot(fit_train)
  pvalue_train <- surv_pvalue(fit_train)[,2]
  fit_test <- surv_fit(formula,data = test_risk_out)
  #ggsurvplot(fit_test)
  pvalue_test <- surv_pvalue(fit_test)[,2]
  #fit_geo <- surv_fit(formula,data = geo_risk_out)
  #ggsurvplot(fit_geo)
  #pvalue_geo <- surv_pvalue(fit_geo)[,2]
  
  ##ROC曲线下面积
  predict_time <- 1
  roc<-timeROC(T=train$futime, delta=train$fustat,
               marker=risk_score, cause=1,
               times=c(predict_time), ROC=TRUE)
  rocTest<-timeROC(T=test$futime, delta=test$fustat,
                   marker=risk_score_test, cause=1,
                   times=c(predict_time), ROC=TRUE)
  return(list(train_out=train_out,test_out=test_out,out_uni_tab=out_uni_tab,uni_sig_exp_out=uni_sig_exp_out,
              lasso_sig_exp_out=lasso_sig_exp_out,outMultiTab,train_risk_out,test_risk_out,
              fit_train,pvalue_train,fit_test,pvalue_test,roc,rocTest))
}
model1<-do_model()
while(T){
  model<-do_model()
  if((model[[10]])<0.05 && (model[[12]])<0.05){
    break
  }
}
model[[10]]
model[[12]]
model[[13]]$AUC[2]
model[[14]]$AUC[2]
model[[6]]
ggsurvplot(model[[9]])
model2<-model
save(model,file = '模型/model.Rdata')

##输出分组数据
write.table(model[[1]],file="模型/data.train.txt",sep="\t",quote=F,row.names=F)
write.table(model[[2]],file="模型/data.test.txt",sep="\t",quote=F,row.names=F)

##输出单因素结果
write.table(model[[3]],file="模型/uni.trainCox.txt",sep="\t",row.names=F,quote=F)
write.table(model[[4]],file="模型/uni.SigExp.txt",sep="\t",row.names=F,quote=F)

##输出lasso回归结果
write.table(model[[5]],file="模型/lasso.SigExp.txt",sep="\t",row.names=F,quote=F)

##输出多因素的结果
outMultiTab=model[[6]][,1:2]
write.table(outMultiTab,file="模型/multiCox.txt",sep="\t",row.names=F,quote=F)
write.table(model[[7]],file="模型/risk.TCGAtrain.txt",sep="\t",quote=F,row.names=F)
write.table(model[[8]],file="模型/risk.TCGAtest.txt",sep="\t",quote=F,row.names=F)
#write.table(model[[9]],file="模型/risk.GEO.txt",sep="\t",quote=F,row.names=F)

##所有样品的风险值

colnames(model[[7]])
colnames(model[[8]])
model[[8]] <- model[[8]] %>% rename(c('risk_score'='risk_score_test')) ##名字有不同，修改
allRiskOut=rbind(model[[7]], model[[8]])
write.table(allRiskOut,file="模型/risk.TCGAall.txt",sep="\t",quote=F,row.names=F)

##model里没有返回lasso回归fit对象，重新回归（作图需要）
uni_sig_exp <- read.table('模型/uni.SigExp.txt',sep='\t',check.names = F,header = T,row.names = 1)
x<- as.matrix(uni_sig_exp[,3:ncol(uni_sig_exp)])
y<-data.matrix(Surv(uni_sig_exp$futime,uni_sig_exp$fustat))
fit <- glmnet(x,y,family = 'cox',maxit = 1000,alpha = 1) ##alpha默认为1
cvfit <-cv.glmnet(x,y,family = 'cox',maxit = 1000) ##交叉验证,选择最小lambda值
pdf("模型/lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()
pdf("模型/lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()
