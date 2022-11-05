rm(list = ls())
##绘制DCA决策曲线
##安装必须R包
library(timeROC)
library(pacman)
library(magrittr)
dir.create('决策曲线')
install.packages('ggDCA')
library(ggDCA)
p_load(survival,survminer,finalfit,tidyverse)
data <- read.csv('列线图临床风险合并数据.csv',header = T,row.names = 1)
data2<- data %>% filter_all(~.!='unknow')
data2$Nomogram<-nomoRisk
data2$Risk<-rt$Risk
data2 %<>% mutate(age=ifelse(age>60,1,0),Nomogram=ifelse(Nomogram>median(Nomogram),1,0))
##对Risk,Nomogram,Age,Gender,Grade,Stage做单因素cox回归
formula <- lapply(c('Risk','Nomogram','age','grade','T','M','N'),function(m)as.formula(paste0('Surv(futime,fustat)~',m)))
cox<- lapply(formula,function(m)coxph(formula = m,data=data2))
names(cox) <-c('Risk','Nomogram','age','grade','T','M','N') 
d_train=dca(cox[[1]],cox[[2]],cox[[3]],cox[[4]],cox[[5]],cox[[6]],cox[[7]],times=1) ##time表示预测时间为一年
d_train2 <- as.data.frame(d_train) ##转换成dataframe绘图
d_train2 %>% ggplot(aes(thresholds,NB))+geom_line(aes(color=model),lty=1,size=1)+theme_classic()+
  scale_y_continuous(limits = c(-0.1,0.1),name=bquote('Net Benifit'))+xlab(bquote('Risk Threshhold'))+
  scale_color_discrete(labels=c('Risk','Nomogram','age','grade','T','M','N','ALL','None'),name='')+
  theme(axis.title = element_text(size=14),axis.text = element_text(size=12),legend.title = element_text(size = 14),
        legend.text = element_text(size=12))
ggsave('决策曲线/DAC.pdf',height = 6,width = 6)

##绘制临床数据的ROC曲线,注意这里要把字符串变量转为数值
##重新读入数据

data3=cbind(data, Nomogram=nomoRisk)
data3 %<>% mutate(grade=case_when(grade=='G1'~1,grade=='G2'~2,grade=='G3'~3,grade=='GX'~4),
                  T=case_when(T=='T1'~1,T=='T2'~2,T=='T3'~3,T=='T4'~4,T=='TX'~5),
                  N=case_when(N=='N0'~0,N=='N1'~1,N=='NX'~2),
                  M=case_when(M=='M0'~0,M=='M1'~1,M=='MX'~2)
                 )
i=3
bioCol=rainbow(ncol(data3)-2, s=0.9, v=0.9)
pdf(file="决策曲线/cliROC2.pdf", width=6, height=6)
ROC_rt=timeROC(T=data3$futime,
               delta=data3$fustat,
               marker=data3$risk_score, cause=1,
               weighting='aalen',
               times=c(1),ROC=TRUE)
plot(ROC_rt, time=1, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

for(i in 3:(ncol(data3)-1)){
  print(i)
  ROC_rt=timeROC(T=data3$futime,
                 delta=data3$fustat,
                 marker=data3[,i], cause=1,
                 weighting='aalen',
                 times=c(1),ROC=TRUE)
  
  plot(ROC_rt, time=1, col=bioCol[i-1], title=FALSE, lwd=2, add=TRUE)
  aucText=c(aucText, paste0(colnames(data3)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(data3)-1)])
dev.off()
aucText

ROC_rt=timeROC(T=data3$futime,
               delta=data3$fustat,
               marker=data3[,'Nomogram'], cause=1,
               weighting='aalen',
               times=c(1),ROC=TRUE)

