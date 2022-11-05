##根据模型预测的风险指数，进行生存差异分析
rm(list = ls())
dir.create('生存分析')
library(pacman)
p_load(survival,survminer,dplyr,finalfit)
library(magrittr)
##读入数据
tcga_train <- read.table('模型/risk.TCGAtrain.txt',sep='\t',header = T,check.names = F,row.names = 1)
tcga_train %<>% mutate(Risk=factor(Risk,levels = c('low','high'))) ##修改因子让高风险组比低风险组
tcga_test <-read.table('模型/risk.TCGAtest.txt',sep='\t',header = T,check.names = F,row.names = 1) 
tcga_test %<>% mutate(Risk=factor(Risk,levels = c('low','high')))
tcga_all<-read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names =F,row.names = 1)
tcga_all %<>% mutate(Risk=factor(Risk,levels = c('low','high')))
#geo<-read.table('model/risk.GEO.txt',sep='\t',header = T,check.names = F)
#geo %<>% mutate(Risk=factor(Risk,levels = c('low','high')))
##
do_surve<-function(dat){
  data<-get(dat)
  formula <- as.formula("Surv(futime,fustat)~Risk")
  fit <- surv_fit(formula,data = data)
  pval<-surv_pvalue(fit)[,2]
  if(pval<0.01){
    pval<-formatC(pval,digits = 2,format = 'E')
  }else pval<- as.character(signif(pval,3))
  ##cox回归，计算HR值
  dependent <- "Surv(futime,fustat==1)"
  explanory <- 'Risk'
  cox <- coxphuni(dependent = dependent,explanatory = explanory,.data = data)
  HR<-fit2df(cox,condense=F)[,2]
  HR<-signif(HR,3)
  ggsurvplot(fit)
  p = ggsurvplot(fit, data = data, xlab = "Time(Years)", ylab = "Survival Probability (%)",  # ylab = "Survival probability (%)"
                 size = 1.2, palette = c("navyblue", "firebrick1"), # 曲线设置
                 ylim = c(0, 101),
                 legend.title = 'Risk', legend.labs = c("low risk","high risk"), font.legend = 14, legend = c(0.8,0.8), # 图例设置
                 pval = F, # 用annotate设置文本
                 break.time.by = 2, # x轴刻度间隔
                 risk.table = T, risk.table.height = 0.2, risk.table.fontsize = 5, risk.table.col = "strata", tables.y.text = F, # 风险表设置
                 font.x = 16, font.y = 16, font.xtickslab = 13, font.ytickslab = 13,
                 fun = function(y) y*100)  # 纵坐标刻度设置为百分数
  p$plot<-p$plot+annotate('text',x = 3.5,y=15,label=bquote(log-rank~italic(pvalue)~.(pval)),size=4)
  p$plot<-p$plot+annotate('text',x=3,y=5,label=bquote(HR~value~.(HR)),size=4)
  p$table<-p$table+ theme(axis.ticks = element_blank(), axis.line = element_blank(), # 去除轴和刻度
                          axis.title.x = element_blank(), axis.text.x = element_blank(), # 去除标题和刻度值
                          axis.title.y = element_blank()) #, axis.text.y = element_blank() 与 tables.y.text = T 搭配使用
  #pdf('survival/TCGA_train_risk.pdf',height = 6,width = 6)
  pdf(paste0('生存分析/',dat,'_risk.pdf'),height = 6,width = 6)
  print(p)
  dev.off()
}
do_surve('tcga_all')
