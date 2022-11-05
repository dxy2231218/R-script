rm(list = ls())
library(pacman)
p_load(survival,survminer,finalfit)
dir.create('临床分组生存分析')
##临床分组生存分析
clini <- read.table('clinical.txt',sep='\t',header = T)
##gender和stage无差异，删除这两列
clini %<>% select(-c(gender,stage))
##整理数据
clini %<>% mutate(futime=futime/365,age=case_when(age>60 ~'>60',age<=60 ~'<=60',TRUE ~'unknow'),
                  grade=case_when(grade %in% c('G1','G2') ~'G1/G2',grade %in% c('G3','G4') ~'G3/G4',TRUE ~'unknow'),
                  T=case_when(str_detect(T,'^(T1|T2)') ~'T1/T2',str_detect(T,'^(T3|T4)') ~'T3/T4',TRUE~'unknow'),
                  M=ifelse(M %in% c('M0','M1'),M,'unknow'),
                  N=ifelse(N %in% c('N0','N1'),N,'unknow'))
rownames(clini) <- clini$Id
clini %<>% select(-1)

##读入风险比例数据
risk <- read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names = F) %>% column_to_rownames('id')
risk %<>% select(-c(1,2))
##合并数据集
common <- intersect(rownames(clini),rownames(risk))
data<- cbind(clini[common,],risk[common,])

do_surve<-function(data,title){
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
  p$plot<-p$plot+annotate('text',x = 4,y=15,label=bquote(log-rank~italic(pvalue)~.(pval)),size=4)
  p$plot<-p$plot+annotate('text',x=4,y=5,label=bquote(HR~value~.(HR)),size=4)
  p$table<-p$table+ theme(axis.ticks = element_blank(), axis.line = element_blank(), # 去除轴和刻度
                          axis.title.x = element_blank(), axis.text.x = element_blank(), # 去除标题和刻度值
                          axis.title.y = element_blank()) #, axis.text.y = element_blank() 与 tables.y.text = T 搭配使用
  #pdf('CliGroupSurv/TCGA_train_risk.pdf',height = 6,width = 6)
  pdf(paste0('临床分组生存分析/',title,'_risk.pdf'),height = 6,width = 6)
  print(p)
  dev.off()
}
cli <- c('grade','T','M','N')
for(col in cli){
  tab<-table(data[,col])
  for(name in names(tab)){
    if(name != 'unknown'){
      data_for_plot <- data[data[,col]==name,] %>% mutate(Risk=factor(Risk,levels = c('low','high')))
      name<-str_sub(name,1,2)
      do_surve(data_for_plot,name)
    }
  }
}
data_for_plot <- data %>% filter(age=='>60') %>% mutate(Risk=factor(Risk,levels = c('low','high')))
do_surve(data_for_plot,'age_over60.pdf')
