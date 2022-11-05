##绘制列线图
rm(list = ls())
dir.create('列线图')
library(rms) 
install.packages('regplot') ##绘制列线图好用的包
library(regplot)
##读入临床数据
clinic <- read.table('clinical.txt',sep='\t',check.names = F,header = T,row.names = 1)
##读入风险数据
risk <- read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)

##整理数据
common <- intersect(rownames(clinic),rownames(risk))
data <- risk[common,] %>% dplyr::select(futime,fustat,risk_score) %>% cbind(.,clinic[common,-c(1,2,4,6)])
table(data$N)
##去除unknow行
data %<>% mutate(
                 T=case_when(str_detect(T,'^T1')~'T1',str_detect(T,'T2')~'T2',
                             str_detect(T,'^T3')~'T3',str_detect(T,'^T4')~'T4',
                             T=='TX'~'TX',
                             TRUE~'unknow'),
                 M=case_when(M=='M0'~'M0',M=='M1'~'M1',M=='MX'~'MX',TRUE~'unknow'),
                 N=case_when(N=='N0'~'N0',N=='N1'~'N1',N=='NX'~'NX',TRUE~'unknow')
                 )
data %<>% filter_all(~. != 'unknow')
data %<>% filter(grade != 'G4')
##调整因子水平

data %<>%  mutate(grade=factor(grade,levels = c('G1','G2','G3','GX')),
                  T=factor(T,levels = c('T1','T2','T3','T4','TX')),
                  N=factor(N,levels = c('N0','N1','NX')),age=as.numeric(age),
                  M=factor(M,levels = c('M0','M1','MX')))
write.csv(data,file = '列线图临床风险合并数据.csv',row.names = F)

fit <- coxph(Surv(futime,fustat)~.,data=data)

nom1<-regplot(fit,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=data[1,],
              rank="sd",
              failtime = c(1,3,5),
              prfail = F)


#列线图风险打分
nomoRisk=predict(fit, data=data, type="risk")
rt=cbind(risk[rownames(data),], Nomogram=nomoRisk)
outTab=rbind(ID=colnames(rt), rt)
write.table(outTab, file="列线图/nomoRisk.txt", sep="\t", col.names=F, quote=F)

#校准曲线
pdf(file="列线图/calibration.pdf", width=5, height=5)
#1年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
#3年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=100)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
#5年校准曲线
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=100)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

##分开校准曲线

#1年校准曲线
pdf(file="列线图/calibration_1_year.pdf", width=5, height=5)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=1)
cal <- calibrate(f, cmethod="KM", method="boot", u=1, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
dev.off()
#3年校准曲线
pdf(file="列线图/calibration_3_year.pdf", width=5, height=5)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=3)
cal <- calibrate(f, cmethod="KM", method="boot", u=3, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="Nomogram-predicted OS (%)", ylab="bserved OS (%)", lwd=1.5, col="green", sub=F)
dev.off()
#5年校准曲线
pdf(file="列线图/calibration_5_year.pdf", width=5, height=5)
f <- cph(Surv(futime, fustat) ~ Nomogram, x=T, y=T, surv=T, data=rt, time.inc=5)
cal <- calibrate(f, cmethod="KM", method="boot", u=5, m=(nrow(rt)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="Nomogram-predicted OS (%)", ylab="bserved OS (%)",  lwd=1.5, col="green", sub=F)
dev.off()
