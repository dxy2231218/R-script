##绘制单因素森林图
dir.create('临床特征森林图')
rm(list = ls())
data <- read.table('列线图临床风险合并数据.csv',sep=',',header = T,check.names = F)
data %<>% mutate(grade=case_when(grade=='G1'~1,grade=='G2'~2,grade=='G3'~3,grade=='G4'~4),
                 T=case_when(T=='T1'~1,T=='T2'~2,T=='T3'~3,T=='T4'~4),
                 N=ifelse(N=='N0',0,1))
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(data[,3:ncol(data)])){
  cox <- coxph(Surv(futime, fustat) ~ data[,i], data = data)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  sigGenes=c(sigGenes,i)
  outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
 
 write.table(outTab,file="显著基因热图和森林图/tcgaUniCox.txt",sep="\t",row.names=F,quote=F)
#write.table(data,file="显著基因热图和森林图/tcgaUniSigExp.txt",sep="\t",row.names=F,quote=F)


######绘制森林图######
#读取输入文件
outTab[,-1]<-sapply(outTab[,-1],as.numeric)
sapply(outTab,typeof)
gene <- outTab$id
hr <- sprintf("%.3f",outTab$HR)
hrLow  <- sprintf("%.3f",outTab$"HR.95L")
hrHigh <- sprintf("%.3f",outTab$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(outTab$pvalue<0.001, "<0.001", sprintf("%.3f", outTab$pvalue))

#输出图形
pdf(file="临床特征森林图/forest.pdf", width = 6,height = 4.5)
n <- nrow(outTab)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#绘制森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


fit_uni<-coxphuni(.data = data,dependent = 'Surv(futime,fustat)',explanatory = colnames(data)[-c(1,2)])
fit_df_uni<-fit2df(fit_multi,condense=F)

fit_multi<-coxphmulti(.data = data,dependent = 'Surv(futime,fustat)',explanatory = c('risk_score','T','N'))
fit_df<-fit2df(fit_multi,condense=F)

write.csv(fit_df_uni,'临床特征森林图//临床特征单因素风险.csv',row.names = F)
write.csv(fit_df,'临床特征森林图//临床特征多因素风险.csv',row.names = F)
