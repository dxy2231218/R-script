rm(list = ls())
dir.create('显著基因热图和森林图')
##构建单因素显著基因森林图和热图
##读入数据
uniSigExp <- read.table('模型/uni.SigExp.txt',sep='\t',header = T,check.names = F)
lncRNAExp <- read.table('lncRNA_logcpm.csv',sep=',',header = T,check.names = F,row.names = 1)
clini <- read.table('clinical.txt',sep='\t',header = T,check.names = F)
clini %<>% select(1:3) %>% mutate(futime=futime/365)
table(clini$Id) %>% table()
rownames(clini)<- clini$Id
clini %<>% select(-1)
##整理数据
data_for_heatmap <- lncRNAExp[colnames(uniSigExp)[-c(1,2,3)],]
#data_for_heatmap <- log2(data_for_heatmap+1)
lncRNA_exp %<>% select_if(str_detect(colnames(lncRNA_exp),'01A$'))
lncRNA_exp %<>% rename_all(~str_sub(.,1,12)) %>% t() %>% as.data.frame()
common <- intersect(rownames(clini),rownames(lncRNA_exp))
data <- cbind(clini[common,],lncRNA_exp[common,])

##绘制热图
type <- c(rep('N',3),rep('T',306))
names(type)<-colnames(data_for_heatmap)
type %<>% as.data.frame()
names(type)<-'type'
pdf('显著基因热图和森林图/SigGenesHeatmap.pdf',height = 10,width = 14)
pheatmap::pheatmap(data_for_heatmap, 
                   annotation=type, 
                   color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
                   cluster_cols =F,
                   cluster_rows = T,
                   scale="row",
                   show_colnames = F,
                   show_rownames = T,
                   fontsize = 8,
                   fontsize_row=8,
                   fontsize_col=8)
dev.off()
##绘制森林图
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(uniSigExp[,4:ncol(uniSigExp)])){
  cox <- coxph(Surv(futime, fustat) ~ uniSigExp[,i], data = uniSigExp)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.05){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
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
pdf(file="显著基因热图和森林图/forest.pdf", width = 6,height = 4.5)
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
