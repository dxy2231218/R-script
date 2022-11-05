##药物敏感预测
rm(list = ls())
setwd('~/R/project1/')
library(pacman)
p_load(oncoPredict,gtools)
##测试oncoPredict
th<-theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
#dir<-'../oncoPredict/DataFiles/Training Data/'
#GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
#GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
#GDSC2_Res <- exp(GDSC2_Res) 
#dim(GDSC2_Res)
#testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
#testExpr[1:4,1:4]  
#colnames(testExpr)=paste0('test',colnames(testExpr))
#dim(testExpr) 
#calcPhenotype(trainingExprData = GDSC2_Expr,
#              trainingPtype = GDSC2_Res,
#              testExprData = testExpr,
#              batchCorrect = 'eb',  #   "eb" for ComBat  
#              powerTransformPhenotype = TRUE,
#              removeLowVaryingGenes = 0.2,
#              minNumSamples = 10, 
#              printOutput = TRUE, 
#              removeLowVaringGenesFrom = 'rawData' )



##
dir.create('药物预测')
load('clean_rna_exp.Rdata')
risk <- read.table('模型/risk.TCGAall.txt',sep='\t',header = T,check.names = F,row.names = 1)

exp %<>% rename_all(~str_sub(.,1,12))
common <- intersect(rownames(risk),colnames(exp))
exp <- exp[,common]
exp <- as.matrix(exp)
calcPhenotype(trainingExprData = GDSC2_Expr,
                            trainingPtype = GDSC2_Res,
                            testExprData = exp,
                            batchCorrect = 'eb',  #   "eb" for ComBat  
                            powerTransformPhenotype = TRUE,
                            removeLowVaryingGenes = 0.2,
                            minNumSamples = 10, 
                            printOutput = TRUE, 
                            removeLowVaringGenesFrom = 'rawData' )
drug_pre <- read.csv('calcPhenotype_Output/DrugPredictions.csv',header = T,check.names = F,row.names = 1)              
drug_pre %<>% mutate(Risk=risk$Risk) %>% relocate(Risk)
str_subset(colnames(drug_pre),'Veli')
durp_response<-drug_pre %>% group_by(Risk) %>% summarise_all(~mean(.))
durp_response %>% t()

##使用GDSC1再次预测
calcPhenotype(trainingExprData = GDSC1_Expr,
              trainingPtype = drugs_Res,
              testExprData = exp,
              batchCorrect = 'eb',  #"eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

drugs <-c('Veliparib_1018','Navitoclax_1011','Rucaparib_1175','Axitinib_1021',
          'AZD8055_1059','Bryostatin 1_197','DMOG_165','CHIR-99021_154','CHIR-99021_1241',
          'BX795_1037','Embelin_172','Palbociclib_1054','Crizotinib_37')
drugs_Res <- GDSC1_Res[,drugs]

drug_pre <- read.csv('calcPhenotype_Output/DrugPredictions.csv',header = T,check.names = F,row.names = 1)
drug_pre <- cbind(Risk=risk$Risk,drug_pre)
drug_pre[,2:14]<-log(drug_pre[,2:14])
drug_pre %>% ggplot(aes(Risk,Rucaparib_1175))+geom_boxplot(aes(fill=Risk))+
  ggsignif::geom_signif(comparisons = list(c('low','high')))+
  labs(y='Rucaparib sensitivity IC(50)')+theme_bw()+
  theme(axis.title = element_text(size=14),axis.text = element_text(size=12),
        legend.title = element_text(size=14),legend.text = element_text(size=12))
ggsave('Rucaparib.pdf',height = 5,width = 5)

drugs <- str_split(drugs,'_',simplify = T)[,1]
for(i in 2:ncol(drug_pre)){
  drug_pre %>% ggplot(aes(Risk,drug_pre[,i]))+geom_boxplot(aes(fill=Risk))+
    ggsignif::geom_signif(comparisons = list(c('low','high')))+
    labs(y=paste0(drugs[i-1],' sensitivity IC(50)'))+theme_bw()+
    theme(axis.title = element_text(size=14),axis.text = element_text(size=12),
          legend.title = element_text(size=14),legend.text = element_text(size=12))
  ggsave(paste0(drugs[i-1],'.pdf'),height = 5,width = 5)
}
