####GSE13850???ݼ?????
setwd("E:/FZZK-20104-10/00_rawdata/GSE13850")
rm(list = ls())
library(readxl)
library(tidyverse)
library(GEOquery)
library(data.table)
library(limma)
# ??ȡ series_matrix?ļ? ----
eset <- getGEO(filename = "01_GSE13850_series_matrix.txt.gz", getGPL = F)

probes_expr <- exprs(eset)
dim(probes_expr)
data_use <- as.data.frame(probes_expr)
boxplot(data_use,outline=F,notch=F,las=2)

####??ȡע???ļ?
anno <- fread('02_GPL96-57554.txt',data.table = FALSE)
anno <- anno[,c("ID","Gene Symbol")]
same <- match(rownames(data_use),anno$ID)
data_use$Gene_symbol <- anno[same,c("Gene Symbol")]
colnames(data_use)


data_use <- data_use[!is.na(data_use$Gene_symbol)&!grepl('///',data_use$Gene_symbol),]#ɾ????һ???ݶ?Ӧ??GeneID??????(ע??????///)
data_use <- data_use[data_use$Gene_symbol!="",]#ɾ??û??Gene_symbol??????
data_use<- as.data.frame(avereps(data_use[,-41],ID = data_use$Gene_symbol) )#??ͬ????ȡƽ??ֵ
data_use <- na.omit(data_use)#5932
View(data_use)

####??ȡ?ٴ???Ϣ
phenoDat <- pData(eset)
dim(phenoDat)
sample <- row.names(phenoDat)
View(phenoDat)
phenoDat <- cbind(sample,phenoDat)
colnames(phenoDat)
group <- phenoDat[,c("sample","title")]
View(group)
table(group$title)
group$group <- c(rep("Control",10),rep("OP",10),rep("Control",10),rep("OP",10))
group <- group[,-2]
group <- group[order(group$group,decreasing = F),]

data1 <- t(data_use)
sam <- group$sample
data1 <- data1[sam,]


write.csv(data1,file = "03_expr.csv",row.names = T,quote = F)
write.csv(group,file = "04_group.csv",row.names = F,quote = F)



####GSE56815???ݼ?????
setwd("E:/FZZK-20104-10/?½??ļ???/00_rawdata/GSE56815")
rm(list = ls())
library(readxl)
library(tidyverse)
library(GEOquery)
library(data.table)
library(limma)
# ??ȡ series_matrix?ļ? ----
eset <- getGEO(filename = "01_GSE56815_series_matrix.txt.gz", getGPL = F)

probes_expr <- exprs(eset)
dim(probes_expr)
data_use <- as.data.frame(probes_expr)
boxplot(data_use,outline=F,notch=F,las=2)

####??ȡע???ļ?
anno <- fread('02_GPL96-57554.txt',data.table = FALSE)
anno <- anno[,c("ID","Gene Symbol")]
same <- match(rownames(data_use),anno$ID)
data_use$Gene_symbol <- anno[same,c("Gene Symbol")]
colnames(data_use)


data_use <- data_use[!is.na(data_use$Gene_symbol)&!grepl('///',data_use$Gene_symbol),]#ɾ????һ???ݶ?Ӧ??GeneID??????(ע??????///)
data_use <- data_use[data_use$Gene_symbol!="",]#ɾ??û??Gene_symbol??????
data_use<- as.data.frame(avereps(data_use[,-81],ID = data_use$Gene_symbol) )#??ͬ????ȡƽ??ֵ
data_use <- na.omit(data_use)#5932
View(data_use)

####??ȡ?ٴ???Ϣ
phenoDat <- pData(eset)
dim(phenoDat)
sample <- row.names(phenoDat)
View(phenoDat)
phenoDat <- cbind(sample,phenoDat)
colnames(phenoDat)
group <- phenoDat[,c("sample","title")]
View(group)
table(group$title)
group$group <- c(rep("Control",40),rep("OP",40))
group <- group[,-2]
#group <- group[order(group$group,decreasing = F),]


write.csv(data_use,file = "03_expr.csv",row.names = T,quote = F)
write.csv(group,file = "04_group.csv",row.names = F,quote = F)


################################????????
setwd("E:/FZZK-20104-10/01_DEG")


library(edgeR)
library(limma)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(ggfun)
library(grid)
library(ggplot2)
library(data.table)

rm(list = ls())
expr <- fread("E:/FZZK-20104-10/00_rawdata/GSE56815/03_expr.csv")
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- expr[,-1]
#expr <- t(expr)
group <- read.csv("E:/FZZK-20104-10/00_rawdata/GSE56815/04_group.csv",header = T,row.names = 1,check.names = F)


design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr)
contrast.matrix <- makeContrasts(OP - Control,levels = design)
fit <- lmFit(expr,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG <- topTable(fit2,coef = 1,n = Inf)


#??ɽͼ####
DEGs_file <- DEG
DEGs_file$sig <- ifelse(DEGs_file$P.Value < 0.05 & abs(DEGs_file$logFC) >= 0,
                        ifelse(DEGs_file$logFC > 0,'up','down'),'none')
table(DEGs_file$sig)

volcano_color <- c( "#FF9999", "#99CC00", "gray80")
#?ѻ?ͼ??????װ?ɺ???
volcano_plot <- function(dataframe){
  ggplot(data = dataframe,aes(logFC,-log10(P.Value),color = sig))+geom_point(size = 2.2)+
    scale_colour_manual(values = volcano_color)+
    geom_point(data = top,aes(logFC,-log10(P.Value)),color =  "#6A5ACD",size=3,alpha=0.9)+
    labs(title="DEGs_volcano",x="log2(fold change)",y="-log10(P.Value)")+ 
    geom_hline(yintercept = c(-log10(0.05)),lwd= 1.3,color= "orange",lty= "dashed")+
    geom_vline(xintercept = c(-1,1),lwd= 1.5,color= "orange",lty= "dashed")+
    geom_label_repel(data= top,aes(logFC,-log10(P.Value),label=gene),nudge_x= 1,nudge_y= 5,
                     color= "white",alpha= 0.9,point.padding = 0.5,size= 6,fill= "#6A5ACD",
                     segment.size = 0.5,segment.color = "grey50",direction= "y",hjust= 0.5)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 24),
          plot.title = element_text(size = 24))+
    guides(color = guide_legend(override.aes = list(size = 5),title = NULL))
}
DEGs_file$sig <- factor(DEGs_file$sig,levels = c('up','down','none'),ordered = T)
#??|log2FC|??????չʾǰ8??????
top_up <- DEGs_file %>% filter(sig == "up") %>% top_n(4, abs(logFC))
top_down <- DEGs_file %>% filter(sig == "down") %>% top_n(4, abs(logFC))
top <- rbind(top_up,top_down)
top$gene <- rownames(top)
# pdf("../data/01_DEGs/volcano_miRNA.pdf",width=10,height=8,onefile=FALSE)
p1 <- volcano_plot(DEGs_file)
ggsave(filename = 'volcano.pdf',plot = p1,width = 8,height = 6)
#ggsave(filename = 'volcano.png',plot = p1,width = 8,height = 6)


#####??ͼ
library(ComplexHeatmap)
library(circlize)

deg <- DEGs_file %>% filter(sig != "none")
DEG_id1 <- rownames(deg)
DEG_exp <- expr[DEG_id1,]
hmexp <- na.omit(DEG_exp)

eset_dat <- hmexp
index <- DEG_id1
index_matrix <- t(scale(t(eset_dat[index,])))##??һ??
index_matrix[index_matrix >= 1] = 1
index_matrix[index_matrix <= -1] = -1
hmexp <- index_matrix
col_fun <- colorRamp2(
  c(-1, 0, 1), 
  c("blue", "#EEEEEE", "red"),space = "RGB"
)
pdf("pheatmap.pdf", width=7, height=6)

ha1 = HeatmapAnnotation(group = as.factor(c(rep("Control", 40),rep("OP",40))),col =list(group = c("OP" = "pink", "Control" = "#95d6f0")))
densityHeatmap(hmexp, ylim = c(-2, 2),quantile_gp = gpar(fontsize = 7),title = " ", 
               ylab = "Expression",top_annotation = ha1)%v%
  #HeatmapAnnotation(Expression = anno_barplot(1:length(colnames(hmexp)))) %v%
  Heatmap(hmexp, name = "expression",col =col_fun,show_row_names = FALSE,show_column_names = FALSE, cluster_rows =TRUE,height = unit(9, "cm"))
dev.off()

write.csv(DEGs_file,file = "all_gene.csv")
write.csv(deg,file = "DEG.csv")


############Mito_DEG??PCD_DEGȡ????
setwd("E:/FZZK-20104-10/02_venn")

library(ggVennDiagram)
library(hrbrthemes)
library(ggplot2)
library(tidyverse)
library(openxlsx)

rm(list = ls())
mito <- read.table("E:/FZZK-20104-10/00_rawdata/Mito-RGs/Mito-RGs.txt",header = T)
pcd <- read.xlsx("E:/FZZK-20104-10/00_rawdata/PCD/PCD.xlsx")
deg <- read.csv("E:/FZZK-20104-10/01_DEG/DEG.csv")

names(mito)[1] <- "ID"
names(pcd)[1] <- "ID"
names(deg)[1] <- "ID"

set.seed(20210419)
x <- list(A=sample(mito$ID),
          B=sample(pcd$ID),
          C=sample(deg$ID))


pdf("inter_DEGs_venn.pdf",height = 6,width = 6)
ggVennDiagram(x, category.names = c("Mito-RGs","PCD-RGs","DEGs"),
              size=1, label_alpha = 0,lty="longdash",color="gray60",label = "count") + 
  scale_fill_gradient(name="Count",low="white",high = "red") +
  hrbrthemes::theme_ipsum(base_family = "sans") +
  theme_bw() + theme(panel.grid=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
                                  size = 20, margin = margin(t = 1, b = 12)),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()


####??????ȡ
venn <- Venn(x)
df <- process_data(venn)
inter_method <-df@region[["item"]][[7]]
write.table(inter_method,file = "inter_DEGs.txt",sep = "\t")

################GO/KEGG 
setwd("E:/FZZK-20104-10/03_enrichment")


rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(stringr)
library(tidyr)
library(GOplot)


DiffGene <- read.table("E:/FZZK-20104-10/02_venn/inter_DEGs.txt",sep = "\t")
DEG <- read.csv("E:/FZZK-20104-10/01_DEG/DEG.csv",header = T,row.names = 1,check.names = F)
DEG_18 <- DEG[DiffGene$x,]
write.csv(DEG_18,file = "01_DEG_18.csv")


Gene = bitr(rownames(DEG_18),fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")  #??SYMBOLת??ΪENTREZIDGene <- na.omit(Gene) 
#FClist <- t(DiffGene)  #ֻ??gene??FC????Ϣ
#FClist <- FClist[1,]
write.csv(as.data.frame(Gene),"02_gene_symboltoentrezid.csv",row.names =F)  #????????


##KEGG????????############
kegg <- enrichKEGG(gene=Gene$ENTREZID, 
                   organism = 'hsa',
                   keyType = 'kegg',
                   pvalueCutoff =0.05, 
                   qvalueCutoff=1,
                   pAdjustMethod = 'BH'
)
#gene ID ransform
##clusterprofile?Դ???ת??
kegg <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
head(kegg)
dim(kegg)

data_kegg <- as.data.frame(kegg)
sum(data_kegg$pvalue<0.05)
write.csv(as.data.frame(kegg),"03_KEGG_enrich.csv",row.names =F)

###########????ͼ################################################################

KEGG <- data_kegg[,c(1,3,4,8,10)]
KEGG$geneID <- str_replace_all(KEGG$geneID,"/",",")##??/??Ϊ,


colnames(KEGG) <- c("category","ID","term","adj_pval","genes")

DEG_18$ID <- rownames(DEG_18)
gene <- DEG_18[,c(8,1)]
names(gene)[2] <- "logFC"
rownames(gene) <- NULL
circ <- circle_dat(KEGG,gene)
pdf(file="04_KEGG_circ.pdf",width=10,height=8,onefile=FALSE)
GOCircle(circ,nsub = 8)
dev.off()


##1.3GO????????#################
ego<-enrichGO(gene=Gene $ENTREZID, 
              OrgDb=org.Hs.eg.db, 
              keyType='ENTREZID', 
              ont="ALL", 
              pAdjustMethod="BH", 
              pvalueCutoff=0.05, 
              qvalueCutoff=1, 
              readable=TRUE) 
head(ego)
write.csv(as.data.frame(ego),"06_GOenrich.csv",row.names =F)
sum(ego$ONTOLOGY=="BP") #Biological process??????????????????·????????
sum(ego$ONTOLOGY=="CC") #Cellular component??????????ϸ????????λ??
sum(ego$ONTOLOGY=="MF") #Molecular function???????????Ӳ??εĹ???
data_GO <- as.data.frame(ego)

#######????ͼ###################################################################
go <- data_GO[,c(1,2,3,8,9)]
go$geneID <- str_replace_all(go$geneID,"/",",")

colnames(go) <- c("category","ID","term","adj_pval","genes")


circ <- circle_dat(go,gene)
pdf(file="07_GO_circ.pdf",width=10,height=8,onefile=FALSE)
GOCircle(circ,nsub = 8)
dev.off()



#################################????ѧϰ##########################################
setwd("E:/FZZK-20104-10/05_machine learn/lasso")

rm(list = ls())
library(glmnet)
library(tidyr)
library(pROC)
library(data.table)

expr <- fread("E:/FZZK-20104-10/00_rawdata/GSE56815/03_expr.csv")
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- expr[,-1]
expr <- t(expr)


Gene <- read.table("E:/FZZK-20104-10/02_venn/inter_DEGs.txt")
Gene <- Gene$x
Gene <- intersect(Gene,colnames(expr))
Gene_expr <- expr[,Gene]
write.csv(Gene_expr,file ="GSE56815_genexpr.csv" )


group <- c(rep("0",40), rep("1", 40)) %>%
  factor(., levels = c("0", "1"), ordered = F)

design <- model.matrix(~factor(group))  #design֮ǰ??DEGʱ?ͳ????ˣ?????????Ϊû?б??????????¼???һ??
gene_expr <- cbind(design[, 2], Gene_expr)
colnames(gene_expr)[1] <- "group"
rownames(gene_expr) <- rownames(Gene_expr)
gene_expr <- as.matrix(gene_expr)
data <- as.data.frame(gene_expr)
write.csv(data,file = "01_input.csv")


#data <- read.csv("E:/CDDZ-0404-8/05_machine learn/lasso/01_input.csv",header = T,row.names = 1,check.names = F)
#data[,1] <- as.factor(data[,1])#??????group


###????lasso?ع?
lasso.v <- as.matrix(data)
set.seed(44)
fit_stage.v <- cv.glmnet(lasso.v[, -1], lasso.v[, "group"], family = "binomial",
                         type.measure = "default", nfolds = 3, gamma = 1)
fit.v <- fit_stage.v$glmnet.fit
pdf(file = 'lasso.Binomial.Deviance.pdf',height = 8,width = 8)
par(mgp = c(2.5, 1, 0), mar = c(4, 5, 3, 3))
plot(fit_stage.v, xlab = "Log Lambda", cex.lab = 2.4) + text(x = log(fit_stage.v$lambda.min),
                                                             y = 0.8, paste("Lambda.min\n", round(fit_stage.v$lambda.min, 3)), cex = 2, adj = 0.5) +
  text(x = log(fit_stage.v$lambda.1se), y = 0.9, paste("Lambda.lse\n", round(fit_stage.v$lambda.1se,
                                                                             3)), cex = 2)
dev.off()


pdf(file = 'lasso.voefficients.venalty.pdf',height = 8,width = 8)
par(mgp = c(4, 1, 0), mai = c(2, 2, 1, 1))
plot(fit.v, xvar = "lambda", cex.lab = 3) + abline(v = c(log(fit_stage.v$lambda.min),
                                                         log(fit_stage.v$lambda.1se)), lty = 2) + text(x = log(fit_stage.v$lambda.min),
                                                                                                       y = 0.5, paste("Lambda.min\n", round(fit_stage.v$lambda.min, 4)), cex = 1.6,
                                                                                                       adj = 0.9) + text(x = log(fit_stage.v$lambda.1se), y = 1, paste("Lambda.lse\n",
                                                                                                                                                                       round(fit_stage.v$lambda.1se, 4)), cex = 1.6, adj = 0.9)
dev.off()

###??һ??lassoԤ??Ч??

library(dplyr, warn.conflicts = F)
library(kableExtra, warn.conflicts = F)
lam = fit_stage.v$lambda.min
pre_train <- predict(fit_stage.v, type = "class", newx = lasso.v[, -1], s = lam)
acc_table_tr <- table(pre_train, lasso.v[, "group"])
colnames(acc_table_tr) <- c("Control", "OP")
rownames(acc_table_tr) <- c("Control", "OP")
suppressMessages(library(kableExtra))
kable(acc_table_tr, "html", caption = "<center>**??. LASSO(Training set) **</center>") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

###??ȡ????????
coefficients <- coef(fit_stage.v, s = lam)
Active.Index <- coefficients@i
Active.voefficients <- coefficients@x
Active.name.train <- colnames(lasso.v[, -1])[Active.Index[-1]]
write.table(Active.name.train,file = "fea_gene.txt",sep = "\t")

###ROC????
probability_train<-predict(fit_stage.v,type = "response", newx=lasso.v[,-1], s=lam)
roc_train<-roc(as.vector(lasso.v[,'group']),as.vector(probability_train),percent=T)

pdf("roc_pre.pdf",height = 10,width = 8)
gbm.ROC <- plot.roc(lasso.v[,'group'], probability_train,ylim=c(0,1),xlim=c(1,0),
                    smooth=F, #????ƽ??????
                    ci=TRUE,
                    main="",
                    col="red",#?ߵ???ɫ
                    lwd=2, #?ߵĴ?ϸ
                    legacy.axes=T,
                    print.auc = T)
dev.off()



##################################svm-rfe##########################################################################

setwd("E:/FZZK-20104-10/05_machine learn/svm-rfe")
rm(list = ls())

library(tidyverse)
library("e1071")
library(caret)
library(ggplot2)
library(lattice)
#install.packages("mlbench")
library(mlbench)
#BiocManager::install("sigFeature")
library(sigFeature) 
#install.packages("bRacatus")
library(bRacatus)
source("01_msvmRFE.r")



data <- read.csv("E:/FZZK-20104-10/05_machine learn/lasso/01_input.csv",header = T,row.names = 1,check.names = F)
data[,1] <- as.factor(data[,1])#??????group


set.seed(123)

svmRFE(data, k = 10, halve.above = 100)

nfold =10
nrows = nrow(data)
folds = rep(1:nfold, len = nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, data, k = 10, halve.above = 100)


# save(results,file = 'result.c.rda')
top.features = WriteFeatures(results, data, save = F)

write.csv(top.features,"02_feature_svm.csv") 

#install.packages("doparallel")
library("doParallel")

#install.packages("parallel")
library(parallel)
cl<- makeCluster(5)
registerDoParallel(cl)
#########################????ʱ???̾ã?????ʮ????????
mydata1c <- foreach(
  i=1:18,          #?????ȴ??????Ĳ???
  .combine=rbind,  #???ؽ?????????
  .packages = c("e1071")
  #???????̹?????ϵͳ????
) %dopar%  FeatSweep.wrap(i, results, data)
stopCluster(cl)
no.info = min(prop.table(table(data[,1])))
errors = as.data.frame(mydata1c)
errors = errors$error

#???ƻ???SVM-REF?㷨?Ĵ?????????ͼ 
pdf(file = "03_svm_error_10??.pdf",height = 5,width = 5)
PlotErrors(unlist(errors), no.info=0.5)
dev.off()

#???ƻ???SVM-REF?㷨????ȷ??????ͼ

#dev.new(width=4, height=4, bg='white')

#pdf("04_svm-accuracy.pdf",width = 5,height = 5)

#accuracy <- as.data.frame(1-unlist(errors))
#colnames(accuracy)[1] <- "accuracy"
#plotAccuracy(accuracy,regional = F) #?鿴׼ȷ??

#dev.off() 

fea_gene <- top.features[1:10, 1]
set.seed(43)
cancerdata <- data
cancerdata$group <- as.factor(data$group)
cancerdata <- cancerdata[, c(1, which(colnames(cancerdata) %in% fea_gene))]

library(kableExtra, warn.conflicts = F)
can_svm <- svm(group ~ ., data = cancerdata, kernel = "linear")
summary(can_svm)
train_pre <- as.character(predict(can_svm, cancerdata))
sprintf("base SVM acc: %f", recall(cancerdata$group, as.factor(train_pre)))
acc <- table(cancerdata$group, train_pre)
colnames(acc) <- c("Control", "RPL")
rownames(acc) <-c("Control", "RPL")
kable(acc, "html", caption = "<center>**??. SVM(Training set) **</center>") %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

write.table(fea_gene, file = "04_fea_gene_10??.txt", sep = "\t", row.names = T, quote = F)



#########################��???㷨ȡ????############################################
setwd("E:/FZZK-20104-10/05_machine learn/venn")
rm(list = ls())


method1 <- read.table("E:/FZZK-20104-10/05_machine learn/lasso/fea_gene.txt",sep = "\t")
method2 <- read.table("E:/FZZK-20104-10/05_machine learn/svm-rfe/04_fea_gene_10??.txt",sep = "\t")


library(ggVennDiagram)
library(hrbrthemes)
library(ggplot2)
library(tidyverse)

names(method1)[1] <- "ID"
names(method2)[1] <- "ID"


set.seed(20210419)
x <- list(A=sample(method1$ID),
          B=sample(method2$ID))


pdf("method_venn.pdf",height = 6,width = 6)
ggVennDiagram(x, category.names = c("LASSO","SVM-RFE"),
              size=1, label_alpha = 0,lty="longdash",color="gray60",label = "count") + 
  scale_fill_gradient(name="Count",low="white",high = "#FC8C5A") +
  hrbrthemes::theme_ipsum(base_family = "sans") +
  theme_bw() + theme(panel.grid=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
                                  size = 20, margin = margin(t = 1, b = 12)),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
dev.off()

####??????ȡ
venn <- Venn(x)
df <- process_data(venn)
inter_method <-df@region[["item"]][[3]]
write.table(inter_method,file = "inter_method.txt",sep = "\t")





################ѵ��???б?????֤

setwd("E:/FZZK-20104-10/06_expression")
rm(list = ls())
library(tidyr)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(limma)
library(data.table)

expr <- fread("E:/FZZK-20104-10/00_rawdata/GSE56815/03_expr.csv")
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- expr[,-1]
expr <- t(expr)


Gene <- read.table("E:/FZZK-20104-10/05_machine learn/venn/inter_method.txt")
Gene <- Gene$x
Gene <- intersect(Gene,colnames(expr))
Gene_expr <- expr[,Gene]
write.csv(Gene_expr,file ="GSE56815_genexpr.csv" )


data <- data.frame(Gene_expr,"group"=rep(c("Control","OP"),c(40,40)))
mydata<-data %>% 
  ## ????????????gather,gather?ķ?ΧӦ????
  gather(key="gene",value="Expression",na.rm = T,-group)

str(mydata)
#mydata$Expression <- as.numeric(mydata$Expression)

P<-mydata%>% 
  ## ȷ??x,y
  ggplot(aes(x = gene, y = Expression, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=c("blue","red"))+
  scale_y_continuous(name = "Expression value")+
  scale_x_discrete(name = "gene") +
  ggtitle("GSE56815") +
  theme_bw() +
  theme(plot.margin=unit(rep(3,4),'lines'),
        plot.title = element_text(size = 14,hjust=0.5),
        text = element_text(size = 16),
        axis.text.x=element_text(face = "bold",size = 12,angle=90),
        axis.title.y=element_text(face = "bold",size = 14),#y????ǩ?Ӵּ???????С
        axis.title.x=element_text(face = "bold",size = 12),#X????ǩ?Ӵּ???????С
        axis.text.y = element_text(face = "bold",size = 12))#y???????̶ȱ?ǩ?Ӵ?
pdf(file="GSE56815_boxplot.pdf",width=8,height=8,onefile=F)
p<-P+stat_compare_means(label = "p.signif", method = "wilcox.test")
print(p)
dev.off()



##############??֤???б?????֤

rm(list = ls())
expr <- fread("E:/FZZK-20104-10/?½??ļ???/00_rawdata/GSE13850/03_expr.csv")
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- expr[,-1]
#expr <- t(expr)

Gene <- read.table("E:/FZZK-20104-10/05_machine learn/venn/inter_method.txt")
Gene <- Gene$x
Gene <- intersect(Gene,colnames(expr))
Gene_expr <- expr[,Gene]
write.csv(Gene_expr,file ="GSE13850_genexpr.csv" )


data <- data.frame(Gene_expr,"group"=rep(c("Control","OP"),c(10,10)))
mydata<-data %>% 
  ## ????????????gather,gather?ķ?ΧӦ????
  gather(key="gene",value="Expression",na.rm = T,-group)

str(mydata)
#mydata$Expression <- as.numeric(mydata$Expression)

P<-mydata%>% 
  ## ȷ??x,y
  ggplot(aes(x = gene, y = Expression, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values=c("blue","red"))+
  scale_y_continuous(name = "Expression value")+
  scale_x_discrete(name = "gene") +
  ggtitle("GSE13580") +
  theme_bw() +
  theme(plot.margin=unit(rep(3,4),'lines'),
        plot.title = element_text(size = 14,hjust=0.5),
        text = element_text(size = 16),
        axis.text.x=element_text(face = "bold",size = 12,angle=90),
        axis.title.y=element_text(face = "bold",size = 14),#y????ǩ?Ӵּ???????С
        axis.title.x=element_text(face = "bold",size = 12),#X????ǩ?Ӵּ???????С
        axis.text.y = element_text(face = "bold",size = 12))#y???????̶ȱ?ǩ?Ӵ?
pdf(file="GSE13580_boxplot.pdf",width=12,height=8,onefile=F)
p<-P+stat_compare_means(label = "p.signif", method = "wilcox.test")
print(p)
dev.off()


##################nomograph
setwd("E:/FZZK-20104-10/08_nomo")

library(ggDCA)
library(rms)
library(nricens)
library(foreign)
library(rmda)
library(regplot)

rm(list = ls())
data3 <- read.csv( "01_input.csv",header = T,check.names = F,row.names = 1)#????????
#1.????????ͼ????
dd <- datadist(data3)
options(datadist="dd")
f <- lrm(group~ MCL1+BIK+PMAIP1+TRAP1, x=T, y=T, data=data3)
print(f)
pdf("nomogram.pdf",height = 8,width = 8)
nom <- nomogram(f, fun=plogis,
                lp=F, 
                fun.at = c(0.1,0.4,0.8,0.99),
                funlabel = "Risk")#ִ??Nomogram????
par(mar = c(1, 3, 2, 2))
plot(nom,xfrac=.5,cex.axis=0.8,cex.var=1)
dev.off()
summary(f)

#?߶?ʽ??̬ŵĪͼ
#??ŵĪͼҲ???Գ?Ϊ????ʽ????ͼ?????? regplot ?????н? clickable ????Ϊ TRUE(clickable = T ) ???????ǾͿ?????ͼ?ϵ???ѡ??ÿ????��?Ĳ?ͬ????????ͬ??ֵ??��????ͬ???????߶?Ӧ???ܰͽ?ת?Ʒ??գ?ʵ?ֽ??????ܡ????⣬??????ѡ????ͬ????ʾ???񣨼?Figure 4???????뼰???????£?

#pdf("nomo_trend.pdf",height = 10,width = 10)
#regplot(f,
#observation = data3[1,],#ѡ?????ݼ???ĳһ?о?????
#interval = "confidence",#????????????
#title = "Nomogram",#ͼƬ????
#clickable = T)#??ʾ?????Ƿ??ֶ?ѡ????TΪ??
#dev.off()



#2.??��У׼???߲?????????ͼ
pdf("calibrate.pdf",height = 8,width = 8)
cal1<-calibrate(f,method="boot",B=1000)
par(mar=c(6,5,1,2),cex = 1.0)#ͼ?β???
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0),
     xlab = "Predicted Pr (PCOS=1)", ylab = "Actual Probablity")
dev.off()


#3.???????߻???
library(rmda)
set.seed(123)

MCL1<- decision_curve(group~ MCL1,data = data3, 
                     family = binomial(link ='logit'),#ģ?????ͣ??????Ƕ?????
                     thresholds= seq(0,1, by = 0.01),
                     confidence.intervals =0.95,#95????????
                     study.design = 'cohort')#?о????ͣ??????Ƕ????о?
TRAP1<- decision_curve(group~ TRAP1,data = data3, family = binomial(link ='logit'),
                        thresholds= seq(0,1, by = 0.01),
                        confidence.intervals =0.95,study.design ='cohort')
BIK<- decision_curve(group~ BIK,data = data3, family = binomial(link ='logit'),
                     thresholds= seq(0,1, by = 0.01),
                     confidence.intervals =0.95,study.design ='cohort')
PMAIP1<- decision_curve(group~ PMAIP1,data = data3, family = binomial(link ='logit'),
                     thresholds= seq(0,1, by = 0.01),
                     confidence.intervals =0.95,study.design ='cohort')
nomogram<- decision_curve(group  ~ MCL1+TRAP1+BIK+PMAIP1,data = data3,
                          family = binomial(link='logit'),
                          thresholds= seq(0,1, by = 0.01),
                          confidence.intervals =0.95,study.design ='cohort')

pdf("DCA.pdf",height = 8,width = 8)
List<-list(MCL1,TRAP1,BIK,PMAIP1,nomogram)
plot_decision_curve(List,curve.names= c('MCL1','TRAP1','BIK','PMAIP1','nomogram'),
                    cost.benefit.axis =T,col = c('red','blue','green','yellow','brown'),
                    confidence.intervals =FALSE,standardize = F,
                    legend.position = "topright")#legend.position = "none"
dev.off()

#4.?????ٴ?Ӱ?????ߣ?Clinical Impact Curve????
pdf("clinical_impact.pdf",height = 8,width = 8)
plot_clinical_impact(nomogram,population.size= 1000,
                     cost.benefit.axis = T,
                     n.cost.benefits= 8,
                     col =c('red','green'),
                     confidence.intervals= T,
                     ylim=c(0,1000),
                     legend.position="none")#legend.position="bottomleft"
dev.off()


nomo <- nomogram(f)
#install.packages("nomogramFormula")
library(nomogramFormula)##????nomogramFormul(nomogram=nomo)
data3$points<-points_cal(formula = results$formula,rd=data3)##????ÿ??????????
ppdf("roc.pdf",height = 10,width = 8)
plot.roc(data3$group, pre,ylim=c(0,1),xlim=c(1,0),
                    smooth=F, #????ƽ??????
                    ci=TRUE,
                    main="",
                    col="red",#?ߵ???ɫ
           lwd=2, #?ߵĴ?ϸ
           legacy.axes=T,
                    print.auc = T)
dev.off()

###################################GSEA(????KEGG

########BIK
stewd(d("E:/FZZK-20104-10/09_GSEA/go")
rm(list = ls())

library(data.table)
expr <- fread("E:/XA-0403-7/00_rawdata/GSE92566/03_expr.csv",header = T)
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- t(expr[,-1])

batch_cor <- function(gene){
  y <- as.numeric(expr[gene,])
  rownames <- rownames(expr)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(expr[x,]),y,type="spearman")
    data.frame(gene=gene,Symbol=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

library(future.apply)
#plan(multiprocess)
#nbrOfWorkers()
system.time(dd <- batch_cor("BIK")) 
gene <- dd$Symbol


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(RColorBrewer)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## ȥ??
gene <- dply## ȥ??
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(cor=dd$cor,
                      SYMBOL = dd$Symbol)
gene_df <- merge(gene_df,gene,by="SYMBOL")

## 
##
geneList <- gene_df$cor
## 
names(geneList) = gene_df$ENTREZID
## 
geneList = sort(geneList, decreasing = TRUE)

#2 
##2.1
#

#devtools::install_github("junjunlab/GseaVis")
library(GseaVis)

kegg_gmt <- read.gmt("01_c2.cp.kegg.v7.5.1.entrez.gmt")
gsea <- GSEA(geneList,
             TERM2GENE = kegg_gmt) (gsea, 'org.Hs.eg.db', 'ENTREZID')
head(gsea)[1:10]
write.csv(as.data.frame(gsea),"08_GSEA_BIK.csv",row.names = F)

geneSetID = c('KEGG_RIBOSOME',
              'KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION',
              'KEGG_OLFACTORY_TRANSDUCTION',
              'KEGG_CALCIUM_SIGNALING_PATHWAY',
              'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')
# all plot
pdf(file ="09_GSEA_kegg_BIK.pdf",width=12,height=8,onefile=FALSE )
gseaNb(object = gsea,
       geneSetID = geneSetID,termWidth = 35,
       curveCol = brewer.pal(5,'Paired'))
dev.off()



###################################GSEA(????GO

########BIK
stewd(d("E:/FZZK-20104-10/09_GSEA/go")
rm(list = ls())

library(data.table)
expr <- fread("E:/FZZK-20104-10/00_rawdata/GSE56815/03_expr.csv")
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- expr[,-1]

batch_cor <- function(gene){
  y <- as.numeric(expr[gene,])
  rownames <- rownames(expr)
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(expr[x,]),y,type="spearman")
    data.frame(gene=gene,Symbol=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

library(future.apply)
#plan(multiprocess)
#nbrOfWorkers()
system.time(dd <- batch_cor("BIK")) 
gene <- dd$Symbol
## ת??
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(RColorBrewer)
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## ȥ??
gene <- dplyr::inct(gene,SYMBOL,.keep_all=TRUE)

gene_df <- data.frame(cor=dd$cor,
                      SYMBOL = dd$Symbol)
gene_df <- merge(gene_df,gene,by="SYMBOL")

## geneList ??????
## 1gene$cor
## 2.????
names(geneLisene_df$SYMBOL
## 3.????????Ҫ
geneLiseList, decreasing = TRUE)

#2 GSEA????#############//_github("junjunlab/GseaVis")
library(GseaVis)

go_gmt <- read.gmt("01_c5.go.bp.v2023.2.Hs.symbols.gmt") #??gmt?ļ?
gsea <- Gst,
             TERM2GENE = go_gmt) #GSEA????
#gsea <- se(gsea, 'org.Hs.eg.db', 'ENTREZID')
#head(gsea)[1:10]
write.csv(as.data.frame(gsea),"08_GSEA_BIK.csv",row.names = F)

geneSetID = c('GOBP_CYTOPLASMIC_TRANSLATION',
              'GOBP_TRANSLATIONAL_INITIATION',
              'GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION',
              'GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS',
              'GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS')
# all plot
pdf(file ="09_GSEA_go_BIK.pdf",width=10,height=8,onefile=FALSE )
gseaNb(object = gsea,
       geneSetID = geneSetID,termWidth = 35,
       curveCol = brewer.pal(5,'Paired'))
dev.off()

#########################immune
setwd("E:/FZZK-20104-10/10_immune")
rm(list = ls())

library(data.table)
library(dplyr)
library(RColorBrewer)
library(tidyr)
#install.packages("radiant.data")
library(radiant.data)
library(aplot)
library(cowplot)
library(pheatmap)
library(ggpubr)
library(ggplot2)
library(GSEABase)
library(GSVA) #????ssGSEA scores

xpr <- fread("E:/FZZK-20104-10/00_rawdata/GSE56815/03_expr.csv")
expr <- as.data.frame(expr)
rownames(expr) <- expr$V1
expr <- expr[,-1]
#expr <- t(expr)
group <- read.csv("E:/FZZK-20104-10/00_rawdata/GSE56815/04_group.csv",header = T,row.names = 1,check.names = F)

#????????ϸ??
gmtFile="01_mmc3.gmt"

#2 ????????ϸ????��#########
#???뱳??????
geneSet=getGmt(gmtFile,
               geneIdType=SymbolIdentifier())

##2.1??ʼssGSEA????,????һ??Ҫ?Ǿ??󣬲??????ݿ?###########
ssgseaScore=gsva(as.matrix(expr), geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

range(ssgseaScore)
ssgseaOut <- as.data.frame(ssgseaScore)
write.csv(ssgseaOut,file = "02_immunecell_ssgseaout.csv")


library(pheatmap)
pdf(file="03_immunecell_pheatmap.pdf",width=10,height=8,onefile=FALSE)
annotation_col=data.frame(Group=rep(c("Control(n=40)","OP(n=40)"),c(40,40)))
annColors <- list(Group = c("Control(n=40)" = "#1CB4B8","OP(n=40)" ="#EB7369"))
rownames(annotation_col)=colnames(ssgseaOut)
pheatmap(ssgseaOut,
         annotation_col = annotation_col,
         annotation_colors = annColors,
         color = colorRampPalette(c("#DD2634", "white", "#206FB0"))(100),
         cluster_rows=TRUE,
         show_rownames=TRUE,
         fontsize_row=10,
         show_colnames=FALSE,
         scale="row",
         cluster_cols=FALSE,
         main="")
dev.off()


##2.5??????ϸ???Ĳ???###################
library(ggpubr)
library(ggplot2)
library(grid)


#ͼƬ?��?
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#ȥ??????
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#ȥx????ǩ
          axis.title.y=element_text(face = "bold",size = 14),#y????ǩ?Ӵּ???????С
          axis.title.x=element_text(face = "bold",size = 14),#X????ǩ?Ӵּ???????С
          axis.text.y = element_text(face = "bold",size = 12),#y???????̶ȱ?ǩ?Ӵ?
          axis.text.x = element_text(face = "bold",size = 10, vjust = 1, hjust = 1, angle = 45),#x???????̶ȱ?ǩ?Ӵ?
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position=c(0.9, 0.8),#ͼ???ڻ?ͼ??????λ??
          # legend.position="top",
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}


#rownames(group) <- group[,1]
data <- merge(group,t(ssgseaOut),by="row.names")
rownames(data) <- data[,1]
data <- data[,-1]

re1 <- data[order(data[,1],decreasing = F),]

table(re1$group)  
#ȥ????????Ϣ
re2 = re1[,-1]       #???з?????Ϣ?????ݡ?
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
dat_cell <- re2 %>% as.data.frame() %>%rownames_to_column("Sample") %>%gather(key = Cell_type,value = Proportion,-Sample)
#??ȡ????
dat_group = gather(re1,Cell_type,Proportion,-group )
#?ϲ?????
dat = cbind(dat_cell,dat_group$group)
colnames(dat)[4] <- "group"

e1 <- ggplot(dat,aes(x=Cell_type,y=Proportion),palette = "jco", add = "jitter")+
  geom_boxplot(aes(fill=group),position=position_dodge(0.6),width=0.6)+
  labs(x = "Cell Type", y = "Immune cell content")+
  scale_fill_manual(values = c("#1CB4B8","#EB7369")) +
  theme_zg()+theme(legend.position = "top") 
e1 = e1 + stat_compare_means(aes(group = group),label = "p.signif",method = "wilcox.test")
e1
#????ͼƬ
ggsave('04_ssGSEAbox_ImmuneGeme_wt.pdf', plot = e1,width=15,height = 6)


###########hub??????????ϸ????????#########################################
gene <- "BIK"
diff_immunecell <- c("CD56dim natural killer cell","Central memory CD4 T cell")
immunecell_expr <- ssgseaOut[diff_immunecell,]


library(psych)


d <- corr.test(t(gene_expr),t(immunecell_expr),use="complete",method = 'spearman')

r <- d$r
p <- d$p

write.csv(r,file = "05_cell_gene_cor.csv")
write.csv(p,file = "06_cell_gene_p.csv")


library(corrplot)
pdf(file="07_cor_pheatmap.pdf",width=10,height=8,onefile=FALSE)
corrplot(r, p.mat = p, method = "circle", insig = "label_sig",tl.col="black",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")

dev.off()


























