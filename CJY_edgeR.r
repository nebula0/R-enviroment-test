# Input: GelataTuber.featureCounts.txt
# Output: DEG list

# 0. Environment ----------------------------------------------------------------------------------------

library(edgeR)
library(WGCNA)
library(RUVSeq)
library(RColorBrewer)
library(tidyverse)

rm(list = setdiff(ls(), lsf.str()))
setwd("/Users/miachen/Desktop")

# 1. Loading data -------------------------------------------------------------------------------------------

input<- read.table("GelataTuber.featureCounts.txt",header=T, row.name=1,check.names=F)

#input<- read.table("../IMB_N317.featureCounts.txt",header=T, row.name=1)
# reorder for easy processing

#data<-input[,grep("_R_",colnames(input))]

#1.2  Labeling ---------------------------------------------------------------------------------------------

Labels=unlist(strsplit(colnames(input),"[_]"))
Labels

Sample_number<-factor(Labels[keep=c(T,F,F)],levels=c("S4","S5","S6","S7","S8","S9"))
Sample_number
Part<-factor(Labels[keep=c(F,F,T)],levels=c("up","basal"))
Part
rep<-factor(rep(c("1","2"),each=3))
rep

## genotype_condition for ploting
group<-factor(paste(Sample_number,Part,sep='_'))
group

sample.info=data.frame(group,Sample_number,Part)
sample.info

raw <- DGEList(counts=input,group=rep,remove.zeros = T)
nrow(raw) #16604

collapse<-collapseRows(t(cpm(raw)),rowGroup=rep,rowID=colnames(raw),method = "Average")
#collapse<-collapseRows(t(cpm(raw)),rowGroup=rep,rowID=colnames(raw),method = "Average")
keep<-rowSums(t(collapse$datETcollapsed)>=2)>=1
filtered <- raw[keep, keep.lib.sizes=F] 
nrow(filtered) #13755

table(is.na(filtered$counts))
genes <- rownames(filtered)


# 2. Normalization -------------------------------------------------------------------------------

genes <- rownames(filtered)

#write.table(genes, quote=FALSE, row.names=F, col.names=F,sep = "\t",
#            file="./List/Tsay-Rice-TNGTCN-root.universe-16780geneID.txt")

##Before normalization
##轉換檔案形式
set <- newSeqExpressionSet(as.matrix(filtered), 
                           phenoData = data.frame(Sample_number,Part,
                                                  row.names =colnames(filtered)))

## After normalization
##挑選上部分（活躍的數值）;上－部分異動概率較低
set1 <- betweenLaneNormalization(set, which="upper")

set1

#RUVs: Estimating the factors of unwanted variation using replicate samples
##分組，告知RUVs哪些為一組（數值會較相似）
differences<-makeGroups(Part)

differences

set3 <- RUVs(set1, genes, k=1, differences)

pData(set3)

par(bg = "white",cex.axis=0.8,cex.lab = 0.8,mar=c(5.1, 3.9, 3.9,11), xpd=TRUE)
plotRLE(set3, outlone=F,ylim=c(-3,3))

plotPCA(set3,cex=1.5,labels=F,pch=c(2,17)[Part])

# 3. DE analysis ---------------------------------------------------------------------------

## 3.1 DE: 

sample.info
#Control

design1 <- model.matrix(~Part+W_1,data=pData(set3)) 
design1
colnames(design1)

##檔案轉換
y <- DGEList(counts=counts(set3),group=Part)
y <- calcNormFactors(y, method = 'TMM')
y$samples

y<-estimateDisp(y,design1,robust = T)
y$common.dispersion #0.01538511

##檔案轉換，可查看edgeR手冊
fit <-glmQLFit(y, design1)

##挑選感興趣之性狀
colnames(design1)
qlf<-glmQLFTest(fit,coef=2)

##將檔案轉換成table
top <- topTags(qlf, n=nrow(set3))$table

###filter
sum(top$FDR<0.05) #7274
sum(top$FDR<0.01) #5148

de <- topTags(qlf, n=sum(top$FDR<0.05))
de.gene<-rownames(de)
de.geneinfo<-de$table

#write.table(de.gene, quote=F, row.names = F, col.names = F,sep="\t",
            file=("GelataTuber_7179.txt")

