# Goal: This script analyzes Gastrodia tuber PE RNA-seq data generated by Lee lab
# Input: feacutreCounts read counts; including 3 top and 3 bottom (infected by armillaria gallica) 
# Output: DEG list
# Conclusion: 
# Ref: 
# Author: Chia-Yi Cheng
# Last updated: 2022-07-02

# 0. Environment ---------------------------------------------------------------

setwd("C:/R_CJY/WGCNA_test")

library(edgeR)
library(WGCNA)
library(RUVSeq)
library(RColorBrewer)


rm(list = setdiff(ls(), lsf.str()))


# 1. Loading data --------------------------------------------------------------
input <- read.table("GelataTuber.featureCounts.txt",header=T, row.name=1,check.names=F)
colnames(input)

#1.1  Labeling -----------------------------------------------------------------

condition <- factor(rep(c("up","basal"),each=3),levels=c("up","basal"))
condition

rep <- rep(c("1","2"),each=3) #? rep meaning?
rep

sample.info = data.frame(condition,rep)
sample.info

#1.2 Filter out lowly expressed genes -----

## convert input into a DGEList
raw <- DGEList(counts=input,group=rep,remove.zeros = T)
nrow(raw) #16604
head(raw)
head(cpm(raw))
## Filgering criteria: the average of 3 rep is above 2 CPM

#- use `cpm` function divide each row with a number, which i don't know how cpm 
#- function calculate that number. after that, exchange row and column by `t`, 
#- now row are samples and column are genes.
#- `collapseRows` function select one representative row per group, in this case, 
#- our representative row is the average gene expression of each group(up 1 or 
#- bottom 2). 
collapse <- collapseRows(t(cpm(raw)),rowGroup=rep,rowID=colnames(raw),method = "Average")
head(t(collapse$datETcollapsed))

#- we want to filter out gene which has low expression level in all sample, and 
#- we set a subjective value `lowest_cpm` to 2. If one of the group average is 
#- bigger than `lowest_cpm`, we will keep that gene.  
#======
lowest_cpm = 2
#======
keep <- rowSums(t(collapse$datETcollapsed)>=lowest_cpm) >= 1
keep
filtered <- raw[keep, keep.lib.sizes=F] 
nrow(filtered) #13755, original 16604

## making sure there's no NA
table(is.na(filtered$counts))

# 2. Normalization -------------------------------------------------------------

## this is the universe (# of expressed genes) of this data set
genes <- rownames(filtered)

write.table(genes, quote=FALSE, row.names=F, col.names=F,sep = "\t",
            file="./List/Lee-Gelata-tuber.universe-13755geneID.txt")

# Before normalization: convert filtered into a Expression Set
set <- newSeqExpressionSet(as.matrix(filtered), 
                           phenoData = data.frame(
                               condition, rep, row.names = colnames(filtered)))

# After normalization
set1 <- betweenLaneNormalization(set, which="upper")

#RUVs: Estimating the factors of unwanted variation using replicate samples
differences<-makeGroups(rep)
differences

set3 <- RUVs(set1, genes, k=1, differences)

pData(set3)

par(bg="white",cex.axis=0.8,cex.lab = 0.8,mar=c(5.1, 3.9, 3.9,11), xpd=TRUE)

plotRLE(set3, outlone=F,ylim=c(-3,3))

plotPCA(set3,cex=1.5,labels=F,pch=c(2,17)[condition])
legend("right",cex=1,ncol=1,bty="n", border=F,inset=c(-0.5,0),
       legend=levels(condition),
       pch=c(2,17))

# 3. DE analysis ---------------------------------------------------------------
## 3.1 design matrix: ~condition

sample.info

design1 <- model.matrix(~condition+W_1,data=pData(set3)) 
design1

colnames(design1)

y <- DGEList(counts=counts(set3),group=condition)
y <- calcNormFactors(y, method = 'TMM')
y$samples

y<-estimateDisp(y,design1,robust = T)
y$common.dispersion #0.01538511

fit <-glmQLFit(y, design1)

## condition effect
colnames(design1)
qlf<-glmQLFTest(fit,coef=2) # coef=2 "conditionbasal"

top <- topTags(qlf, n=nrow(set3))$table



## DEG criteria

## only FDR: too many genes, so will require logFC FDR:false discovery rates  stattquest!
sum(top$FDR<0.05) #7274
sum(top$FDR<0.01) #5148

## FDR plus logFC
sum(top$FDR<0.05 & abs(top$logFC)>1) #3125
sum(top$FDR<0.01 & abs(top$logFC)>1) #2831

## de1: top$FDR<0.01 & abs(top$logFC)>1
de <- top$FDR<0.01 & abs(top$logFC)>1
table(de)
de.gene<-rownames(top[de,])
de.geneinfo<-top[de,]

write.table(de.gene, quote=F, row.names = F, col.names = F,sep="\t",
            file="./List/Lee-Gelata-tuber.edgeR-FDR01logFC1-2831geneID.txt")

## Generate cpm output
logcpm<-cpm(y,prior.count=1, log=TRUE)
write.table(logcpm, quote=FALSE, row.names=T, col.names=T,sep = "\t",
            file="./List/Lee-Gelata-tuber.universe-13755logcpm.txt")




#4. Enrichment analysis -------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))
#setwd("C:/Users/chiayi/Desktop/Student/Gastrodia/")

library(clusterProfiler)

#4.1 Load data
gene<-read.table("List/Lee-Gelata-tuber.edgeR-FDR01logFC1-2831geneID.txt",header=F)

gene=as.factor(gene[,1])

universe = read.table("List/Lee-Gelata-tuber.universe-13755geneID.txt",header=F)
universe = as.character(universe[,1])

#4.2 Enrichment: KEGG
term2gene <- read.table("enrichment_input/GWHBDNU00000000_KEGG.tsv", 
                        header=F,comment.char = "",sep="\t",quote="")

term2gene=term2gene[,c(2,1)]
names(term2gene)=c("KEGG","geneID")
head(term2gene)

term2name<-read.table("enrichment_input/ko_definition_cleaned.txt",
                      header=F ,comment.char = "",sep="\t",quote="")

names(term2name)=c("KO","Description")
head(term2name)

ego<-enricher(gene=gene, universe=universe,
              TERM2GENE=term2gene, 
              TERM2NAME=term2name)

KEGG.out<-as.data.frame(ego)

write.table(KEGG.out, sep="\t", file="enrichment_output/KEGG_output.tsv")
write.csv(KEGG.out, file="enrichment_output/KEGG_output.csv")


#4.3 Enrichment: GO
term2gene <- read.table("enrichment_input/GWHBDNU00000000_go_two_column.tsv", 
                        header=F,comment.char = "",sep="\t",quote="")
term2gene <- term2gene[, c(2,1)]
names(term2gene)=c("GO","geneID")
head(term2gene)

term2name<-read.table("enrichment_input/GO_term_name_type_all_cleaned_2022_07.tsv",
                      header=F,comment.char = "",sep="\t",quote="")

names(term2name)=c("GO","Description", "type")
head(term2name)

GO.ego<-enricher(gene=gene, universe = universe,
                 TERM2GENE = term2gene, 
                 TERM2NAME = term2name[, c(1, 2)])

GO.out<-as.data.frame(GO.ego)
write.table(GO.out, sep="\t", file="enrichment_output/GO_output.tsv")
write.csv(GO.out, file="enrichment_output/GO_output.csv")

