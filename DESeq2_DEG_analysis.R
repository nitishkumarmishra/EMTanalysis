# ###################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2021
# Email:  nitishimtech@gmail.com
# Date: 2021-09-22
# Script Name: code_name.R
#####################################################
# #############  Script Description: ################
#
############### Notes: ##############################
#
#####################################################
########## SET WORKING DIRECTORY ####################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################
#setwd("C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts")

###### Read RSEM output files #####
library(tximport)
files = list.files(pattern = "*.genes.results")
names(files) <- substr(files, 1, 15)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$counts <- round(txi.rsem$counts)
write.table(txi.rsem$counts, file="RSEM_Counts.tsv", sep="\t", quote=FALSE)


###### Make GTF for annotation #####
library('rtracklayer')
GTF="C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/gencode.vM27.annotation.gtf"
GTF <- import(GTF)
gencode.vM27.gtf <- as.data.frame(GTF)
library(dplyr)
gencode.vM27.gtf.selected <- gencode.vM27.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type,ENSG=gene_id, GeneType=gene_type, Symbol=gene_name) %>%
  select(Chr, Type, GeneType, ENSG, Symbol) #%>%
#  distinct(Symbol, .keep_all= TRUE) ## Remove duplicated, check at final stage.

### Count duplicated genes, just to explore.
Duplicated.genes <- gencode.vM27.gtf.selected %>%
  group_by(Symbol) %>%
  filter(n()>1) %>% summarize(n=n())

## Execute RSEM_AWK.sh code and read Expected counts/TPM/FPKM data
#system('sh RSEM_AWK.sh') #This command will generate Counts/TPM/FPKM files in current directory
#RSEM_ExpCounts <- read.table("RSEM_Expected_Counts.tsv", header=TRUE, row.names=1)


# Remove all gene with length zero
# https://bioinformatics.stackexchange.com/questions/13521/deseqdatasetfromtximport-alllengths-0-is-not-true
txi.rsem$abundance <- txi.rsem$abundance[apply(txi.rsem$length, 1, function(row) all(row !=0 )),]
txi.rsem$counts <- txi.rsem$counts[apply(txi.rsem$length, 1, function(row) all(row !=0 )),]
txi.rsem$length <- txi.rsem$length[apply(txi.rsem$length, 1, function(row) all(row !=0 )),]


################################################################
################################################################
################################################################
########### Filtering Step/remove probes #############
Gene_filtering <- function(mat, Count_filter=0){
  not_all_na <- function(x) any(!is.na(x))
  if(!isTRUE(is.data.frame(mat))) # If not dataFrame convert it.
    mat <- as.data.frame(mat)
  mat <- mat %>%
    dplyr::select(where(not_all_na))
  mat <- subset(mat, apply(mat, 1, sum) != 0)
  #mat <- subset(mat, apply(mat, 2, sum) != 0)# not_all_na already did this
  mat <- mat[apply(mat,1,function(x) sum(x<=0))< ncol(mat)*0.75,]##Remove genes which have over 25% zero
  mat <- subset(mat, apply(mat, 1, sum) >= Count_filter)
  return(mat)}

txi.rsem$counts <- as.matrix(Gene_filtering(txi.rsem$counts))
txi.rsem$abundance <- txi.rsem$abundance[rownames(txi.rsem$counts),]
txi.rsem$length <- txi.rsem$length[rownames(txi.rsem$counts),]

########### DEG analysis by using DESeq2 #########
library(DESeq2)
sampleTable <- data.frame(condition = factor(rep(c("Control", "TGFB", "TGFb_CX5461"), each = 3)))
rownames(sampleTable) <- colnames(txi.rsem$counts)
dds <- DESeq2::DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= 3
dds <- dds[keep,]

# dds is now ready for DESeq() see DESeq2 vignette
#dds <- DESeq2::DESeq(dds)
#res <- DESeq2::results(dds)
#res <- DESeq2::results(dds, contrast=c("condition","TGFb_CX5461","TGFB"))
dds1 <- DESeq(dds, test="LRT", reduced=~1) ## LRT giving more DEGs

### TGFbeta Vs. Control analysis
res1 <- DESeq2::results(dds1, contrast=c("condition","TGFB", "Control"))
DESeq2.res1 <- as.data.frame(res1)
DESeq2.res1$ENSG <- rownames(DESeq2.res1)
resdata1 <- left_join(DESeq2.res1,gencode.vM27.gtf.selected, by ="ENSG")
resdata1 <- resdata1 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
resdata1_adjP_0.05 <- resdata1 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= 0.05 & abs(log2FoldChange) >= 1)


### TGFbeta+CX5461 Vs Control #### This analysis don't required.
res2 <- DESeq2::results(dds1, contrast=c("condition","TGFb_CX5461","Control"))
DESeq2.res2 <- as.data.frame(res2)
DESeq2.res2$ENSG <- rownames(DESeq2.res2)
resdata2 <- left_join(DESeq2.res2,gencode.vM27.gtf.selected, by ="ENSG")
resdata2 <- resdata2 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
resdata2_adjP_0.05 <- resdata2 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= 0.05 & abs(log2FoldChange) >= 1)


### TGFbeta+CX5461 Vs. TGFbeta
res3 <- DESeq2::results(dds1, contrast=c("condition","TGFb_CX5461","TGFB"))
DESeq2.res3 <- as.data.frame(res3)
DESeq2.res3$ENSG <- rownames(DESeq2.res3)
resdata3 <- left_join(DESeq2.res3,gencode.vM27.gtf.selected, by ="ENSG")
resdata3 <- resdata3 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
resdata3_adjP_0.05 <- resdata3 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= 0.05 & abs(log2FoldChange) >= 1)


###########################################################
###########################################################
###########################################################
########## Save workspace ############
save.image("DESeq2_LRT.RData")

###########################################################
###########################################################
###########################################################
######################## DEG table ########################


colnames(resdata1)[5:10] <- 
  paste0(colnames(resdata1[,5:10]), "_TGFB_Vs_","Control")
colnames(resdata3)[5:10] <- 
  paste0(colnames(resdata3[,5:10]), "_TGFb_CX5461_Vs_","TGFB")

resdata.merge <- merge(resdata1, resdata3, all=TRUE)
rownames(resdata.merge) <- resdata.merge$ENSG
resdata.merge <- resdata.merge[resdata.merge$GeneType=="protein_coding",]

FDR.thresh=0.05; th_log2fc=log(1.2)


# logi.sig  <-  resdata.merge$padj_TGFB_Vs_Control <= FDR.thresh
# logi.up  <-  df$log2FC > th_log2fc
# logi.dn  <-  df$log2FC < -th_log2fc
# df$direc  <-  "notSig"
# df$direc[logi.sig & logi.up]  <-  "up"
# df$direc[logi.sig & logi.dn]  <-  "down"

TGFB_Vs_Ctrl.status <- ifelse(resdata.merge$log2FoldChange_TGFB_Vs_Control >= th_log2fc & 
         resdata.merge$padj_TGFB_Vs_Control <= FDR.thresh, "up", 
       ifelse(resdata.merge$log2FoldChange_TGFB_Vs_Control <= -th_log2fc &
                resdata.merge$padj_TGFB_Vs_Control <= FDR.thresh, "down", "notSig"))

CX5461_Vs_TGFB.status <- ifelse(resdata.merge$log2FoldChange_TGFb_CX5461_Vs_TGFB >= th_log2fc & 
                                resdata.merge$padj_TGFb_CX5461_Vs_TGFB <= FDR.thresh, "up", 
                              ifelse(resdata.merge$log2FoldChange_TGFb_CX5461_Vs_TGFB <= -th_log2fc &
                                       resdata.merge$padj_TGFb_CX5461_Vs_TGFB <= FDR.thresh, "down", "notSig"))
direction <- ifelse(TGFB_Vs_Ctrl.status=="up" & CX5461_Vs_TGFB.status=="up", "upUp", 
       ifelse(TGFB_Vs_Ctrl.status=="down" & CX5461_Vs_TGFB.status=="down", "downDown", 
              ifelse(TGFB_Vs_Ctrl.status=="up" & CX5461_Vs_TGFB.status=="down", "upDown",
                     ifelse(TGFB_Vs_Ctrl.status=="down" & CX5461_Vs_TGFB.status=="up", "downUp", "notSig"))))

resdata.merge$TGFB_Vs_Ctrl.status <- TGFB_Vs_Ctrl.status
resdata.merge$CX5461_Vs_TGFB.status <- CX5461_Vs_TGFB.status
resdata.merge$direction <- direction

write.csv(resdata.merge, "result_Merged.csv")

# fdr.12  <-  paste( "FDR" ,  conditions[1] , "VS" , conditions[2] , sep="" )
# fdr.23  <-  paste( "FDR" ,  conditions[2] , "VS" , conditions[3] , sep="" )
# l2fc.12  <-  paste( "log2FC" ,  conditions[1] , "VS" , conditions[2] , sep="" )
# l2fc.23  <-  paste( "log2FC" ,  conditions[2] , "VS" , conditions[3] , sep="" )
# 
# logi.sig2  <-  df[[fdr.12]] <= FDR.thresh  &  df[[fdr.23]] <= FDR.thresh
# 
# logi.upDown  <-  df[[l2fc.12]] > th_log2fc  &  df[[l2fc.23]] < -th_log2fc
# logi.downUp  <-  df[[l2fc.12]] < -th_log2fc  &  df[[l2fc.23]] > th_log2fc
# 
# df$direc  <-  "notSig"
# df$direc[logi.sig2 & logi.upDown]  <-  "upDown"
# df$direc[logi.sig2 & logi.downUp]  <-  "downUp"


###########################################################
###########################################################
###########################################################
##### Just to check S values in in lfshrinkage ########
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#altshrink
#https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#alternative-shrinkage-estimators
library(apeglm)
sampleTable <- data.frame(condition = factor(rep(c("Sample", "TGFB", "TGFb_CX5461"), each = 3)))
rownames(sampleTable) <- colnames(txi.rsem$counts)
sampleTable$condition <- relevel(sampleTable$condition, "TGFB")
dds <- DESeq2::DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
dds <- DESeq2::DESeq(dds, test="LRT", reduced=~1)
resApeT <- lfcShrink(dds, coef="condition_TGFb_CX5461_vs_TGFB", type="apeglm", lfcThreshold=1) 
plotMA(resApeT, ylim=c(-3,3), cex=.8)
abline(h=c(-1,1), col="dodgerblue", lwd=2)


###########################################################
###########################################################
###########################################################
########## GSEA Analysis ###########
# https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html?q=signature#msigdb-analysis
# https://github.com/nitishkumarmishra/PDAC/blob/main/TCGA_Ductal_Clinical.R
# Filter GSEA by using gage package
# https://bioinformaticsbreakdown.com/how-to-gsea/
# http://crazyhottommy.blogspot.com/2016/08/gene-set-enrichment-analysis-gsea.html
#resdata1 <- na.omit(resdata1)
testdata=resdata3 ## Change resdata1 with other analysis results
nv <- sign(testdata$log2FoldChange)*(-log10(testdata$padj))
names(nv) <- testdata$Symbol;
f <- is.finite(nv); t <- min(nv[f]); nv[nv < t] <- t*1.1; t <- max(nv[f]); nv[nv > t] <- t*1.1 #From Hyunsoo Kim fig_riboseq_gsea.ipynb

library(msigdbr); library(clusterProfiler)
#msigdbr_species() ## Check the list of organism in msigdbr
#msigdbr_collections() ## Check the category in msigdbr

gmt_mus <- msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)

### 50 Hallmark pathway analysis; http://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H
h_gene_sets = msigdbr(species = "mouse", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

set.seed(1234)
gsea_all <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                 nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 TERM2GENE=gmt_mus, TERM2NAME=NA, seed=FALSE, by="fgsea")

head(gsea_all); dim(gsea_all)

gsea_h_gene_sets <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                 nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 TERM2GENE=h_gene_sets, TERM2NAME=NA, seed=FALSE, by="fgsea")
head(gsea_h_gene_sets); dim(gsea_h_gene_sets)

gsea_h_gene_sets@result %>% 
  dplyr::select(enrichmentScore,  NES, p.adjust) %>% 
  arrange(desc(enrichmentScore)) %>% 
  write.csv(., file = "TGFBeta Vs TGFBeta Plus CX5461.csv")


###########################################################
###########################################################
###########################################################
##################### Plot GSEA analysis output ##############
library(enrichplot)
p1 <- gseaplot(gsea_all, geneSetID = 1, by = "runningScore", title = gsea_all$Description[1])
p2 <- gseaplot(gsea_all, geneSetID = 1, by = "preranked", title = gsea_all$Description[1])
p3 <- gseaplot(gsea_all, geneSetID = 1, title = gsea_all$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])


p1 <- gseaplot(gsea_all, geneSetID = 2, by = "runningScore", title = gsea_all$Description[2])
p2 <- gseaplot(gsea_all, geneSetID = 2, by = "preranked", title = gsea_all$Description[2])
p3 <- gseaplot(gsea_all, geneSetID = 2, title = gsea_all$Description[2])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])


enrichplot::gseaplot2(gsea_all, geneSetID = 1, title = gsea_all$Description[1])