# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2021
# Email:  nitishimtech@gmail.com
# Date: 19 October, 2021
# Script Name: DEG_Analysis_EMT.R
#####################################################################
# #################  Script Description: ############################
#
################## Notes: ###########################################
# This code for DEG analysis
#####################################################################
# SET WORKING DIRECTORY -----------------------------
getwd()
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")

#####################################################################
#####################################################################
############# Source code for basic R functions #####################
#####################################################################

suppressMessages(suppressWarnings(source("C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Common_Functions.R")))

#####################################################################
########################### required library ########################
#####################################################################
library(tximport)
library(DESeq2); library(limma)
library(dplyr); library(tidyverse)
library(msigdbr); library(clusterProfiler)
library(enrichplot)

#####################################################################
#####################################################################
################ Load GTF and selected GFT files ####################
#####################################################################
load("C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/gencode.vM27.gtf.selected.rds")
load("C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/gencode.vM27.gtf.rds")


Duplicated.genes <- gencode.vM27.gtf.selected %>%
  group_by(Symbol) %>%
  filter(n()>1) %>% summarize(n=n())

#####################################################################
#####################################################################
################### Make RSEM Gene Counts Matrix ####################
#####################################################################

## Option1 :: Read RSEM file
#STAR_Counts <- STAR_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/STAR-COUNTS/", Stranded = "reverse")
#RSEM_Counts <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "expected_count")

## Option2 :: Read RSEM file
files = list.files(path = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/",pattern = "*.genes.results")
names(files) <- substr(files, 1, 15)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$counts <- round(txi.rsem$counts)

## Option3 :: Read RSEM file
## Execute RSEM_AWK.sh code and read Expected counts/TPM/FPKM data
#system('sh RSEM_AWK.sh') #This command will generate Counts/TPM/FPKM files in current directory
#RSEM_ExpCounts <- read.table("RSEM_Expected_Counts.tsv", header=TRUE, row.names=1)


#####################################################################
#####################################################################
################### Filtering Step/remove probes ####################
#####################################################################

# Remove all gene with length zero
# https://bioinformatics.stackexchange.com/questions/13521/deseqdatasetfromtximport-alllengths-0-is-not-true
txi.rsem$abundance <- txi.rsem$abundance[apply(txi.rsem$length, 1, function(row) all(row !=0 )),]
txi.rsem$counts <- txi.rsem$counts[apply(txi.rsem$length, 1, function(row) all(row !=0 )),]
txi.rsem$length <- txi.rsem$length[apply(txi.rsem$length, 1, function(row) all(row !=0 )),]

##Remove genes which have over 25% zero
txi.rsem$counts <- as.matrix(as.data.frame(Gene_filtering(txi.rsem$counts)))
txi.rsem$abundance <- txi.rsem$abundance[rownames(txi.rsem$counts),]
txi.rsem$length <- txi.rsem$length[rownames(txi.rsem$counts),]


#####################################################################
#####################################################################
#################### Parameters for Analysis  #######################
FDR.thresh=0.05; th_log2fc=log2(1.2)

#####################################################################
#####################################################################
################## DEG analysis by using DESeq2 #####################
#####################################################################

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
CtrlVsTGFb <- resdata1 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR.thresh & abs(log2FoldChange) >= th_log2fc)

CtrlVsTGFb %>%
  group_by(Symbol) %>%
  filter(n()>1) %>% summarize(n=n()) # Septin2 is duplicated
CtrlVsTGFb_Exp <- CtrlVsTGFb[!duplicated(CtrlVsTGFb$Symbol),]

# ### TGFbeta+CX5461 Vs Control #### This analysis don't required.
# res2 <- DESeq2::results(dds1, contrast=c("condition","TGFb_CX5461","Control"))
# DESeq2.res2 <- as.data.frame(res2)
# DESeq2.res2$ENSG <- rownames(DESeq2.res2)
# resdata2 <- left_join(DESeq2.res2,gencode.vM27.gtf.selected, by ="ENSG")
# resdata2 <- resdata2 %>% 
#   select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
# CtrlVsCX5461 <- resdata2 %>%
#   filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR.thresh & abs(log2FoldChange) >= th_log2fc)


### TGFbeta+CX5461 Vs. TGFbeta
res3 <- DESeq2::results(dds1, contrast=c("condition","TGFb_CX5461","TGFB"))
DESeq2.res3 <- as.data.frame(res3)
DESeq2.res3$ENSG <- rownames(DESeq2.res3)
resdata3 <- left_join(DESeq2.res3,gencode.vM27.gtf.selected, by ="ENSG")
resdata3 <- resdata3 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
TGFbVsCX5461 <- resdata3 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR.thresh & abs(log2FoldChange) >= th_log2fc)
TGFbVsCX5461 %>%
  group_by(Symbol) %>%
  filter(n()>1) %>% summarize(n=n()) # Septin2 is duplicated
TGFbVsCX5461_Exp <- TGFbVsCX5461[!duplicated(TGFbVsCX5461$Symbol),]


#####################################################################
#####################################################################
########################### DEG table ###############################
#####################################################################

colnames(resdata1)[5:10] <- 
  paste0(colnames(resdata1[,5:10]), "_TGFB_Vs_","Control")
colnames(resdata3)[5:10] <- 
  paste0(colnames(resdata3[,5:10]), "_TGFb_CX5461_Vs_","TGFB")

resdata.merge <- merge(resdata1, resdata3, all=TRUE)
rownames(resdata.merge) <- resdata.merge$ENSG
resdata.merge <- resdata.merge[resdata.merge$GeneType=="protein_coding",]

resdata.merge <-  resdata.merge[!duplicated(resdata.merge$Symbol),]


TGFB_Vs_Ctrl_Exp.status <- ifelse(resdata.merge$log2FoldChange_TGFB_Vs_Control >= th_log2fc & 
                                resdata.merge$padj_TGFB_Vs_Control <= FDR.thresh, "up", 
                              ifelse(resdata.merge$log2FoldChange_TGFB_Vs_Control <= -th_log2fc &
                                       resdata.merge$padj_TGFB_Vs_Control <= FDR.thresh, "down", "notSig"))

CX5461_Vs_TGFB_Exp.status <- ifelse(resdata.merge$log2FoldChange_TGFb_CX5461_Vs_TGFB >= th_log2fc & 
                                  resdata.merge$padj_TGFb_CX5461_Vs_TGFB <= FDR.thresh, "up", 
                                ifelse(resdata.merge$log2FoldChange_TGFb_CX5461_Vs_TGFB <= -th_log2fc &
                                         resdata.merge$padj_TGFb_CX5461_Vs_TGFB <= FDR.thresh, "down", "notSig"))
Exp_direction <- ifelse(TGFB_Vs_Ctrl_Exp.status=="up" & CX5461_Vs_TGFB_Exp.status=="up", "upUp", 
                    ifelse(TGFB_Vs_Ctrl_Exp.status=="down" & CX5461_Vs_TGFB_Exp.status=="down", "downDown", 
                           ifelse(TGFB_Vs_Ctrl_Exp.status=="up" & CX5461_Vs_TGFB_Exp.status=="down", "upDown",
                                  ifelse(TGFB_Vs_Ctrl_Exp.status=="down" & CX5461_Vs_TGFB_Exp.status=="up", "downUp", "notSig"))))

resdata.merge$TGFB_Vs_Ctrl_Exp.status   <- TGFB_Vs_Ctrl_Exp.status
resdata.merge$CX5461_Vs_TGFB_Exp.status <- CX5461_Vs_TGFB_Exp.status
resdata.merge$Exp_direction <- Exp_direction

#write.csv(resdata.merge, "result_Merged.csv")


#####################################################################
#####################################################################
############# Differential translation analysis #####################
#####################################################################

Ribo_STAR_genes <- STAR_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/STAR-Counts/", Stranded = "unstranded")
## For HTSeq-counts and featureCounts I have to decide strandness in advanced. They will provide only one read counts. I am using "reverse" by default, i real life "reverse" is still "unstranded" read counts.
Ribo_HTSeq_Counts_genes <- HTSeq_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/HTSeq-Counts/Gene/", Stranded = "reverse", counts = "genes")
Ribo_FeatureCounts_genes <- FeatureCounts_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/FeatureCounts/Gene/", Stranded = "reverse", counts = "genes")

sampleTable <- data.frame(condition = factor(rep(c("Control", "TGFB", "TGFb_CX5461"), each = 3)))
### Set the name of input file here
matrix_input_Ribo <- Ribo_FeatureCounts_genes
rownames(sampleTable) <- colnames(matrix_input_Ribo)
dds.Ribo <- DESeq2::DESeqDataSetFromMatrix(matrix_input_Ribo, sampleTable, ~condition)
keep <- rowSums(counts(dds.Ribo)) >= 2
dds.Ribo <- dds.Ribo[keep,]

dds1.ribo <- DESeq(dds.Ribo, test="LRT", reduced=~1) ## LRT giving more DEGs


### TGFbeta Vs. Control analysis
res1.ribo <- DESeq2::results(dds1.ribo, contrast=c("condition","TGFB", "Control"))
DESeq2.res1.ribo <- as.data.frame(res1.ribo)
DESeq2.res1.ribo$ENSG <- rownames(DESeq2.res1.ribo)
resdata1.ribo <- left_join(DESeq2.res1.ribo,gencode.vM27.gtf.selected, by ="ENSG")
resdata1.ribo <- resdata1.ribo %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
Ribo_CtrlVsTGFb <- resdata1.ribo %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR.thresh & abs(log2FoldChange) >= th_log2fc)

Ribo_CtrlVsTGFb %>%
  group_by(Symbol) %>%
  filter(n()>1) %>% summarize(n=n()) # Septin2 is duplicated
Ribo_CtrlVsTGFb <- Ribo_CtrlVsTGFb[!duplicated(Ribo_CtrlVsTGFb$Symbol),]


### TGFbeta+CX5461 Vs. TGFbeta
res3 <- DESeq2::results(dds1.ribo, contrast=c("condition","TGFb_CX5461","TGFB"))
DESeq2.res3.ribo <- as.data.frame(res3)
DESeq2.res3.ribo$ENSG <- rownames(DESeq2.res3.ribo)
resdata3.ribo <- left_join(DESeq2.res3.ribo,gencode.vM27.gtf.selected, by ="ENSG")
resdata3.ribo <- resdata3.ribo %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
Ribo_TGFbVsCX5461 <- resdata3.ribo %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR.thresh & abs(log2FoldChange) >= th_log2fc)
Ribo_TGFbVsCX5461 %>%
  group_by(Symbol) %>%
  filter(n()>1) %>% summarize(n=n()) # Septin2 is duplicated
Ribo_TGFbVsCX5461_Exp <- Ribo_TGFbVsCX5461[!duplicated(Ribo_TGFbVsCX5461$Symbol),]


#####################################################################
#####################################################################
############################# DEG table #############################
#####################################################################

colnames(resdata1.ribo)[5:10] <- 
  paste0(colnames(resdata1.ribo[,5:10]), "_TGFB_Vs_","Control_Ribo")
colnames(resdata3.ribo)[5:10] <- 
  paste0(colnames(resdata3.ribo[,5:10]), "_TGFb_CX5461_Vs_","TGFB_Ribo")

resdata.merge.Ribo <- merge(resdata1.ribo, resdata3.ribo, all=TRUE)
rownames(resdata.merge.Ribo) <- resdata.merge.Ribo$ENSG
resdata.merge.Ribo <- resdata.merge.Ribo[resdata.merge.Ribo$GeneType=="protein_coding",]

resdata.merge.Ribo <-  resdata.merge.Ribo[!duplicated(resdata.merge.Ribo$Symbol),]


TGFB_Vs_Ctrl_Ribo.status <- ifelse(resdata.merge.Ribo$log2FoldChange_TGFB_Vs_Control_Ribo >= th_log2fc & 
                                    resdata.merge.Ribo$padj_TGFB_Vs_Control_Ribo <= FDR.thresh, "up", 
                                  ifelse(resdata.merge.Ribo$log2FoldChange_TGFB_Vs_Control_Ribo <= -th_log2fc &
                                           resdata.merge.Ribo$padj_TGFB_Vs_Control_Ribo <= FDR.thresh, "down", "notSig"))

CX5461_Vs_TGFB_Ribo.status <- ifelse(resdata.merge.Ribo$log2FoldChange_TGFb_CX5461_Vs_TGFB_Ribo >= th_log2fc & 
                                      resdata.merge.Ribo$padj_TGFb_CX5461_Vs_TGFB_Ribo <= FDR.thresh, "up", 
                                    ifelse(resdata.merge.Ribo$log2FoldChange_TGFb_CX5461_Vs_TGFB_Ribo <= -th_log2fc &
                                             resdata.merge.Ribo$padj_TGFb_CX5461_Vs_TGFB_Ribo <= FDR.thresh, "down", "notSig"))

Ribo_direction <- ifelse(TGFB_Vs_Ctrl_Ribo.status=="up" & CX5461_Vs_TGFB_Ribo.status=="up", "upUp", 
                        ifelse(TGFB_Vs_Ctrl_Ribo.status=="down" & CX5461_Vs_TGFB_Ribo.status=="down", "downDown", 
                               ifelse(TGFB_Vs_Ctrl_Ribo.status=="up" & CX5461_Vs_TGFB_Ribo.status=="down", "upDown",
                                      ifelse(TGFB_Vs_Ctrl_Ribo.status=="down" & CX5461_Vs_TGFB_Ribo.status=="up", "downUp", "notSig"))))

resdata.merge.Ribo$TGFB_Vs_Ctrl_Ribo.status <- TGFB_Vs_Ctrl_Ribo.status
resdata.merge.Ribo$CX5461_Vs_TGFB_Ribo.status <- CX5461_Vs_TGFB_Ribo.status
resdata.merge.Ribo$Ribo_direction <- Ribo_direction

#####################################################################
#####################################################################
########### Merge RNAseq and Riboprofiling final results ############
#####################################################################
Results <- merge(resdata.merge,resdata.merge.Ribo, all=TRUE)


#####################################################################
#####################################################################
################### Data for Star Burst plot of DE ##################
#####################################################################

############################ mRNA ###################################
list_genes <- list()
rownames(Results) <- Results$Symbol

#mRNA_up
sym_mrna_up <- rownames(Results[grep("up", Results$TGFB_Vs_Ctrl_Exp.status),])
list_genes[['mrna_up']] <- sym_mrna_up

#mRNA_down
sym_mrna_dn <- rownames(Results[grep("down", Results$TGFB_Vs_Ctrl_Exp.status),])
list_genes[['mrna_dn']] <- sym_mrna_dn

##mRNA_total DE
sym_mrna <- union(sym_mrna_up, sym_mrna_dn)
length(sym_mrna)

#mRNA_no
df_mrna <- resdata.merge %>%
   filter(stringr::str_detect(GeneType, "protein_coding"))
sym_rnaseq <- df_mrna$Symbol
sym_mrna_no <- setdiff(sym_rnaseq, sym_mrna)
length(sym_mrna_no)


######################### Riboprofiling #############################
#Ribo_up
sym_ribo_up <- rownames(Results[grep("up", Results$TGFB_Vs_Ctrl_Ribo.status),])
list_genes[['ribo_up']] <- sym_ribo_up

#Ribo_down
sym_ribo_dn <- rownames(Results[grep("down", Results$TGFB_Vs_Ctrl_Ribo.status),])
list_genes[['ribo_dn']] <- sym_ribo_dn

##Ribo_total DE
sym_ribo <- union(sym_ribo_up, sym_ribo_dn)
length(sym_ribo)

#Ribo_no
df_ribo <- resdata.merge.Ribo %>%
  filter(stringr::str_detect(GeneType, "protein_coding"))
sym_riboseq <- df_ribo$Symbol
sym_ribo_no <- setdiff(sym_riboseq, sym_ribo)
length(sym_mrna_no)


####################### Data For Display ############################
sym_mrna_only <- setdiff(sym_mrna, sym_ribo)
sym_mrna_only <- setdiff(sym_mrna_only, excludes)
sym_mrna_up_ribo_up <- setdiff(sym_mrna_up_ribo_up, excludes)
sym_mrna_dn_ribo_dn <- setdiff(sym_mrna_dn_ribo_dn, excludes)
sym_mrna_no_ribo_up <- setdiff(sym_mrna_no_ribo_up, excludes)
sym_mrna_no_ribo_dn <- setdiff(sym_mrna_no_ribo_dn, excludes)
sym_mrna_up_ribo_dn <- setdiff(sym_mrna_up_ribo_dn, excludes)
sym_mrna_dn_ribo_up <- setdiff(sym_mrna_dn_ribo_up, excludes)


#####################################################################
#####################################################################
####################### GSEA/Pathway analysis #######################
####################### For Transcription ###########################
#####################################################################

########## Control Vs. TFGbeta pathway analysis ###########
testdata=resdata1 ## Change resdata1 with other results
nv <- sign(testdata$log2FoldChange)*(-log10(testdata$padj))
names(nv) <- testdata$Symbol;
f <- is.finite(nv); t <- min(nv[f]); nv[nv < t] <- t*1.1; t <- max(nv[f]); nv[nv > t] <- t*1.1 #From Hyunsoo Kim fig_riboseq_gsea.ipynb

gmt_mus <- msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, gene_symbol)

### 50 Hallmark pathway analysis; http://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H
h_gene_sets = msigdbr(species = "mouse", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

set.seed(1234)
gsea_all_EXP_12 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                 nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 TERM2GENE=gmt_mus, TERM2NAME=NA, seed=FALSE, by="fgsea")

head(gsea_all_EXP_12); dim(gsea_all_EXP_12)

gsea_h_gene_sets_eXP_12 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                         nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         TERM2GENE=h_gene_sets, TERM2NAME=NA, seed=FALSE, by="fgsea")
head(gsea_h_gene_sets_eXP_12); dim(gsea_h_gene_sets_eXP_12)

gsea_h_gene_sets_eXP_12@result %>% 
  dplyr::select(enrichmentScore,  NES, p.adjust) %>% 
  arrange(desc(enrichmentScore)) %>% 
  write.csv(., file = "eXPRESSION Control Vs TGFBeta.csv")

#####################################################################
##################### Plot GSEA analysis output #####################
p1 <- gseaplot(gsea_h_gene_sets_eXP_12, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_eXP_12$Description[1])
p2 <- gseaplot(gsea_h_gene_sets_eXP_12, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_eXP_12$Description[1])
p3 <- gseaplot(gsea_h_gene_sets_eXP_12, geneSetID = 1, title = gsea_h_gene_sets_eXP_12$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

p1 <- gseaplot(gsea_h_gene_sets_eXP_12, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_eXP_12$Description[2])
p2 <- gseaplot(gsea_h_gene_sets_eXP_12, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_eXP_12$Description[2])
p3 <- gseaplot(gsea_h_gene_sets_eXP_12, geneSetID = 1, title = gsea_h_gene_sets_eXP_12$Description[2])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

enrichplot::gseaplot2(gsea_h_gene_sets_eXP_12, geneSetID = 1, title = gsea_h_gene_sets_eXP_12$Description[1])
enrichplot::gseaplot2(gsea_h_gene_sets_eXP_12, geneSetID = 8, title = gsea_h_gene_sets_eXP_12$Description[8])

#####################################################################
###################### GSEA/Pathway analysis ########################
###################### CX5461 Vs TFGBeta ############################

testdata=resdata3 ## Change resdata1 with other analysis results
nv <- sign(testdata$log2FoldChange)*(-log10(testdata$padj))
names(nv) <- testdata$Symbol
f <- is.finite(nv); t <- min(nv[f]); nv[nv < t] <- t*1.1; t <- max(nv[f]); nv[nv > t] <- t*1.1 #From Hyunsoo Kim fig_riboseq_gsea.ipynb

set.seed(1234)
gsea_all_Exp_23 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                 nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 TERM2GENE=gmt_mus, TERM2NAME=NA, seed=FALSE, by="fgsea")

head(gsea_all_Exp_23); dim(gsea_all_Exp_23)

gsea_h_gene_sets_Exp_23 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                         nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         TERM2GENE=h_gene_sets, TERM2NAME=NA, seed=FALSE, by="fgsea")
head(gsea_h_gene_sets_Exp_23); dim(gsea_h_gene_sets_Exp_23)

gsea_h_gene_sets_Exp_23@result %>% 
  dplyr::select(enrichmentScore,  NES, p.adjust) %>% 
  arrange(desc(enrichmentScore)) %>% 
  write.csv(., file = "Expression TGFBeta Vs CX5461.csv")

#####################################################################
##################### Plot GSEA analysis output #####################
p1 <- gseaplot(gsea_h_gene_sets_Exp_23, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_Exp_23$Description[1])
p2 <- gseaplot(gsea_h_gene_sets_Exp_23, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_Exp_23$Description[1])
p3 <- gseaplot(gsea_h_gene_sets_Exp_23, geneSetID = 1, title = gsea_h_gene_sets_Exp_23$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

p1 <- gseaplot(gsea_h_gene_sets_Exp_23, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_Exp_23$Description[2])
p2 <- gseaplot(gsea_h_gene_sets_Exp_23, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_Exp_23$Description[2])
p3 <- gseaplot(gsea_h_gene_sets_Exp_23, geneSetID = 1, title = gsea_h_gene_sets_Exp_23$Description[2])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

enrichplot::gseaplot2(gsea_h_gene_sets_Exp_23, geneSetID = 1, title = gsea_h_gene_sets_Exp_23$Description[1])

#####################################################################
#####################################################################
####################### GSEA/Pathway analysis #######################
####################### For Translation #############################
#####################################################################

########## Control Vs. TFGbeta pathway analysis ###########
testdata=resdata1.ribo ## Change resdata1 with other results
nv <- sign(testdata$log2FoldChange)*(-log10(testdata$padj))
names(nv) <- testdata$Symbol;
f <- is.finite(nv); t <- min(nv[f]); nv[nv < t] <- t*1.1; t <- max(nv[f]); nv[nv > t] <- t*1.1 #From Hyunsoo Kim fig_riboseq_gsea.ipynb


set.seed(1234)
gsea_all_Ribo_12 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                 nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 TERM2GENE=gmt_mus, TERM2NAME=NA, seed=FALSE, by="fgsea")

head(gsea_all_Ribo_12); dim(gsea_all_Ribo_12)

gsea_h_gene_sets_Ribo_12 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                         nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         TERM2GENE=h_gene_sets, TERM2NAME=NA, seed=FALSE, by="fgsea")
head(gsea_h_gene_sets_Ribo_12); dim(gsea_h_gene_sets_Ribo_12)

gsea_h_gene_sets_Ribo_12@result %>% 
  dplyr::select(enrichmentScore,  NES, p.adjust) %>% 
  arrange(desc(enrichmentScore)) %>% 
  write.csv(., file = "Riboprofiling Control Vs TGFBeta.csv")

#####################################################################
##################### Plot GSEA analysis output #####################
p1 <- gseaplot(gsea_h_gene_sets_Ribo_12, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_Ribo_12$Description[1])
p2 <- gseaplot(gsea_h_gene_sets_Ribo_12, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_Ribo_12$Description[1])
p3 <- gseaplot(gsea_h_gene_sets_Ribo_12, geneSetID = 1, title = gsea_h_gene_sets_Ribo_12$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

p1 <- gseaplot(gsea_h_gene_sets_Ribo_12, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_Ribo_12$Description[2])
p2 <- gseaplot(gsea_h_gene_sets_Ribo_12, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_Ribo_12$Description[2])
p3 <- gseaplot(gsea_h_gene_sets_Ribo_12, geneSetID = 1, title = gsea_h_gene_sets_Ribo_12$Description[2])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

enrichplot::gseaplot2(gsea_h_gene_sets_Ribo_12, geneSetID = 1, title = gsea_all$Description[2])

#####################################################################
###################### GSEA/Pathway analysis ########################
###################### CX5461 Vs TFGBeta ############################

testdata=resdata3.ribo ## Change resdata1 with other analysis results
nv <- sign(testdata$log2FoldChange)*(-log10(testdata$padj))
names(nv) <- testdata$Symbol
f <- is.finite(nv); t <- min(nv[f]); nv[nv < t] <- t*1.1; t <- max(nv[f]); nv[nv > t] <- t*1.1 #From Hyunsoo Kim fig_riboseq_gsea.ipynb


set.seed(1234)
gsea_all_Ribo_23 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                 nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                 TERM2GENE=gmt_mus, TERM2NAME=NA, seed=FALSE, by="fgsea")

head(gsea_all_Ribo_23); dim(gsea_all_Ribo_23)

gsea_h_gene_sets_Ribo_23 <- GSEA(sort(nv[!is.na(nv)], decreasing=T), exponent =1,
                         nPerm = 10000, minGSSize = 5, maxGSSize = 500,   
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                         TERM2GENE=h_gene_sets, TERM2NAME=NA, seed=FALSE, by="fgsea")
head(gsea_h_gene_sets_Ribo_23); dim(gsea_h_gene_sets_Ribo_23)

gsea_h_gene_sets_Ribo_23@result %>% 
  dplyr::select(enrichmentScore,  NES, p.adjust) %>% 
  arrange(desc(enrichmentScore)) %>% 
  write.csv(., file = "Riboprofiling TGFBeta Vs CX5461.csv")

#####################################################################
##################### Plot GSEA analysis output #####################
p1 <- gseaplot(gsea_h_gene_sets_Ribo_23, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_Ribo_23$Description[1])
p2 <- gseaplot(gsea_h_gene_sets_Ribo_23, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_Ribo_23$Description[1])
p3 <- gseaplot(gsea_h_gene_sets_Ribo_23, geneSetID = 1, title = gsea_h_gene_sets_Ribo_23$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])

p1 <- gseaplot(gsea_h_gene_sets_Ribo_23, geneSetID = 1, by = "runningScore", title = gsea_h_gene_sets_Ribo_23$Description[2])
p2 <- gseaplot(gsea_h_gene_sets_Ribo_23, geneSetID = 1, by = "preranked", title = gsea_h_gene_sets_Ribo_23$Description[2])
p3 <- gseaplot(gsea_h_gene_sets_Ribo_23, geneSetID = 1, title = gsea_h_gene_sets_Ribo_23$Description[2])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])


#####################################################################
#####################################################################
######################### Save workspace ############################
save.image("../../DESeq2_LRT.RData")

#####################################################################
#####################################################################
########################### Analysis done ###########################