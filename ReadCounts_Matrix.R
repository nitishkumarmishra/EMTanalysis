# ###################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2021
# Email:  nitishimtech@gmail.com
# Date: 2021-09-22
# Script Name: code_name.R
#####################################################
# Script Description: R code for differential expression analysis
# Notes:
### This code is to make expression matrix from RSEM/STAR count files
# For recursive I can use https://community.rstudio.com/t/read-table-of-list-files/24004/2
#####################################################
# SET WORKING DIRECTORY -----------------------------
getwd()
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################
### PATH of STAR and RSEM counts files ########
#directory <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/STAR-COUNTS/"
#directory <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/"

############## STAR Count matrix ##############
#This sequencing is based on reverse-stand protocol. Reads are calculated base on reverse strand, so I have to use column 4 of geneCount
# https://www.biostars.org/p/360610/
STAR_Matrix <- function(directory="", Stranded="farward")
  {
  basedir <- directory  
  setwd(basedir)
  files <- dir(pattern="_ReadsPerGene.out.tab", path = basedir)
  counts <- c()
  for( i in seq_along(files) ){
    x <- read.table(file=files[i], sep="\t", header=F, as.is=T, skip = 4, row.names = 1) # First four lines are statistics
    colnames(x) <- c("unstranded", "farward", "reverse")
    counts <- cbind(counts, x[,grep(Stranded, colnames(x))]) ## In STAR column 1) ID 2) Unstranded 3) Farward Strand 4) Reverse Strand
  }
  rownames(counts) <- rownames(x) # set the row names
  colnames(counts) <- sub("_ReadsPerGene.out.tab","",files) # set the column names based on input file names, with pattern removed
  return(counts)
}
STAR_Counts <- STAR_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/STAR-COUNTS/", Stranded = "reverse")


############## RSEM Count matrix ##############
RSEM_Matrix <- function(directory="", count_type="FPKM")
{
  basedir <- directory
  setwd(basedir)
  #count_type="FPKM" # "expected_count"   "TPM"              "FPKM"
  #This sequencing is based on reverse-stand protocol. Reads are calculated base on reverse strand, so I have to use column 4 of geneCount
  # https://www.biostars.org/p/360610/
  
  files <- dir(pattern=".RSEM.genes.results", path = basedir)
  counts <- c()
  for( i in seq_along(files) ){
    x <- read.table(file=files[i], sep="\t", header=T, as.is=T, row.names = 1) # First four lines are statistics
    counts <- cbind(counts, x[,grep(count_type, colnames(x))]) ## In STAR column 1) ID 2) Unstranded 3) Farward Strand 4) Reverse Strand
  }
  rownames(counts) <- rownames(x) # set the row names
  colnames(counts) <- sub(".RSEM.genes.results","",files) # set the column names based on input file names, with pattern removed
  return(counts)
}

RSEM_Counts <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "expected_count")
RSEM_TPM <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "TPM")
RSEM_FPKM <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "FPKM")


#####################################################################
#####################################################################
#####################################################################
############### GTF file, remove duplicated genes ###################
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


###########################################################
###########################################################
###########################################################
## This part when rtracklayer was not working
#gene_lengths <- read.csv("gencode.v37.table.txt", header = TRUE, sep = "\t") 
gene_lengths <- gencode.vM27.gtf %>%
  filter(type=="gene") %>%
  rename(Chromosome= seqnames, Type= type, Geneid = gene_id, Class=gene_type, GeneSymbol= gene_name, Length=width, Start = start, End=end, Strand=strand) %>%
  select(Geneid, GeneSymbol, Chromosome, Start, End, Class,  Strand, Length) %>%
  filter(!grepl("chrM",Chromosome))

rownames(gene_lengths) <- gene_lengths$Geneid
#exp <- STAR_Counts

#####################################################################
#####################################################################
#####################################################################
################### Parameters for Analysis  ########################
FoldChange = log2(1.2)
FDR = 0.05

#####################################################################
#####################################################################
#####################################################################
################### Limma+Voom DEG Analysis  ########################

library(edgeR)
library(limma)

# define group
groupLabels <- rep(c("Control", "TGFB", "CX5461"), each = 3)
TS <- factor(groupLabels, levels = c("Control", "TGFB", "CX5461"))

design <- model.matrix(~ 0 + TS )
colnames(design) <- levels(TS)
dge <- DGEList(counts = STAR_Counts, group = groupLabels)

# filter out low expressed genes
cpm_cutoff=5; min.expr.counts=1; min.expr.num.samples=3

cutoff <- as.vector(cpm(cpm_cutoff,mean(dge$samples$lib.size) ) )
#keep <- rowSums(cpm(dge) > cutoff) >= min(as.numeric(table(groupLabels)))
keep <- rowSums(getCounts(dge) >= min.expr.counts)  >=  min.expr.num.samples
dge <- dge[keep, keep.lib.sizes = FALSE]
normMethod = "TMM"
dge <- calcNormFactors(dge, method=normMethod)
#getCounts(dge) #get expression matrix from 
v <- voom(dge, design, plot = T)
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(TGFBvsControl = (TGFB - Control), CX5461vsControl = (CX5461 - Control), CX5461vsTGFB = (CX5461 - TGFB), levels = design)
fitcon <- contrasts.fit(fit, cont.matrix)
fitcon <- eBayes(fitcon)

results1 <- topTable(fitcon, n = Inf, sort.by="P", coef="TGFBvsControl")
#outFile1 <- paste0("TGFBvsControl", "Limma_Voom_diff.txt")
#write.table(results1, file = outFile1, sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)


results2 <- topTable(fitcon, n = Inf, sort.by="P", coef="CX5461vsControl")
#outFile2 <- paste0("CX5461vsControl", "Limma_Voom_diff.txt")
#write.table(results2, file = outFile2, sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)


results3 <- topTable(fitcon, n = Inf, sort.by="P", coef="CX5461vsTGFB")
#outFile3 <- paste0("CX5461vsTGFB", "Limma_Voom_diff.txt")
#write.table(results3, file = outFile3, sep = "\t", col.names = TRUE, row.names = TRUE, quote = F)


### TGFbeta Vs Normal
results1$ENSG <- rownames(results1)
resdata1 <- left_join(results1,gencode.vM27.gtf.selected, by ="ENSG")
resdata1 <- resdata1 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
results1_adjP_0.05 <- resdata1 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & adj.P.Val <= FDR & abs(logFC) >= FoldChange)


# CX5461 Vs TGFbeta
results3$ENSG <- rownames(results3)
resdata3 <- left_join(results3,gencode.vM27.gtf.selected, by ="ENSG")
resdata3 <- resdata3 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
results3_adjP_0.05 <- resdata3 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & adj.P.Val <= FDR & abs(logFC) >= FoldChange)


###########################################################
###########################################################
###########################################################
############### DEG analysis by using DESeq2 ##############

library(DESeq2)
sampleTable <- data.frame(condition = factor(rep(c("Control", "TGFB", "TGFb_CX5461"), each = 3)))
rownames(sampleTable) <- colnames(STAR_Counts)
dds <- DESeq2::DESeqDataSetFromMatrix(round(RSEM_Counts), sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= 5
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
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR & abs(log2FoldChange) >= FoldChange)


### TGFbeta+CX5461 Vs Control #### This analysis don't required.
res2 <- DESeq2::results(dds1, contrast=c("condition","TGFb_CX5461","Control"))
DESeq2.res2 <- as.data.frame(res2)
DESeq2.res2$ENSG <- rownames(DESeq2.res2)
resdata2 <- left_join(DESeq2.res2,gencode.vM27.gtf.selected, by ="ENSG")
resdata2 <- resdata2 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
resdata2_adjP_0.05 <- resdata2 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR & abs(log2FoldChange) >= FoldChange)


### TGFbeta+CX5461 Vs. TGFbeta
res3 <- DESeq2::results(dds1, contrast=c("condition","TGFb_CX5461","TGFB"))
DESeq2.res3 <- as.data.frame(res3)
DESeq2.res3$ENSG <- rownames(DESeq2.res3)
resdata3 <- left_join(DESeq2.res3,gencode.vM27.gtf.selected, by ="ENSG")
resdata3 <- resdata3 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
resdata3_adjP_0.05 <- resdata3 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= FDR & abs(log2FoldChange) >= FoldChange)



###########################################################
###########################################################
###########################################################
### Differential expression analysis with FPKM ############

# filter out low expressed genes
count_matrix <- log2(RSEM_FPKM+1)
dataset1 <- count_matrix[apply(count_matrix,1,function(x) sum(x==0))<ncol(count_matrix)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
#keep <- rowSums(dataset1) >= 0.5 ## Remove very low readcount genes.
#dataset1 <- dataset1[keep,]

fit <- lmFit(dataset1, design)
cont.matrix <- makeContrasts(TGFBvsControl = (TGFB - Control), CX5461vsControl = (CX5461 - Control), CX5461vsTGFB = (CX5461 - TGFB), levels = design)
fitcon <- contrasts.fit(fit, cont.matrix)
fitcon <- eBayes(fitcon)

### TGFbeta Vs Normal
results1 <- topTable(fitcon, n = Inf, sort.by="P", coef="TGFBvsControl")
results1$ENSG <- rownames(results1)
resdata1 <- left_join(results1,gencode.vM27.gtf.selected, by ="ENSG")
resdata1 <- resdata1 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
results1_adjP_0.05 <- resdata1 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & adj.P.Val <= FDR & abs(logFC) >= FoldChange)

# CX5461 Vs TGFbeta
results3 <- topTable(fitcon, n = Inf, sort.by="P", coef="CX5461vsTGFB")
results3$ENSG <- rownames(results3)
resdata3 <- left_join(results3,gencode.vM27.gtf.selected, by ="ENSG")
resdata3 <- resdata3 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
results3_adjP_0.05 <- resdata3 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & adj.P.Val <= FDR & abs(logFC) >= FoldChange)

###########################################################
###########################################################
###########################################################
### Differential expression analysis with TPM #############

# filter out low expressed genes
count_matrix <- log2(RSEM_TPM+1)
dataset1 <- count_matrix[apply(count_matrix,1,function(x) sum(x==0))<ncol(count_matrix)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
#keep <- rowSums(dataset1) >= 0.5 ## Remove very low readcount genes.
#dataset1 <- dataset1[keep,]

fit <- lmFit(dataset1, design)
cont.matrix <- makeContrasts(TGFBvsControl = (TGFB - Control), CX5461vsControl = (CX5461 - Control), CX5461vsTGFB = (CX5461 - TGFB), levels = design)
fitcon <- contrasts.fit(fit, cont.matrix)
fitcon <- eBayes(fitcon)

### TGFbeta Vs Normal
results1 <- topTable(fitcon, n = Inf, sort.by="P", coef="TGFBvsControl")
results1$ENSG <- rownames(results1)
resdata1 <- left_join(results1,gencode.vM27.gtf.selected, by ="ENSG")
resdata1 <- resdata1 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
results1_adjP_0.05 <- resdata1 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & adj.P.Val <= FDR & abs(logFC) >= FoldChange)

# CX5461 Vs TGFbeta
results3 <- topTable(fitcon, n = Inf, sort.by="P", coef="CX5461vsTGFB")
results3$ENSG <- rownames(results3)
resdata3 <- left_join(results3,gencode.vM27.gtf.selected, by ="ENSG")
resdata3 <- resdata3 %>% 
  select(c("ENSG","Symbol","Chr","GeneType"), everything()) ## Reoder result file
results3_adjP_0.05 <- resdata3 %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & adj.P.Val <= FDR & abs(logFC) >= FoldChange)



#####################################################################
#####################################################################
#####################################################################
################## Filtering Step/remove probes #####################
Gene_filtering <- function(mat, Count_filter=0){
  not_all_na <- function(x) any(!is.na(x))
  mat <- mat %>%
    dplyr::select(where(not_all_na))
  mat <- subset(mat, apply(mat, 1, sum) != 0)
  #mat <- subset(mat, apply(mat, 2, sum) != 0)# not_all_na already did this
  mat <- mat[apply(mat,1,function(x) sum(x<=0))< ncol(mat)*0.75,]##Remove genes which have over 25% zero
  mat <- subset(mat, apply(mat, 1, sum) >= Count_filter)
  return(mat)}

#exp <- Gene_filtering(as.data.frame(STAR_Counts))
############ Function for the FPKM & FPKM UQ conversion #############
GDC_FPKM <- function(mat, gene_lengths, method = "FPKM") 
{
  genes <- intersect(rownames(mat), gene_lengths$Geneid)
  gene_lengths <- gene_lengths[genes,]
  mat <- mat[genes,]
  protein_coding <- gene_lengths[grep("protein_coding", gene_lengths$Class),]
  RC_g = mat
  RC_pc <- colSums(mat[rownames(protein_coding),])
  #mat.pc <- mat[rownames(protein_coding),]
  eff_leng <- gene_lengths$Length
  names(eff_leng) <- gene_lengths$Geneid
  #eff_leng <- t(eff_leng)
  if (method == "FPKM") {
    fpkm <- do.call(cbind, lapply(1:ncol(RC_g), function(i) {
      (((RC_g[, i]) * 1e+09) / (eff_leng *RC_pc[i]))
    }))
    colnames(fpkm) <- colnames(RC_g)
    rownames(fpkm) <- rownames(RC_g)
    return(fpkm)}
  
  if(method=="FPKM-UQ") {
    uqs <- apply(RC_g, 2, quantile, 0.75)
    fpkm <- do.call(cbind, lapply(1:ncol(RC_g), function(i) {
      (((RC_g[, i]) * 1e+09) / (eff_leng * uqs[i]))
    }))
    colnames(fpkm) <- colnames(RC_g)
    rownames(fpkm) <- rownames(RC_g)
    return(fpkm)}
}

######################################################
# https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e
# https://www.reneshbedre.com/blog/expression_units.html
# This RPKM, TPM, CPM is same as Python 
rpkm <- function(counts, lengths) {
  rate <- (counts / lengths) * 1e9
  rate / sum(counts) }

tpm <- function(counts, lengths) {
  rate <- (counts / lengths)*1e3
  (rate / sum(rate, na.rm = TRUE)) * 1e6}

cpm <- function(counts, lengths)
{  (counts * 1e6) /sum(counts)}


######################################################
####### This part is for getting common genes ########
common_Gene <- function(mat, geneLength){
  genes <- intersect(rownames(mat), geneLength$Geneid)
  gene_lengths <- geneLength[genes,]
  mat <- mat[genes,]
  return(list(Expression=mat, GeneLength=gene_lengths))}

######################################################
########### Filter and save normalized data ##########
ExpData <- Gene_filtering(mat = exp) ## Save after filtering
Common_List <- common_Gene(ExpData, gene_lengths)
rpkms <- apply(Common_List$Expression, 2, function(x) rpkm(x, Common_List$GeneLength$Length))
tpms <- apply(Common_List$Expression, 2, function(x) tpm(x, Common_List$GeneLength$Length))
cpms <- apply(Common_List$Expression, 2, function(x) cpm(x, Common_List$GeneLength$Length))
mat_FPKM <- as.data.frame(GDC_FPKM(mat = Common_List$Expression, gene_lengths = Common_List$GeneLength, method = "FPKM"))
mat_FPKM.UQ <- as.data.frame(GDC_FPKM(mat = Common_List$Expression, gene_lengths = Common_List$GeneLength, method = "FPKM-UQ"))

mat_FPKM.UQ.filter <- Gene_filtering(mat = mat_FPKM.UQ, Count_filter = 10)


###########################################################
###########################################################
###########################################################

#Hollingsworh Cell Cycle Progression (CCP) paper :: https://www.sciencedirect.com/science/article/pii/S1535610818305828?via%3Dihub
#Hypoxia Biomaker https:: //www.sciencedirect.com/science/article/pii/S0936655515003118
#Cell cycle Biomarker https :: //github.com/hbc/tinyatlas
## Column based scaling
#http://www.sthda.com/english/wiki/rna-sequencing-data-analysis-counting-normalization-and-differential-expression
## Why should normalize with DESeq2. Next we can normalize with SD and mean

dds <- DESeq2::DESeqDataSetFromMatrix(round(RSEM_Counts), sampleTable, ~condition)
dds <- estimateSizeFactors( dds )
normalized_counts <- counts(dds, normalized=TRUE)

logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
logcounts <- Gene_filtering(mat = as.data.frame(logcounts), Count_filter = 1)

Selected_Gene <- function(mat, geneLength){
  genes <- intersect(rownames(mat), geneLength$Geneid)
  gene_lengths <- geneLength[genes,]
  gene_lengths$ENSG_ID <- gsub("\\..*","", gene_lengths$Geneid)
  mat <- mat[genes,]
  return(list(Expression=mat, GeneLength=gene_lengths))}



#scale(log2(RSEM_FPKM+1), center = TRUE, scale = TRUE) ##https://medium.com/swlh/data-normalisation-with-r-6ef1d1947970
#exp.scaled <- (scale(log2(t(RSEM_FPKM+1)), center = TRUE, scale = TRUE))

Selected_Gene_Exp <- Selected_Gene(mat = logcounts, geneLength = gene_lengths)

# https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html
mouse.cell.cycle <- read.csv("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
# www.pnas.org/cgi/doi/10.1073/pnas.1818210116 #Stemness marker
#https://satijalab.org/seurat/reference/cc.genes.updated.2019.html#source
#Seurat::cc.genes.updated.2019
median(colMeans(logcounts)) ## omit genes below the median of average expression across samples Biranda et al. PNAS Stemness paper

cc.S <- mouse.cell.cycle[mouse.cell.cycle$phase=="S",]
cc.G2 <- mouse.cell.cycle[mouse.cell.cycle$phase=="G2/M",]
CCP_gene <- stringr::str_to_sentence(c("PLK1", "RRM2", "CDCA8", "CDC20", "CDCA3", "FOXM1", "BIRC5", "PBK", "CDKN3", "CENPM", "MCM10", "TK1", "DTL", "ASF1B", "NUSAP1", "PRC1", "BUB1B", "DLGAP5", "TOP2A", "ASPM", "CENPF", "KIF20A", "PTTG1", "RAD54L", "ORC6", "CDK1", "CEP55", "KIF11", "RAD51", "KIAA0101", "SKA1"))
Housekeeping_gene <- stringr::str_to_sentence(c("RPL38", "UBA52", "PSMC1", "RPL4", "RPL37", "RPS29", "SLC25A3", "CLTC", "TXNL1", "PSMA1", "RPL8", "MMADHC", "RPL13A", "PPP2CA", "MRFAP1"))
Hypoxia_Gene <- stringr::str_to_sentence(c("BNIP3", "F3", "LOX", "TNF", "TH", "SLC2A1", "PGK1", "NDRG1", "GAL", "BNIP3L", "ANG", "P4HA1", "ADM", "AK3", "PDK1", "ERO1L", "ALDOC", "PLOD2", "P4HA2", "MXI1", "DDIT4", "ANGPTL4"))
Mesenchymal_Gene <- stringr::str_to_sentence(c("VIM", "CDH2", "FOXC2", "SNAI1", "SNAI2", "TWIST1", "FN1", "ITGB6", "MMP2", "MMP3", "MMP9", "SOX10", "GCS"))
Epithelial_Gene <- stringr::str_to_sentence(c("CDH1", "DSP", "OCLN"))
STEM_Gene <- stringr::str_to_sentence(c("DNMT3B", "PFAS", "XRCC5", "HAUS6", "TET1", "IGF2BP1", "PLAA", "TEX10", "MSH6", "DLGAP5", "SKIV2L2", "SOHLH2", "RRAS2", "PAICS", "CPSF3", "LIN28B", "IPO5", "BMPR1A", "ZNF788", "ASCC3", "FANCB", "HMGA2", "TRIM24", "ORC1", "HDAC2", "HESX1", "INHBE", "MIS18A", "DCUN1D5", "MRPL3", "CENPH", "MYCN", "HAUS1", "GDF3", "TBCE", "RIOK2", "BCKDHB", "RAD1", "NREP", "ADH5", "PLRG1", "ROR1", "RAB3B", "DIAPH3", "GNL2", "FGF2", "NMNAT2", "KIF20A", "CENPI", "DDX1", "XXYLT1", "GPR176", "BBS9", "C14orf166", "BOD1", "CDC123", "SNRPD3", "FAM118B", "DPH3", "EIF2B3", "RPF2", "APLP1", "DACT1", "PDHB", "C14orf119", "DTD1", "SAMM50", "CCL26", "MED20", "UTP6", "RARS2", "ARMCX2", "RARS", "MTHFD2", "DHX15", "HTR7", "MTHFD1L", "ARMC9", "XPOT", "IARS", "HDX", "ACTRT3", "ERCC2", "TBC1D16", "GARS", "KIF7", "UBE2K", "SLC25A3", "ICMT", "UGGT2", "ATP11C", "SLC24A1", "EIF2AK4", "GPX8", "ALX1", "OSTC", "TRPC4", "HAS2", "FZD2", "TRNT1", "MMADHC", "SNX8", "CDH6", "HAT1", "SEC11A", "DIMT1", "TM2D2", "FST", "GBE1"))




S_Score <- Selected_Gene_Exp$Expression%>% filter(Selected_Gene_Exp$GeneLength$ENSG_ID %in% cc.S$geneID) %>% select(everything()) %>% colSums()
G2M_Score <- Selected_Gene_Exp$Expression%>% filter(Selected_Gene_Exp$GeneLength$ENSG_ID %in% cc.G2$geneID) %>% select(everything()) %>% colSums()
CCP <- Selected_Gene_Exp$Expression %>% filter(Selected_Gene_Exp$GeneLength$GeneSymbol %in% CCP_gene) %>% select(everything()) %>% colSums()
Housekeeping_Score <- Selected_Gene_Exp$Expression %>% filter(Selected_Gene_Exp$GeneLength$GeneSymbol %in% Housekeeping_gene) %>% select(everything()) %>% colSums()
Hypoxia_Score <- Selected_Gene_Exp$Expression %>% filter(Selected_Gene_Exp$GeneLength$GeneSymbol %in% Hypoxia_Gene) %>% select(everything()) %>% colSums()
Mesenchymal_Score <- Selected_Gene_Exp$Expression %>% filter(Selected_Gene_Exp$GeneLength$GeneSymbol %in% Mesenchymal_Gene) %>% select(everything()) %>% colSums()
Epithelial_Score <- Selected_Gene_Exp$Expression %>% filter(Selected_Gene_Exp$GeneLength$GeneSymbol %in% Epithelial_Gene) %>% select(everything()) %>% colSums()

Score_Total <- as.data.frame(t(rbind(S_Score, G2M_Score, CCP, Housekeeping_Score, Hypoxia_Score, Mesenchymal_Score, Epithelial_Score)))

Score_Total <- Score_Total %>%
  select(everything()) %>%
  mutate(EMT_Score=scales::rescale(Mesenchymal_Score-Epithelial_Score, to = c(0,1)), 
         CCP_Score=scales::rescale(CCP-Housekeeping_Score, to=c(0,1)), 
         S_Score=scales::rescale(S_Score, to = c(0,1)),
         G2M_Score=scales::rescale(G2M_Score, to = c(0,1)),
         Hypoxia_Score=scales::rescale(Hypoxia_Score, to = c(0,1)),
         Sample=rep(c("Control", "TGFB", "TGFb_CX5461"), each = 3)) %>%
  select(c("Sample","EMT_Score","CCP_Score","S_Score","G2M_Score", "Hypoxia_Score"), everything())

#Score_Total %>% 
#  mutate(Scaled_EMT=scales::rescale(Score_Total$EMT_Score, to = c(-1,1)))
write.csv(Score_Total, file = "EMT_Score_File.csv")

