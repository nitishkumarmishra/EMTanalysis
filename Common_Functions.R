# ###################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2021
# Email:  nitishimtech@gmail.com
# Date: 2021-09-22
# Script Name: Common_Functions.R
#####################################################
# #############  Script Description: ################
# R code for basic functions which we need on regular basis
# Reading STAR Conts file and RSEM files and make matrix for DEG analysis
# Filtering of genes which are lowely expressed
# Convert Read counts in FPKM and FPKM-UQ based on gene length
###########################################################
###########################################################
load("C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/gene_length.rds") ##gene length for FPKM/UQ
library(dplyr); library('rtracklayer')
###########################################################
###########################################################

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


#STAR_Counts <- STAR_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/STAR-COUNTS/", Stranded = "reverse")


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

#RSEM_Counts <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "expected_count")
#RSEM_TPM <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "TPM")
#RSEM_FPKM <- RSEM_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS/Counts/", count_type = "FPKM")


############## HTSeq counts matrix ############
HTSeq_Matrix <- function(directory="", Stranded="reverse", counts="CDS")
{
  if (counts=="genes") {
    basedir <- directory
    setwd(basedir)
    files <- dir(pattern=".htseq.gene.counts", path = basedir)
    counts <- c()
    for( i in seq_along(files) ){
      x <- read.table(file=files[i], sep="\t", header=F, as.is=T, row.names = 1) # First four lines are statistics
      colnames(x) <- c("reverse")
      x <- head(x, -5)
      counts <- cbind(counts, x[,grep(Stranded, colnames(x))]) ## In STAR column 1) ID 2) Unstranded 3) Farward Strand 4) Reverse Strand
    }
    rownames(counts) <- rownames(x) # set the row names
    colnames(counts) <- sub(".htseq.gene.counts","",files) # set the column names based on input file names, with pattern removed
    return(counts)
  }
  if (counts=="CDS") {
    basedir <- directory
    setwd(basedir)
    files <- dir(pattern=".htseq.CDS.counts", path = basedir)
    counts <- c()
    for( i in seq_along(files) ){
      x <- read.table(file=files[i], sep="\t", header=F, as.is=T, row.names = 1) # First four lines are statistics
      colnames(x) <- c("reverse")
      x <- head(x, -5)
      counts <- cbind(counts, x[,grep(Stranded, colnames(x))]) ## In STAR column 1) ID 2) Unstranded 3) Farward Strand 4) Reverse Strand
    }
    rownames(counts) <- rownames(x) # set the row names
    colnames(counts) <- sub(".htseq.CDS.counts","",files) # set the column names based on input file names, with pattern removed
    return(counts)
  }
  
}

#HTSeq_Counts_genes <- HTSeq_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/HTSeq-Counts/Gene/", Stranded = "reverse", counts = "genes")
#HTSeq_Counts_CDS <- HTSeq_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/HTSeq-Counts/CDS/", Stranded = "reverse", counts = "CDS")


############## FeatureCounts matrix ############
FeatureCounts_Matrix <- function(directory="", Stranded="reverse", counts="CDS")
{
  if (counts=="genes") {
    basedir <- directory
    setwd(basedir)
    files <- dir(pattern="_featureCounts.gene.out", path = basedir)
    counts <- c()
    for( i in seq_along(files) ){
      x <- read.table(file=files[i], sep="\t", header=F, as.is=T, row.names = 1, skip = 2) # First four lines are statistics
      colnames(x) <- c("V2", "V3", "V4", "V5", "V6", "reverse")
      counts <- cbind(counts, x[,grep(Stranded, colnames(x))]) ## In STAR column 1) ID 2) Unstranded 3) Farward Strand 4) Reverse Strand
    }
    rownames(counts) <- rownames(x) # set the row names
    colnames(counts) <- sub("_featureCounts.gene.out","",files) # set the column names based on input file names, with pattern removed
    return(counts)
  }
  if (counts=="CDS") {
    basedir <- directory
    setwd(basedir)
    files <- dir(pattern="_featureCounts.CDS.out", path = basedir)
    counts <- c()
    for( i in seq_along(files) ){
      x <- read.table(file=files[i], sep="\t", header=F, as.is=T, row.names = 1, skip = 2) # First four lines are statistics
      colnames(x) <- c("V2", "V3", "V4", "V5", "V6", "reverse")
      counts <- cbind(counts, x[,grep(Stranded, colnames(x))]) ## In STAR column 1) ID 2) Unstranded 3) Farward Strand 4) Reverse Strand
    }
    rownames(counts) <- rownames(x) # set the row names
    colnames(counts) <- sub("_featureCounts.CDS.out","",files) # set the column names based on input file names, with pattern removed
    return(counts)
  }
  
}
#FeatureCounts_genes <- FeatureCounts_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/FeatureCounts/Gene/", Stranded = "reverse", counts = "genes")
#FeatureCounts_CDS <- FeatureCounts_Matrix(directory = "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/Translation/FeatureCounts/CDS/", Stranded = "reverse", counts = "CDS")

#####################################################################
#####################################################################
################## Filtering Step/remove probes #####################

########### Filtering Step/remove probes #############
Gene_filtering <- function(mat, Count_filter=0, min.prop=0.25){
  not_all_na <- function(x) any(!is.na(x))
  if(!isTRUE(is.data.frame(mat))) # If not dataFrame convert it.
    mat <- as.data.frame(mat)
  mat <- mat %>%
    dplyr::select(where(not_all_na))
  mat <- subset(mat, apply(mat, 1, sum) != 0)
  #mat <- subset(mat, apply(mat, 2, sum) != 0)# not_all_na already did this
  #mat <- mat[apply(mat,1,function(x) sum(x<=0))< ncol(mat)*0.75,]##Remove genes which have over 25% zero
  mat <- mat[apply(mat,1,function(x) sum(x<=0))< ncol(mat)*(1-min.prop),]##Remove genes which have over 25% zero. We can change min.prop accordingly.
  mat <- subset(mat, apply(mat, 1, sum) >= Count_filter)
  return(mat)}

#exp <- Gene_filtering(as.data.frame(STAR_Counts), Count_filter=5, min.prop=0.25)

#####################################################################
#####################################################################
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

# STAR_FPKM <- GDC_FPKM(mat = STAR_Counts, gene_lengths = gene_lengths, method = "FPKM")
# STAR_FPKM_UQ <- GDC_FPKM(mat = STAR_Counts, gene_lengths = gene_lengths, method = "FPKM-UQ")

#####################################################################
#####################################################################
################## Normalized read count matrix #####################
# https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e
# https://www.reneshbedre.com/blog/expression_units.html
# This RPKM, TPM, CPM is same as Python 
rpkm <- function(counts, lengths) {
  rate <- (counts / lengths) * 1e9
  rate / sum(counts, na.rm = TRUE) }

tpm <- function(counts, lengths) {
  rate <- (counts / lengths)*1e3
  (rate / sum(rate, na.rm = TRUE)) * 1e6}

cpm <- function(counts, lengths)
{  (counts * 1e6) /sum(counts, na.rm = TRUE)}

#####################################################################
#####################################################################
################## verb function for name print #####################
verb <- function(...) cat(sprintf(...), sep='', file=stdout())
#####################################################################

#####################################################################
#####################################################################
######################## plotly 3D PCA Plot #########################
plotPCA3D <- function (mat, filter=T, discard=TRUE, sample=NA, GeneSelect=TRUE, ntop = 500, file="3D-PCA-Plot"){
  #mat <- as.matrix(mat)
  not_all_na <- function(x) any(!is.na(x))
  mat <- mat %>%
    select(where(not_all_na))
  mat <- subset(mat, apply(mat, 1, sum) != 0)
  #mat <- subset(mat, apply(mat, 2, sum) != 0)# not_all_na already did this
  if(filter==T){
    mat <- mat[apply(mat,1,function(x) sum(x==0))< ncol(mat)*0.75,]##Remove genes which have over 25% zero
    mat <- subset(mat, apply(mat, 1, sum) >= 10)  }
  
  if((discard==TRUE) && (!is.na(sample)))
  {
    mat <- mat[, !names(mat)%in%sample]
    #sampletable <- sampletable[colnames(mat),]
    sampletable <- sampletable[!rownames(sampletable)%in%sample,]
  }
  rv <- matrixStats::rowVars(as.matrix(mat))
  if(GeneSelect==TRUE && (!is.na(ntop)))
  { 
    select <- order(rv, decreasing = TRUE)
    select <- select[1:ntop]
    mat <- mat[select,] 
  }
  #select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(mat), scale. = TRUE)
  percentVar <- round(pca$sdev^2/sum(pca$sdev^2),3)*100
  group <- sampletable$Class
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
                  group = group, name = colnames(mat))
  message("Generating plotly plot")
  fig <- plotly::plot_ly(data = d, x = ~PC1, y = ~PC2, z = ~PC3, color = group, text = rownames(d))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(
    title = "Layout options in a 3d scatter plot",
    scene = list(
      xaxis = list(title = paste0("Comp 1: ", percentVar[1], "%", sep = "")),
      yaxis = list(title =  paste0("Comp 2: ", percentVar[2], "%", sep = "")),
      zaxis = list(title = paste0("Comp 3: ", percentVar[3], "%", sep = ""))
    ))
  #return(fig)
  htmlwidgets::saveWidget(as_widget(fig), paste0(file, ".html"))
}

#plotPCA3D(mat = mat, file = "3D-PCA-Plot")
#plotPCA3D(mat = mat, GeneSelect = TRUE, ntop = 1000,file = "3D-PCA-Plot-select10000")

  