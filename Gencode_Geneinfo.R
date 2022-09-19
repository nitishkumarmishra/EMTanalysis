# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2022
# Email:  nitishimtech@gmail.com
# Date: 25 January, 2022
# Script Name: Gencode_Geneinfo.R
#####################################################################
# ##################  Script Description: ###########################
# R code to generate gene information file from GTF file
# This gene information files is similar to TCGA
# We will use this file for FPKM and FPKM-UQ calculation
##################### Notes: ########################################
# Length of genes by using all exon as in TCGA information file for FPKM/FPKM-UQ count calculation
# https://www.biostars.org/p/185665/#9507426;  http://rseqc.sourceforge.net/#fpkm-uq-py
# Python code GTFtools will also work http://www.genemine.org/gtftools.php
#####################################################################
################## SET WORKING DIRECTORY ############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/jupyter codes/"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################################
library(GenomicFeatures)
library(dplyr)

################ GRCm39 MM39 GENCODE annotation #####################

txdb <- makeTxDbFromGFF("gencode.vM28.annotation.gtf.gz",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then :: https://www.biostars.org/p/83901/
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
gene.exon.counts <- lengths(exons.list.per.gene)
gene.info <- cbind(exonic.gene.sizes, gene.exon.counts)
colnames(gene.info) <- c("exon_length", "exon_number")


gtf <- rtracklayer::import("gencode.vM28.annotation.gtf.gz")
gencode.vM28.gtf <- as.data.frame(gtf)
gencode.vM28.gtf.selected <- gencode.vM28.gtf %>%
  filter(type=="gene") %>%
  rename(Chr= seqnames, Type= type, Gene = gene_id, GeneType=gene_type, Symbol= gene_name, Havana=havana_gene, Length=width, Start = start, End=end, Strand=strand) %>%
  select(Gene, Symbol,  Chr, Start, End, Strand,  GeneType,  Havana, Length) #%>%
#filter(!grepl("chrM",Chr))

gene.info <- merge(gencode.vM28.gtf.selected, gene.info, by.x="Gene", by.y= "row.names")
write.table(gene.info, file = "geneinfo.gencode.vM28.txt", na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


################################################################
# This part is for Human Gencode
################### Option 1 ###################################
################################################################

txdb <- makeTxDbFromGFF("gencode.v22.annotation.gtf.gz", format="gtf")
exonic <- exonsBy(txdb, by="gene")
red.exonic <- reduce(exonic)
## Length of genes by using all exon as in TCGA information file for FPKM/FPKM-UQ count calculation
# https://www.biostars.org/p/185665/#9507426;  http://rseqc.sourceforge.net/#fpkm-uq-py
exon.lengths <- vapply(width(red.exonic), sum, numeric(1))
### Number of exons for each genes for TCGA feature files
exon.counts <- sapply(exonic,length)
#exon.counts <- lengths(exonic) ## This is alternate to sapply command to counts number of exons
gene.info <- cbind(exon.lengths, exon.counts)

################### Option 2 ###################################
################################################################

txdb <- makeTxDbFromGFF("gencode.v39.annotation.gtf",format="gtf")
# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
gene.exon.counts <- lengths(exons.list.per.gene)
gene.info <- cbind(exonic.gene.sizes, gene.exon.counts)
colnames(gene.info) <- c("exon_length", "exon_number")


gtf <- rtracklayer::import("gencode.v39.annotation.gtf")
gencode.v39.gtf <- as.data.frame(gtf)
gencode.v39.gtf.selected <- gencode.v39.gtf %>%
  filter(type=="gene") %>%
  rename(Chr= seqnames, Type= type, Gene = gene_id, GeneType=gene_type, Symbol= gene_name, Havana=havana_gene, Length=width, Start = start, End=end, Strand=strand) %>%
  select(Gene, Symbol,  Chr, Start, End, Strand,  GeneType,  Havana, Length) #%>%
#filter(!grepl("chrM",Chr))

gene.info <- merge(gencode.v39.gtf.selected, gene.info, by.x="Gene", by.y= "row.names")
write.table(gene.info, file = "geneinfo.gencode.V39.txt", na = "NA", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


