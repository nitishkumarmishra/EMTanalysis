## ---------------------------------------------------------
## Script name: Volcano_mTORC_Proteomics.R
## Purpose of script: make volcano plot for proteomics data
## Author: Nitish Kumar Mishra
## Date Created: 2023-02-23
## Copyright (c) Nitish Mishra, 2023
## Email: nitish.mishra@stjude.org
## ---------------------------------------------------------
## Notes:
## Gene list for volcano plot for Jake Bachelder data presentation 
## Non-significant with "grey", significant with log2/FDR with black and mTORC2 components in purple

## ---------------------------------------------------------


################################################################################
####### ************* load library and set working path ************** ########
################################################################################

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(EnhancedVolcano)))

setwd("/Users/nmishra/Dropbox/Mac/Desktop")


################################################################################
####### ************* Make color file for selected genes ************** ########
################################################################################
mTORC2 <- c("Mapkap1", "Rictor", "mLST8", "mTOR", "Prr5", "Prr5l")
colo <- rep("purple", length(mTORC2))
mTORC2.col <- cbind(mTORC2, colo) 

Ribosomal_protein <- c("Rack1", "Rps8", "Rps26", "Rps16", "Rps3a", "Rps15", "Rps2", "Rps11", "Rps25", "Rps15a", "Rps12", "Rps23", "Rps6", "Rps3", "Rps7", "Rps18", "Rps19", "Rps4x", "Rps5", "Rps13", "Rps9", "Rps21", "Rps28", "Rpsa", "Rps10", "Rps14", "Rps17", "Rps20", "Rpl10", "Rpl18a", "Rpl23a", "Rpl30", "Rpl3", "Rpl14", "Rpl34", "Rpl24", "Rpl10a", "Rpl6", "Rpl9", "Rpl7", "Rpl8", "Rpl17", "Rpl5", "Rpl15", "Rpl18", "Rpl35a", "Rpl4", "Rplp0", "Rpl38", "Rpl26", "Rpl12", "Rpl27a", "Rplp2", "Rpl23", "Rpl28", "Rpl27", "Rpl22l1", "Rpl21", "Rpl36a", "Rpl31", "Rpl7a", "Rpl11", "Rpl37", "Rpl13", "Rpl22", "Rpl32", "Rpl29", "Rpl35", "Rpl19")
colo <- rep("orange", length(Ribosomal_protein))
Ribosomal_protein.col <- cbind(Ribosomal_protein, colo)

TI_factor <- c("eIF4g1", "eIF4e", "eIF4b", "eIF3h", "eIF2s2", "eIF3k", "eIF5a", "eIF2s1", "eIF2s3x", "eIF3c", "eIF3i", "eIF3l", "eIF5b", "eIF3a", "eIF3b", "eIF4g2")
colo <- rep("green", length(TI_factor))
TI_factor.col <- cbind(TI_factor, colo)

TE_factor <- c("eEF1b", "eEF1g", "eEF1d", "eEF1e1", "eEF2", "eEF1a1")
colo <- rep("magenta", length(TE_factor))
TE_factor.col <- cbind(TE_factor, colo)

selected.genes <- rbind(Ribosomal_protein.col, TE_factor.col, TI_factor.col)


################################################################################
########## *************** Read proteomics excel data **************** #########
################################################################################
data <- readxl::read_excel("Blancgrp_122222_sin1_IP_TMT_Project_Report.V1.0.1_total.xlsx", sheet = "Data")

data1 <- data[data$GN!="NA",]

################################################################################
###### *************** Ctrl Vs. WT data and volcano plot **************** ######
################################################################################
data2 <- 
  data1 %>%
  dplyr::distinct(GN, .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "GN") %>%
  dplyr::mutate(`Log2Fold(WT_Vs_Ctrl)` = (-1 * `Log2Fold(WT_Vs_Ctrl)`)) # Make Ctrl Vs. WT by multiplying -1
  

## **** Limit for the X/Y axis *** ##
xlim=max(abs(data2$`Log2Fold(WT_Vs_Ctrl)`))+1
ylim=max(-log10(data2$FDR_WT_vs_cntrl))+1

## ***** Make genes name font in italics **** ##
lab_italics <- paste0("italic('", rownames(data2), "')")
selectLab_italics = paste0("italic('",mTORC2.col[,1]  ,"')")

## ******** EnhancedVolcano plot ********** ##
p <- EnhancedVolcano(data2,
                lab = lab_italics,
                x = 'Log2Fold(WT_Vs_Ctrl)',
                y = 'FDR_WT_vs_cntrl',
                selectLab = selectLab_italics,
                xlab = bquote(~Log[2]~ '(Fold Change)'),
                ylab  = bquote(~Log[10]~'(FDR)'),
                xlim = c(-xlim, xlim),
                ylim = c(0, ylim),
                title = "Ctrl Vs. WT",
                subtitle = "mTOR bound ribosomes proteome",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 4.0,
                labCol = 'purple',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = TRUE,
                col = c('grey', 'grey', 'grey', 'black'),
                colAlpha = 1,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'purple',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                caption = "") +
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))
ggsave(p, file="Ctrl_Vs_WT.pdf", height = 8, width = 8, dpi = 500)

################################################################################
##### *************** TGFb Vs. Ctrl data and volcano plot **************** #####
################################################################################
data2 <- 
  data1 %>%
  dplyr::distinct(GN, .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "GN") %>%
  dplyr::mutate(`Log2Fold(Ctrl_Vs_TGFb)` = (-1 * `Log2Fold(Ctrl_Vs_TGFb)`)) # Make Ctrl Vs. WT by multiplying -1

# create custom key-value pairs for 'Ribosomal', 'Trans. Init.', 'Trans. Elong.' expression by fold-change
# this can be achieved with nested ifelse statements
# Here I will use 

pCutoff <- 0.05
FCcutoff <- 0.5

keyvals <- ifelse(rownames(data2) %in% Ribosomal_protein.col[,1], 'blue',
                  ifelse(rownames(data2) %in% stringr::str_to_title(TE_factor.col[,1]), 'green', 
                         ifelse(rownames(data2) %in% stringr::str_to_title(TI_factor.col[,1]), 'firebrick',
                                ifelse(abs(data2$`Log2Fold(Ctrl_Vs_TGFb)`) >= FCcutoff & data2$FDR_cntrl_vs_TGFb <= pCutoff, 'black','grey'))))
#keyvals[is.na(keyvals)] <- 'grey'


names(keyvals)[keyvals == 'blue'] <- 'Riboproteins'
names(keyvals)[keyvals == 'green'] <- 'TE factors'
names(keyvals)[keyvals == 'firebrick'] <- 'TI factors'
names(keyvals)[keyvals == 'black'] <- 'DE proteins'
names(keyvals)[keyvals == 'grey'] <- 'Not Signi.'           


## **** Limit for the X/Y axis *** ##
xlim=max(abs(data2$`Log2Fold(Ctrl_Vs_TGFb)`))+1
ylim=max(-log10(data2$FDR_cntrl_vs_TGFb))+1

## ***** Make genes name font in italics **** ##
#lab_italics <- paste0("italic('", rownames(data2), "')")
#selectLab_italics = paste0("italic('",mTORC2.col[,1]  ,"')")
selectLab_italics <- NULL


## ******** EnhancedVolcano plot ********** ##
p1 <- EnhancedVolcano(data2,
                lab = rownames(data2),
                x = 'Log2Fold(Ctrl_Vs_TGFb)',
                y = 'FDR_cntrl_vs_TGFb',
                selectLab = c("NA"),
                xlab = bquote(~Log[2]~ '(Fold Change)'),
                ylab  = bquote(~Log[10]~'(FDR)'),
                xlim = c(-xlim, xlim),
                ylim = c(0, ylim),
                title = "TGFb Vs. Ctrl",
                subtitle = "mTOR bound ribosomes proteome",
                #legendPosition = 'right',
                colCustom = keyvals,
                colAlpha = 1,
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = c(ifelse(rownames(data2) %in% stringr::str_to_title(selected.genes[,1]), 2.0, 1)),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                caption = "") +
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5))

ggsave(p1, file="TGFb Vs Ctrl.pdf", height = 8, width = 10, dpi = 500)
