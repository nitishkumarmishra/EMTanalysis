# ###################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2021
# Email:  nitishimtech@gmail.com
# Date: 2021-09-22
# Script Name: HK_EXIT_STRATEGY_ReRun.R
#####################################################
# #############  Script Description: ################
#
############### Notes: ##############################
##Hyunsoo Kim Analyses
# This code is just for checking discripencies in results
### Set path of "HK EXIT STRATEGY" folder, I copied it locally on my Desktop.
#####################################################
########## SET WORKING DIRECTORY ####################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/HK EXIT STRATEGY/jupyter codes"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################



(suppressWarnings(source("./jupyter_common.R")))

load('./rdata/jupyter_common.rdata')
load('./rdata/riboprof_170224.rdna_rn18s_161021.rdna_rn18s_limma-voom.mrna.rdata')



##############################################################
##############################################################
# common parameters

f_display_sym <- FALSE
seed_geom_text_repel <- 40
#biotype_exclude <- NULL
biotype_exclude <- "pseudogene"  # current choice of fig1
#biotype_exclude <- "protein_coding|pseudogene"
#biotype_include <- NULL
biotype_include <- "protein_coding"  # current choice of fig1
fname_appendix <- ""
if (!is.null(biotype_include)) {
  fname_appendix <- sprintf("_%s", gsub("\\|", "_", biotype_include))
} 
if (!is.null(biotype_exclude)) {
  fname_appendix <- sprintf("%s_wo_%s", fname_appendix, gsub("\\|", "_", biotype_exclude))
} 
fname_appendix


if (!is.null(biotype_exclude)) {
  f <- grepl(biotype_exclude, df_all$biotype)
  df_exclude <- df_all[f, ]
  show(head(df_exclude))
  show(dim(df_exclude))
  excludes <- rownames(df_exclude)
  
  sym_mrna_only <- setdiff(sym_mrna_only, excludes)
  sym_mrna_up_ribo_up <- setdiff(sym_mrna_up_ribo_up, excludes)
  sym_mrna_dn_ribo_dn <- setdiff(sym_mrna_dn_ribo_dn, excludes)
  sym_mrna_no_ribo_up <- setdiff(sym_mrna_no_ribo_up, excludes)
  sym_mrna_no_ribo_dn <- setdiff(sym_mrna_no_ribo_dn, excludes)
  sym_mrna_up_ribo_dn <- setdiff(sym_mrna_up_ribo_dn, excludes)
  sym_mrna_dn_ribo_up <- setdiff(sym_mrna_dn_ribo_up, excludes)    
}



if (!is.null(biotype_include)) {
  f <- grepl(biotype_include, df_all$biotype)
  df_include <- df_all[f, ]
  show(head(df_include))
  show(dim(df_include))
  include_genes <- rownames(df_include)
  
  sym_mrna_only <- intersect(sym_mrna_only, include_genes)
  sym_mrna_up_ribo_up <- intersect(sym_mrna_up_ribo_up, include_genes)
  sym_mrna_dn_ribo_dn <- intersect(sym_mrna_dn_ribo_dn, include_genes)
  sym_mrna_no_ribo_up <- intersect(sym_mrna_no_ribo_up, include_genes)
  sym_mrna_no_ribo_dn <- intersect(sym_mrna_no_ribo_dn, include_genes)
  sym_mrna_up_ribo_dn <- intersect(sym_mrna_up_ribo_dn, include_genes)
  sym_mrna_dn_ribo_up <- intersect(sym_mrna_dn_ribo_up, include_genes)    
}



##############################################################
##############################################################
#rpf_only_up_dn

df_fig <- df_all
sym <- rownames(df_fig)
df_fig$fig.type <- NA

verb("sym_mrna_only: %d\n", length(sym_mrna_only))
f <- sym %in% sym_mrna_only
str_mrna_only <- sprintf("RNA Only (n=%d)", length(sym_mrna_only))
df_fig[f, "fig.type"] <- str_mrna_only

verb("sym_mrna_up_ribo_up: %d\n", length(sym_mrna_up_ribo_up))
verb("sym_mrna_dn_ribo_dn: %d\n", length(sym_mrna_dn_ribo_dn))
f <- sym %in% sym_mrna_up_ribo_up | sym %in% sym_mrna_dn_ribo_dn
str_mrna_rpf <- sprintf("RNA and RPF (n=%d)", length(sym_mrna_up_ribo_up)+length(sym_mrna_dn_ribo_dn))
df_fig[f, "fig.type"] <- str_mrna_rpf

verb("sym_mrna_no_ribo_up: %d\n", length(sym_mrna_no_ribo_up))
f <- sym %in% sym_mrna_no_ribo_up
str_rpf_only_up <- sprintf("RPF Only - Up (n=%d)", length(sym_mrna_no_ribo_up))
df_fig[f, "fig.type"] <- str_rpf_only_up

verb("sym_mrna_no_ribo_dn: %d\n", length(sym_mrna_no_ribo_dn))
f <- sym %in% sym_mrna_no_ribo_dn
str_rpf_only_dn <- sprintf("RPF Only - Down (n=%d)", length(sym_mrna_no_ribo_dn))
df_fig[f, "fig.type"] <- str_rpf_only_dn
# verb('sym_no_de_both: %d\n',length(sym_no_de_both)); f <- sym %in%
# sym_no_de_both; df_fig[f, 'fig.type'] <- 'No DE'

f <- !is.na(df_fig$fig.type) & !is.na(df_fig$log2FCunt48VStgfb48.transcription) & 
  !is.na(df_fig$log2FCunt48VStgfb48.translation)
df_fig <- df_fig[f, ]
df_fig$fig.type <- factor(df_fig$fig.type, level = c(str_mrna_only, str_mrna_rpf, str_rpf_only_up, str_rpf_only_dn))
df_fig <- df_fig[order(df_fig$fig.type), ]

nv_color <- c('RNA Only'='#aaaaaa',
              'RNA and RPF'='#7aaa3d',
              'RPF Only - Up'='#d78c5e',
              'RPF Only - Down'='#fcf050')
names(nv_color) <- levels(df_fig$fig.type)



nticks=13; xmax=12; ymax=12; gap_tick_label=4
gg <- ggplot(data=df_fig,
             aes(x=log2FCunt48VStgfb48.transcription, 
                 y=log2FCunt48VStgfb48.translation, colour=fig.type)) + 
  theme_geometry(ticks=nticks, xlim=xmax, ylim=ymax, linesize=0.5, 
                 xlab=expression('RNA-seq log'[2]*'(fold change)'),
                 ylab=expression('RPF-seq log'[2]*'(fold change)'),
                 labsize=3.5, labgap=1.5, epsilon=max(xmax,ymax)/50, gap_tick_label=gap_tick_label) +
  geom_point(alpha=1, size=1) +
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=9),
        legend.background = element_rect(color = NA),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0,"cm"),
        #legend.spacing.y = unit(0.5, "cm"),          
        legend.position = c(0.02, 0.92), legend.justification = c(0, 1) ) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ggtitle("") +
  scale_colour_manual(values = nv_color)

if (f_display_sym) {
  
  df_fig$sig <- 'no'; 
  df_fig$sym <- rownames(df_fig)  
  f_sym <- grepl('Malat', df_fig$sym)
  
  th_log2fc <- log2(1.2)
  f_tgfb48_rnaseq_up_rpf_up <- !is.na(df_fig$fig.type) & df_fig$log2FCunt48VStgfb48.transcription > th_log2fc & df_fig$log2FCunt48VStgfb48.translation > th_log2fc
  f_tgfb48_rnaseq_dn_rpf_up <- !is.na(df_fig$fig.type) & df_fig$log2FCunt48VStgfb48.transcription < -th_log2fc & df_fig$log2FCunt48VStgfb48.translation > th_log2fc
  df_fig$sig[f_tgfb48_rnaseq_up_rpf_up & f_sym] <- 'up_up';
  df_fig$sig[f_tgfb48_rnaseq_dn_rpf_up & f_sym] <- 'dn_up';
  #f <- grepl("Gm|Rik", df_fig$sym)
  #df_fig$sig[f] <- 'no'
  
  require('ggrepel')
  df_fig$nudge_x <- 0.2
  df_fig$nudge_x[f_tgfb48_rnaseq_dn_rpf_up] <- -0.2 
  df_fig$nudge_y <- 0.2
  #f <- grepl("Malat1", df_fig$sym)
  #df_fig[f, 'nudge_x'] <- 0.4  
  #df_fig[f, 'nudge_y'] <- 0.4
  
  # change name
  df_fig$sym <- mgsub::mgsub(df_fig$sym,
                             tolower(c('HALLMARK','_')), c('',' '))
  #df_fig$sym <- str_wrap(df_fig$sym, width=15)
  
  f <- df_fig$sig != 'no'
  if (any(f)) {
    df1 <- df_fig[f,,drop=F]          
    gg <- gg + geom_text_repel(data=df1,
                               aes(label=sym, lineheight=.75),
                               size=3.5, colour='black',
                               force=1, box.padding=0.25, point.padding=0.5,
                               min.segment.length = unit(0, 'lines'),          
                               nudge_x=df1$nudge_x, nudge_y=df1$nudge_y,
                               seed=seed_geom_text_repel)
  }
  
}

print_figure(gg, width=4, height=4,
             file=sprintf("scatter_plot.fig1%s.rpf_only_up_dn", fname_appendix))



##############################################################
##############################################################
#Fig 1
df_fig <- df_all
sym <- rownames(df_fig)
df_fig$fig.type <- NA

verb("sym_mrna_only: %d\n", length(sym_mrna_only))
f <- sym %in% sym_mrna_only
str_mrna_only <- sprintf("RNA Only (n=%s)", format(length(sym_mrna_only), big.mark = ","))
df_fig[f, "fig.type"] <- str_mrna_only

verb("sym_mrna_up_ribo_dn: %d\n", length(sym_mrna_up_ribo_dn))
verb("sym_mrna_dn_ribo_up: %d\n", length(sym_mrna_dn_ribo_up))
f <- sym %in% sym_mrna_up_ribo_dn | sym %in% sym_mrna_dn_ribo_up
str_mrna_rpf_discord <- sprintf("RNA RPF Discord (n=%s)", format(length(sym_mrna_up_ribo_dn) + 
                                                                   length(sym_mrna_dn_ribo_up), big.mark = ","))
df_fig[f, "fig.type"] <- str_mrna_rpf_discord

verb("sym_mrna_up_ribo_up: %d\n", length(sym_mrna_up_ribo_up))
verb("sym_mrna_dn_ribo_dn: %d\n", length(sym_mrna_dn_ribo_dn))
f <- sym %in% sym_mrna_up_ribo_up | sym %in% sym_mrna_dn_ribo_dn
str_mrna_rpf <- sprintf("RNA and RPF (n=%s)", format(length(sym_mrna_up_ribo_up) + 
                                                       length(sym_mrna_dn_ribo_dn), big.mark = ","))
df_fig[f, "fig.type"] <- str_mrna_rpf

verb("sym_mrna_no_ribo_up: %d\n", length(sym_mrna_no_ribo_up))
verb("sym_mrna_no_ribo_dn: %d\n", length(sym_mrna_no_ribo_dn))
f <- sym %in% sym_mrna_no_ribo_up | sym %in% sym_mrna_no_ribo_dn
str_rpf_only <- sprintf("RPF Only (n=%s)", format(length(sym_mrna_no_ribo_up) + length(sym_mrna_no_ribo_dn), 
                                                  big.mark = ","))
df_fig[f, "fig.type"] <- str_rpf_only

f <- !is.na(df_fig$fig.type) & !is.na(df_fig$log2FCunt48VStgfb48.transcription) & 
  !is.na(df_fig$log2FCunt48VStgfb48.translation)
df_fig <- df_fig[f, ]
df_fig$fig.type <- factor(df_fig$fig.type, level = c(str_mrna_only, str_mrna_rpf_discord, 
                                                     str_mrna_rpf, str_rpf_only))
df_fig <- df_fig[order(df_fig$fig.type), ]

nv_color <- c(`RNA Only` = "#aaaaaa", `RNA RPF Discord` = "#fcf050",
              `RNA and RPF` = "#7aaa3d", `RPF Only` = "#d78c5e")
names(nv_color) <- levels(df_fig$fig.type)



nticks=13; xmax=12; ymax=12; gap_tick_label=4
gg <- ggplot(data=df_fig,
             aes(x=log2FCunt48VStgfb48.transcription, 
                 y=log2FCunt48VStgfb48.translation, colour=fig.type)) + 
  theme_geometry(ticks=nticks, xlim=xmax, ylim=ymax, linesize=0.5, 
                 xlab=expression('RNA-seq log'[2]*'(fold change)'),
                 ylab=expression('RPF-seq log'[2]*'(fold change)'),
                 labsize=3.5, labgap=1.5, epsilon=max(xmax,ymax)/50, gap_tick_label=gap_tick_label) +
  geom_point(alpha=1, size=1) +
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=9),
        legend.background = element_rect(color = NA),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.4, "cm"),
        legend.key.width = unit(0,"cm"),
        #legend.spacing.y = unit(0.5, "cm"),          
        legend.position = c(0.06, 0.9), legend.justification = c(0, 1) ) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ggtitle("") +
  scale_colour_manual(values = nv_color)

if (f_display_sym) {
  
  df_fig$sig <- 'no'; 
  df_fig$sym <- rownames(df_fig)  
  # Snord43, Snora23, Snord83b, Snord17, Snord104  
  f_sym <- grepl('Malat|Snord43', df_fig$sym)
  
  th_log2fc <- log2(1.2)
  f_tgfb48_rnaseq_up_rpf_up <- !is.na(df_fig$fig.type) & df_fig$log2FCunt48VStgfb48.transcription > th_log2fc & df_fig$log2FCunt48VStgfb48.translation > th_log2fc
  f_tgfb48_rnaseq_dn_rpf_up <- !is.na(df_fig$fig.type) & df_fig$log2FCunt48VStgfb48.transcription < -th_log2fc & df_fig$log2FCunt48VStgfb48.translation > th_log2fc
  df_fig$sig[f_tgfb48_rnaseq_up_rpf_up & f_sym] <- 'up_up';
  df_fig$sig[f_tgfb48_rnaseq_dn_rpf_up & f_sym] <- 'dn_up';
  #f <- grepl("Gm|Rik", df_fig$sym)
  #df_fig$sig[f] <- 'no'
  
  require('ggrepel')
  df_fig$nudge_x <- 0.2
  df_fig$nudge_x[f_tgfb48_rnaseq_dn_rpf_up] <- -0.2 
  df_fig$nudge_y <- 0.2
  #f <- grepl("Malat1", df_fig$sym)
  #df_fig[f, 'nudge_x'] <- 0.4  
  #df_fig[f, 'nudge_y'] <- 0.4
  
  # change name
  df_fig$sym <- mgsub::mgsub(df_fig$sym,
                             tolower(c('HALLMARK','_')), c('',' '))
  #df_fig$sym <- str_wrap(df_fig$sym, width=15)
  
  f <- df_fig$sig != 'no'
  if (any(f)) {
    df1 <- df_fig[f,,drop=F]          
    gg <- gg + geom_text_repel(data=df1,
                               aes(label=sym, lineheight=.75),
                               size=3.5, colour='black',
                               force=1, box.padding=0.25, point.padding=0.5,
                               min.segment.length = unit(0, 'lines'),          
                               nudge_x=df1$nudge_x, nudge_y=df1$nudge_y,
                               seed=seed_geom_text_repel)
  }
  
}

print_figure(gg, width=4.1, height=4.1,
             file=sprintf("scatter_plot.fig1%s", fname_appendix))


## GO enrichment
#RNA only up
condstr <- "sym_rna_only_up"

sym_rna_only_up <- intersect(sym_mrna_only, sym_mrna_up)
verb("sym_rna_only_up: %d\n", length(sym_rna_only_up))


library(clusterProfiler)
library(DOSE)

entrez.id <- unique(entrezdf$entrez_id[ match(sym_rna_only_up, entrezdf$gene_name) ])

set.seed(40)
egmt_rna_only_up <- enrichGO(entrez.id, OrgDb=org.Mm.eg.db, keyType = "ENTREZID", ont='all',
                             pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,
                             minGSSize = 10, maxGSSize = 500)
egmt_rna_only_up <- setReadable(egmt_rna_only_up, OrgDb=org.Mm.eg.db, keyType="ENTREZID")

head(egmt_rna_only_up)
dim(egmt_rna_only_up)

write.table(egmt_rna_only_up, file = sprintf('table/fig1.%s.go.txt', condstr),
            row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )




##dotplot
n_top_egmt <- 10
gg <- dotplot(egmt_rna_only_up, x = "GeneRatio",
              color = "p.adjust", showCategory = n_top_egmt, split = "ONTOLOGY", font.size = 11,
              title = "RNA Only Up") + facet_grid(ONTOLOGY~., scale="free")

print_figure(gg, width=9, height=8.5,
             file=sprintf("dotplot.%s.enricher.go_top%d", condstr, n_top_egmt))





