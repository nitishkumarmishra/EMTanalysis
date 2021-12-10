#
#

# silencing the warning message of dplyr.summarise
options(dplyr.summarise.inform=F) 

# use pryr::object_size(), pryr:mem_used() rather than loading pryr.
#suppressMessages(suppressWarnings(library(pryr)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(Biostrings)))
# use qdap::mgsub() instead of loading qdap package
#suppressMessages(suppressWarnings(library(qdap)))
# do not use dplyr since pathview/eg2id function did not work with dplyr
# consider using magrittr for pipe %>% and dplyr::function() instead of loading dplyr
#suppressMessages(suppressWarnings(library(dplyr))) # mutate
suppressMessages(suppressWarnings(library(magrittr))) # for pipe %>%
# caution: the order of loading matrixStats after Rfast was intended since rowMaxs() needs to be overlayed.
# solution: try to avoid Rfast, but to use matrixStats
# suppressMessages(suppressWarnings(library(Rfast))) # for colMinsMaxs()
suppressMessages(suppressWarnings(library(matrixStats))) # for rowMaxs()
# id conversion
#suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(org.Mm.eg.db)))
suppressMessages(suppressWarnings(library(AnnotationDbi)))
# basic graph
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))

# enrichment analysis
suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(enrichplot)))

# display
suppressMessages(suppressWarnings(library(IRdisplay)))


### common parameters
dir_jupyter <- "."
f_reload_gmt <- FALSE


### color
c1<-brewer.pal(n=9, name = "Set1")
c1[6]<-'#fee287' # light yellow --> yellow
c2<-brewer.pal(n=8, name = "Set2")
c3<-brewer.pal(n=12, name = "Set3")
p1<-brewer.pal(n=9, name = "Pastel1")
p2<-brewer.pal(n=8, name = "Pastel2")

nv_row_annot_color <- c(

   #'unt48'=c1[2], 'tgfb48'=c1[1],
   #'tgfbWO'=c1[3], 'tgfbCX5461100nm'=c1[4],
   'unt48'="#00BA38",
   'tgfb48'="#F8766D", 
   'tgfbCX5461100nm'="#619CFF",

   'unt'="#00BA38",
   'tgfb'="#F8766D", 
   'tgfbCX'="#619CFF",

   'GNPs'=c1[4],
   'SHH model'="#ff0000",
    
   'Retro-Myc model'="#999900", # dd yellow
   'Retro-Myc model 19568-Cas9-Tumor'="#cccc00", # d yellow
   'CRISPR-Myc model'="#eeee00", # ~yellow
    
   'Utx-E3-Cre+p53DN+Myc'=c1[3], # green
   'Utx-E3-Cre+p53DN+Mycn'="#000000", # black
    
   # dark orange to light orange
   'Retro-Myc+Gfi1 model'=c1[5], # ff7f00,
   'WT-NS-Myc+Gfi1'="#ffaf30",
   'WT-NS-Myc+Gfi1 (secondary from 1232)'="#ffd700", # gold
    
   'WT-Prom1-Myc+Gfi1'=c1[7], # a65628
   'WT-Prom1-Myc+Gfi1 (secondary from 1266)'="#c07050",
   'WT-Prom1-Myc+Gfi1 (secondary from 57163)'="#d09070",
   'WT-NS or Prom1-Myc+Gfi1B'="#cccccc"

)


map_row_annot_color <- function(group_u, type_group="default", verbose=F) {
  group_u <- unique(as.character(group_u))

  group_u <- gsub('HALLMARK_', '', group_u)
  col_group <- nv_row_annot_color[group_u]
  names(col_group) <- group_u;
  idx_na <- which(is.na(col_group));
  n_na <- length(idx_na)
  if (n_na > 0) {
    # when any color was not mapped, use palette
    switch(type_group,
      'up'={col_group[idx_na] <- colorRampPalette(brewer.pal(9,"Reds"))(n_na+3)[(n_na+1):2]},
      'dn'={col_group[idx_na] <- colorRampPalette(brewer.pal(9,"Blues"))(n_na+3)[2:(n_na+1)]},
      {col_group[idx_na] <- colorRampPalette(brewer.pal(9,"Blues"))(n_na+3)[2:(n_na+1)]}
    )
  }

  if (verbose) {
    display(col_group)
  }
  return(col_group)

}    



# add_ortholog_info
# add columns of ortholog information (e.g., HomoloGene.ID, mouse.sym, mouse.eid, human.sym, human.eid)
add_ortholog_info <- function(df, from="mouse", to="human") {

  df_from <- orthodf[orthodf$specie==from,]
  idx <- match(rownames(df), df_from[,"Symbol"])
  df$HomoloGene.ID <- df_from[idx,"HomoloGene.ID"]
  df[, sprintf("%s.sym",from)] <- df_frome[idx,"Symbol"]
  df[, sprintf("%s.eid",from)] <- df_from[idx,"EntrezGene.ID"]

  df_to <- orthodf[orthodf$specie==to,]
  idx <- match(df$HomoloGene.ID , df_to$HomoloGene.ID)
  df[, sprintf("%s.sym",to)] <- df_to[idx,"Symbol"]
  df[, sprintf("%s.eid",to)] <- df_to[idx,"EntrezGene.ID"]

  return(df)

}


### heatmap
suppressMessages(suppressWarnings(library(ComplexHeatmap)))
suppressMessages(suppressWarnings(library(circlize))) # chordDiagram()

get_heatmap <- function(df, vec_color_group, genes=NULL, type_col="zscore", col=NULL,
    platform="rnaseq", type_heatmap="details",
    row_pattern=NULL, column_pattern=NULL,
    row_title=NULL, row_title_side="left",row_title_gp=gpar(fontsize=10),
    column_title=NULL, column_title_side="top", column_title_gp=gpar(fontsize=10),
    show_column_names=T, show_row_names=T, show_heatmap_legend=T,
    column_annot_title_gp=gpar(fontsize=9), column_annot_labels_gp=gpar(fontsize=6),
    fontsize_row_names=6,fontsize_column_names=5,
    cluster_rows=T, clustering_distance_rows="euclidean", clustering_method_rows="complete",
    cluster_columns=F, clustering_distance_columns="euclidean", clustering_method_columns="complete",
    column_annot_title="group", column_annot_ncol=1, sample_names=NULL, column_name_angle=NULL, ...) {

  # special treatment
  colnames(df) <- mgsub::mgsub(colnames(df),
         c("tss", "_plusminus_", "_minus_", "_plus_", "five_prime_utr|5UTR", "CDS", "three_prime_utr|3UTR"),
         c("TSS", " ±", " -", " +", "5'UTR", "CDS", "3'UTR"))
  if (!is.null(row_pattern)) df <- df[grepl(row_pattern, rownames(df), perl=T), , drop=F]
  if (!is.null(column_pattern)) df <- df[, grepl(column_pattern, colnames(df),perl=T), drop=F]

  top_annotation <- NULL
  top_annotation_height <- NULL
  if (!is.null(vec_color_group)) {
    # column annotation
    df_annot <- data.frame(group=vec_color_group, row.names=colnames(df))
    col_group <- map_row_annot_color(vec_color_group)
    top_annotation=HeatmapAnnotation(df=df_annot,col=list(group=col_group),
        annotation_legend_param=list(
            group=list(title=column_annot_title,title_gp=column_annot_title_gp, labels_gp=column_annot_labels_gp, grid_height=unit(3,"mm"), ncol=column_annot_ncol)
         ))
    top_annotation_height <- top_annotation@size
  }

  bottom_annotation <- NULL
  bottom_annotation_height <- NULL
  if ((!is.null(column_name_angle)) && (show_column_names)) {
    # default: angle=90, bottom to top
    just="right"
    bottom_annotation_height <- max_text_width(colnames(df), gp=gpar(fontsize=fontsize_column_names))
    if (column_name_angle == -90) {
       # angle=-90, top to bottom
       just='left'
    } else if (column_name_angle == 0) { just='top';
       bottom_annotation_height <- max_text_height(colnames(df), gp=gpar(fontsize=fontsize_column_names))
    }
    bottom_annotation=HeatmapAnnotation(cn=anno_text(colnames(df),'column', gp=gpar(fontsize=fontsize_column_names), rot=column_name_angle, just=just, offset=unit(1,"npc")) )
    show_column_names <- F
  }

  type_col_short <- tolower(gsub('[()]','',type_col))
  if (!is.null(genes)) {
    f <- !is.na(genes) & genes %in% rownames(df)
    genes <- genes[f]
    df <- df[genes,,drop=F]
  } else if (nrow(df) > 5000) {
    return(list(type_col_short=type_col_short, list_ht=NULL))
  }

  # df_score
  df_score <- df

  switch(type_col_short,
    'zscore' = {
      df_score <- t(scale(t(df_score),center=T,scale=T))
      col_score <- colorRamp2(c(-2,0,2),c("#0000ff", "#ffffcc","#ff0000")) },   
    'log2cpm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },
    'log2fpkm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },         
    '-log10pvalue'={
      type_col <- expression('-log'[10]*'(pvalue)')
      col_score <- colorRamp2(c(0,4),c("#ffffcc","#ff0000")) },
    'dlog10p'={
      col_score <- colorRamp2(c(-4,0,4),c("#0000ff", "#ffffcc","#ff0000")) },   
    'gsva'={
      col_score <- colorRamp2(c(-4,0,4),c("#0000ff", "#ffffcc","#ff0000")) },   
    'ssgsea'={
      col_score <- colorRamp2(c(-0.5,0,0.5),c("#0000ff","#ffffcc","#ff0000")) },
    {}
  )
 
  if (!is.null(col)) col_score <- col
    
  switch(type_heatmap,
    "overview"={
       rect_gp=gpar(col=NA)
     },
     "details"={
       rect_gp=gpar(col="white", lty=1, lwd=1)
     },
     {}
  )
    
  if (!is.null(sample_names)) colnames(df_score) <- sample_names
  ht11=Heatmap(df_score, col=col_score,
        cluster_rows=cluster_rows,
        clustering_distance_rows=clustering_distance_rows,
        clustering_method_rows=clustering_method_rows,
        cluster_columns=cluster_columns,
        clustering_distance_columns=clustering_distance_columns,
        clustering_method_columns=clustering_method_columns,
        show_column_names=show_column_names, show_row_names=show_row_names,
        row_names_gp=gpar(fontsize=fontsize_row_names),
        column_names_gp=gpar(fontsize=fontsize_column_names),
        rect_gp=rect_gp,
        heatmap_legend_param=list(title=type_col,
            title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=9)),
        row_title=row_title,
        row_title_side=row_title_side,
        row_title_gp=row_title_gp,
        row_title_rot = 90,
        column_title=column_title,
        column_title_side="top",
        column_title_gp=column_title_gp,
        top_annotation=top_annotation, top_annotation_height=top_annotation_height,
        bottom_annotation=bottom_annotation,
        bottom_annotation_height=bottom_annotation_height,
        show_heatmap_legend=show_heatmap_legend,...)

  result <- list()
  result$type_col_short <- type_col_short
  result$list_ht <- ht11
  return(result)
}



get_heatmap_row_annotation <- function(df, vec_color_group, nv_genes=NULL,
        vec_gene_type=NULL, type_col="zscore", nv_row_annot=NULL,
        platform="rnaseq", type_heatmap="details",
        row_pattern=NULL, column_pattern=NULL,
        row_title=NULL, row_title_side="left", row_title_gp=gpar(fontsize=10),
        column_title=NULL, column_title_side="top", column_title_gp=gpar(fontsize=10),
        show_column_names=T, show_row_names=T, show_heatmap_legend=T,
        row_annot_bar_width=unit(0.1,"in"), row_annot_max_len=15, row_annot_f_conv=F,
        row_annot_f_remove_first_word=F,row_annot_offset=-0.1,nv_sym_up=NULL,nv_sym_dn=NULL,
        column_annot_title_gp=gpar(fontsize=9), column_annot_labels_gp=gpar(fontsize=6),
        fontsize_row_names=6, fontsize_column_names=5, fontsize_row_annot=8,
        cluster_rows=T, clustering_distance_rows="euclidean", clustering_method_rows="complete",
        cluster_columns=F, clustering_distance_columns="euclidean", clustering_method_columns="complete",
        column_annot_title="group", column_annot_ncol=1, sample_names=NULL,...) {


  # special treatment
  colnames(df) <- mgsub::mgsub(colnames(df),
	 c("tss", "_plusminus_", "_minus_", "_plus_", "five_prime_utr|5UTR", "CDS", "three_prime_utr|3UTR"),
	 c("TSS", " ±", " -", " +", "5'UTR", "CDS", "3'UTR"))
  if (!is.null(row_pattern)) df <- df[grepl(row_pattern, rownames(df), perl=T), , drop=F]
  if (!is.null(column_pattern)) df <- df[, grepl(column_pattern, colnames(df), perl=T), drop=F]
  
  # column annotation
  df_annot <- data.frame(group=vec_color_group, row.names=colnames(df))
  col_group <- map_row_annot_color(vec_color_group)
  top_annotation <- HeatmapAnnotation(df=df_annot,col=list(group=col_group),
    annotation_legend_param=list(
      group=list(title=column_annot_title, title_gp=column_annot_title_gp, labels_gp=column_annot_labels_gp, grid_height=unit(3,"mm"), ncol=column_annot_ncol)
      ))

  type_col_short <- tolower(gsub('[()]','',type_col))
  if (!is.null(nv_genes)) {
    f <- !is.na(nv_genes) & nv_genes %in% rownames(df)
    nv_genes <- nv_genes[f]
    df <- df[nv_genes,,drop=F]
    if (!is.null(vec_gene_type)) vec_gene_type <- vec_gene_type[f]
  } else if (nrow(df) > 5000) {
    return(list(type_col_short=type_col_short, list_ht=NULL))
  }

  # colors for row annotation
  if (!is.null(vec_gene_type)) {
    # gene groups
    df_annot <- data.frame(group=vec_gene_type, row.names=rownames(df))
    group_u <- unique(vec_gene_type); n_group <- length(group_u)
    col_group <- map_row_annot_color(group_u)
    vec_split <- vec_gene_type
    row_title_gp <- gpar(fontsize=0) # suppress the original row title
  } else {
    vec_gene_type <- names(nv_genes)
    df_annot <- data.frame(group=vec_gene_type, row.names=rownames(df))
    col_group_up <- map_row_annot_color(unique(names(nv_sym_up)),'up')
    col_group_dn <- map_row_annot_color(unique(names(rev(nv_sym_dn))),'dn')
    col_group <- c(col_group_up, col_group_dn)
    vec_split <- NULL
  }
    
  # row annotation
  ha12 <- rowAnnotation(df=df_annot, col=list(group=col_group),
    annotation_legend_param=list( group=list(title="gene set",
      title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=5))),
        width=row_annot_bar_width, show_legend=F)
   
  # df_score
  df_score <- df    
    
  switch(type_col_short,
    'zscore' = {
      df_score <- t(scale(t(df_score),center=T,scale=T))
      col_score <- colorRamp2(c(-2,0,2),c("#0000ff", "#ffffcc","#ff0000")) },
    'log2cpm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },
    'log2fpkm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },            
    {}
  )
    
  switch(type_heatmap,
    "overview"={
       rect_gp=gpar(col=NA)
     },
     "details"={
       rect_gp=gpar(col="white", lty=1, lwd=1)
     },
     {}
  )

  if (!is.null(sample_names)) colnames(df_score) <- sample_names
  ht11 <- Heatmap(df_score, col=col_score,
        split=vec_split,
        cluster_rows=cluster_rows,
        clustering_distance_rows=clustering_distance_rows,
        clustering_method_rows=clustering_method_rows,
        cluster_columns=cluster_columns,
        clustering_distance_columns=clustering_distance_columns,
        clustering_method_columns=clustering_method_columns,
        show_column_names=show_column_names, show_row_names=show_row_names,
        row_names_gp=gpar(fontsize=fontsize_row_names),
        column_names_gp=gpar(fontsize=fontsize_column_names),
        rect_gp=rect_gp,
        heatmap_legend_param=list(title=type_col,
            title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=9)),
        row_title=row_title,
        row_title_side=row_title_side,
        row_title_gp=row_title_gp,
        row_title_rot = 90,
        column_title=column_title,
        column_title_side="top",
        column_title_gp=column_title_gp,
        top_annotation=top_annotation,
        show_heatmap_legend=show_heatmap_legend,...)
    
  # row annotation text
  ht_split <- ht11
  if (is.null(vec_split)) ht_split <- NULL

  # row annotation
  row_annot_text <- rep('',1,length(vec_gene_type))
  idx <- match(nv_row_annot, rownames(df_score))
  row_annot_text[idx] <- names(nv_row_annot)

  ha13 <- rowAnnotation(
         id=anno_text(row_annot_text,'row', gp=gpar(fontsize=fontsize_row_annot),
                rot=0, just="left", offset=row_annot_offset),
         #width=unit(0.5,"cm")+max_text_width(row_annot_text, just="left", gp=gpar(fontsize=fontsize_row_annot)),
         width=max_text_width(row_annot_text, just="left", gp=gpar(fontsize=fontsize_row_annot)))

  # result
  result <- list()
  result$type_col_short <- type_col_short
  result$list_ht <- ht11 + ha12 + ha13
  return(result)
    
}




# heatmap_pca
# input:
#   list_pca:
#     df_loading: df_loading <- as.data.frame(prcomp_obj$rotation)
#   str_pc: {'PC1','PC2','PC3',...}
#   df_log2cpm:
#   vec_color_group:
#   n_pos:
#   n_neg:
#   fname_table:
#   f_format:
# output:
#   list_out$df_loading_sorted
# usage:
# list_out <- heatmap_pca(df_loading, 'PC2', df_log2cpm, vec_color_group=project$condition, n_pos=10, n_neg=10, cluster_rows = F, column_title='', column_annot_ncol=3, fontsize_row_names=6, fontsize_column_names=6, column_name_angle=45)
# print_figure(list_out$list_ht, width = 3.25, height = 3.25, file = sprintf("heatmap.pca2_%s", list_out$type_col_short))
# df <- cbind(list_out$df_dn[,'PC2',drop=F], df_log2cpm[rownames(list_out$df_dn),])
# cols <- colnames(df)
# df[,cols] <- sapply(cols, function(x) { formatC(df[,x], digits=2) })
# write.table(df, file = "table/pca_pc2_170224+190122.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
#
#
heatmap_pca <- function(list_pca, str_pc, df_log2cpm, vec_color_group, n_pos=10, n_neg=10, fname_table=NULL, f_formatc=T, ...) {

  if (is.list(list_pca)) {
    df_loading <- list_pca$df_loading
  } else {
    df_loading <- list_pca
  }

  df_loading_sorted <- df_loading[order(df_loading[,str_pc], decreasing=T),]
  pc_pos <- head(rownames(df_loading_sorted), n_pos)
  pc_neg <- tail(rownames(df_loading_sorted), n_neg)

  list_out <- get_heatmap(df_log2cpm, vec_color_group = vec_color_group,
      genes=c(pc_pos, pc_neg), type_col = "zscore", ...)
  list_out$df_loading_sorted <- df_loading_sorted

  if (!is.null(fname_table)) {
    df <- cbind(df_loading_sorted[,str_pc,drop=F], df_log2cpm[rownames(df_loading_sorted),])
    if (f_formatc) {
      cols <- colnames(df)
      df[,cols] <- sapply(cols, function(x) { formatC(df[,x], digits=2) })
    }
    write.table(df, file = fname_table, row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }

  return(list_out)

}



# theme_geometry
# original function from https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2
# modified by H. Kim
# input:
#   xvals: values of x that will be plotted
#   yvals: values of y that will be plotted
#   xgeo: x intercept value for y axis
#   ygeo: y intercept value for x axis
#   color: default color for axis
#   linesize: line size for axis
#   xlab: label for x axis
#   ylab: label for y axis
#   labsize: label size for axis
#   labgap: gap betwen axis and axis label
#   ticks: number of ticks to add to plot in each axis
#   textsize: size of text for ticks
#   xlimit: limit value for x axis
#   ylimit: limit value for y axis
#   epsilon: parameter for small space
#   gap_tick_label: gap betwen axis and tick label
#
# usage:
# see fig2.ipynb or r_ggplot_scatterplot.txt
#
theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
        color = "black", linesize = 1, 
        xlab = "x", ylab = "y", labsize=3.5, labgap=2,
        ticks = 10,
        textsize = 3, 
        xlimit = max(abs(xvals),abs(yvals)),
        ylimit = max(abs(yvals),abs(xvals)),
        epsilon = max(xlimit,ylimit)/50, gap_tick_label=3){

  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))

  #Add axis
  theme.list <- 
  list(
    theme_void(), #Empty the current theme
    geom_line(aes(x = x_ax, y = y_ax), color = color, size = linesize, data = xaxis),
    geom_line(aes(x = x_ax, y = y_ax), color = color, size = linesize, data = yaxis),
    # begin of modification by H. Kim
    #annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
    #annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
    annotate("text", x = 0, y = -ylimit-labgap, label = xlab, size = labsize),
    annotate("text", x = -xlimit-labgap, y = 0, label = ylab, size = labsize, angle=90),
    # end of modification
    xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
    ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
  )

  #Add ticks programatically
  ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)

  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){

    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))

    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))

    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                         data = xtick, size = linesize, 
                                         color = color)

    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - gap_tick_label*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))


    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = linesize, 
                                             color = color)

    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                            x = xgeo - gap_tick_label*epsilon, 
                                            y = ticks_y[k],
                                            size = textsize,
                                            label = paste(ticks_y[k]))
  }

  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}



print_figure <- function(obj, width, height, file=NULL, graph_format="tiff", resolution=300, f_display2screen=T) {

  if (is.null(obj)) return()

  type_obj=class(obj)[1]

  if (!is.null(file)) {
    # save figure to a file
    switch(graph_format,
	"tiff"={ 
    		filename <- sprintf('tiff/%s.tiff', file);
    		tiff(filename, width=width, height=height, units='in', res=resolution, compression='lzw', type='cairo')
	},
	"pdf"={
    		filename <- sprintf('pdf/%s.pdf', file);
    		cairo_pdf(filename, width=width, height=height)      
	},
	{}
    )
    switch(type_obj,
      "Heatmap"={
         obj <- draw(obj, heatmap_legend_side="left",
                          annotation_legend_side="bottom")
      },
      "HeatmapList"={
         obj <- draw(obj, heatmap_legend_side="left",
                         annotation_legend_side="bottom")
      },
      "list"={
         # print gg
         obj <- obj$gg
      },
      { }
    )
    suppressMessages(suppressWarnings(print(obj)))
    dev.off()  
  } 

  if (f_display2screen) {
    options(repr.plot.width=width, repr.plot.height=height)
    switch(type_obj,
      "Heatmap"={
         draw(obj, heatmap_legend_side="left",
                         annotation_legend_side="bottom")
      },
      "HeatmapList"={
         draw(obj, heatmap_legend_side="left",
                         annotation_legend_side="bottom")
      },
      "list"={
         grid.draw(obj$gt)
      },
      { 
	suppressMessages(suppressWarnings(print(obj)))
      }
    )
  } # f_display2screen

  #IRdisplay::display_html('')  
  cat(sprintf(""))
}



##### dge

# usage:
# list_venn <- list()
# list_venn[['ours']] <- sym_up; list_venn[['elife2018']] <- sym_up_nmumg
# vec_color <- c(c1[2],c1[1],c1[4])
# gg <- get_euler_diagram_from_list(list_venn, shape='circle', vec_color=vec_color, fills=list(fill=vec_color, alpha=0.6), fill_alpha=0.6, edges=list(), labels=list(cex=0.8), quantities=list(cex=0.7), strips=list(), legend=FALSE, label_pos_x=NULL, label_pos_y=NULL, main='up genes in in uint48 vs. tgfb48 (transcription)', main.gp=list(cex=0.75))
# print_figure(gg, width=3.1, height=2.8, file='euler_diagram_up_genes')
#
get_euler_diagram_from_list <- function(list_venn, shape='circle', vec_color=NULL, fills=NULL, fill_alpha=0.6, edges=list(), labels=list(cex=0.8), quantities=list(cex=0.7), strips=list(), legend=FALSE, label_pos_x=NULL, label_pos_y=NULL, quantities_pos_x=NULL, quantities_pos_y=NULL, main=NULL, main.gp=list(cex=0.75), ...) {

  mtx <- conv_list2mtx(list_venn)
  if (is.null(fills)) {
    vec_color <- get_colors_from_list(list_venn, vec_color)
    fills <- list(fill=vec_color, alpha=fill_alpha)
  }

  require(eulerr)
  fit <- euler(mtx, shape=shape)
  gg <- plot(fit, fills = fills,
           edges = edges,
           labels = labels,
           quantities = quantities,
           strips = strips, legend = legend, main = main, ...)

  if (!is.null(main)) {
    t <- getGrob(gg,'main.grob')
    gg <- editGrob(gg, 'main.grob', gp=main.gp)
  }

  if (!is.null(label_pos_x)) {
    #grid.ls(gg)
    t <- getGrob(gg,'labels.grob')
    x <- t$x; y <- t$y
    for (i in 1:length(label_pos_x)) {
      x[i] <- t$x[i] + unit(label_pos_x[i], "mm")
      y[i] <- t$y[i] + unit(label_pos_y[i], "mm")
    }
    gg <- editGrob(gg,"labels.grob", x=x, y=y)
  }

  if (!is.null(quantities_pos_x)) {
    #grid.ls(gg)
    t <- getGrob(gg,'quantities.grob')
    x <- t$x; y <- t$y
    for (i in 1:length(quantities_pos_x)) {
      x[i] <- t$x[i] + unit(quantities_pos_x[i], "mm")
      y[i] <- t$y[i] + unit(quantities_pos_y[i], "mm")
    }
    gg <- editGrob(gg,"quantities.grob", x=x, y=y)
  }


  return(gg)
}



#https://github.com/cran/VennDiagram/blob/master/R/hypergeometric.test.R
# This function performs the hypergeometric test on the two categories. Taken from package BoutrosLab.statistics.general
# reference:
# http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
#
# usage:
# (1) under-representation: calculate.overlap.and.pvalue(sym_mrna_up, sym_mrna_up_elife, total.size=n.total.coding.genes, lower.tail = TRUE, adjust=FALSE) # under-representation, use lower.tal=TRUE (default) and adjust=FALSE, since probabilities are P[X ≤ x].
# (2) over-representation: calculate.overlap.and.pvalue(sym_mrna_up, sym_mrna_up_elife, total.size=n.total.coding.genes, lower.tail = FALSE, adjust=TRUE) # over-representation, use lower.tail=FALSE, adjust=TRUE, subtract x by 1, when P[X ≥ x] is needed.
calculate.overlap.and.pvalue = function(list1, list2, total.size, lower.tail = TRUE, adjust = FALSE) {

        # calculate actual overlap
        actual.overlap <- length(intersect(list1, list2));

        # calculate expected overlap
        # need to cast to avoid integer overflow when length(list1) * length(list2) is extremely large
        expected.overlap <- as.numeric(length(list1)) * length(list2) / total.size;

        adjust.value <- 0;

        # adjust actual.overlap to reflect P[X >= x]
        if (adjust & !lower.tail) {
                adjust.value <- 1;
                #warning('Calculating P[X >= x]');
                }

        # calculate significance of the overlap
        overlap.pvalue <- phyper(
                q = actual.overlap - adjust.value,
                m = length(list1),
                n = total.size - length(list1),
                k = length(list2),
                lower.tail = lower.tail
                );

        # return values
        return( c(actual.overlap, expected.overlap, overlap.pvalue) );

}


# https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
fisher.test.overlap.and.pvalue = function(list1, list2, total.size) {

    #contingency table:
    mtx <- matrix(c(total.size - length(union(list1,list2)),
                    length(setdiff(list1,list2)),
                    length(setdiff(list2,list1)),
                    length(intersect(list1,list2))), nrow=2)

    return( fisher.test(mtx, alternative="greater"))

}




##### utility

# row.match
# input:
#   nomatch: the value to be returned in the case when no match is found.
# output:
#   see help(match)
# reference:
# https://github.com/tagteam/prodlim/blob/master/R/row.match.R
#
row.match <- function(x, table, nomatch=NA){

  if (class(table)=="matrix") table <- as.data.frame(table)
  if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
  cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
  ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
  match(cx, ct, nomatch=nomatch)

}

conv_list2mtx <- function(list_strings) {

  sym <- Reduce(union, list_strings)
  n_elements <- length(list_strings)
  mtx <- matrix(FALSE, length(sym), n_elements,
           dimnames=list(sym,names(list_strings)))
  for (i in 1:n_elements) {
     f <- sym %in% list_strings[[i]]
     mtx[,i] <- f
  }

  return(mtx)
}


##' This function retrieves the indices of non-zero elements in sparse matrices
##' of class dgCMatrix from package Matrix. This function is largely inspired from 
##' the package \code{Ringo}
##' 
##' @title Retrieve the indices of non-zero elements in sparse matrices
##' @param x A sparse matrix of class dgCMatrix
##' @return A two-column matrix
##' @author Samuel Wieczorek
##' @examples
##' library(Matrix)
##' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow=TRUE, sparse=TRUE)
##' res <- nonzero(mat)
nonzero <- function(x){
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix
    
    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
        return(matrix(0, nrow=0, ncol=2,
                      dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}

# matchlast
# reference:
# https://stackoverflow.com/questions/37404084/return-last-match-from-vector
matchlast <- function (needles, haystack) {
   length(haystack) + 1L - match(needles, rev(haystack))
}



# split_string_with_length
# reference: 
# https://stackoverflow.com/questions/11619616/how-to-split-a-string-into-substrings-of-a-given-length
split_string_with_length <- function(text, n) {
    substring(text, seq(1, nchar(text)-(n-1), n), seq(n, nchar(text), n))
}



# conv_codon2aa
conv_codon2aa <- function(vec_codon, f_aaa=FALSE, f_append_aa=FALSE, f_append_codon=FALSE) {

  mycode <- GENETIC_CODE
  #Selenocysteine  Sec/U   UGA     https://en.wikipedia.org/wiki/Selenocysteine
  mycode[["TGA"]] <- "U"
  #Pyrrolysine     Pyl/O   UAG     https://en.wikipedia.org/wiki/Pyrrolysine
  mycode[["TAG"]] <- "O"

  suppressWarnings( vec_aa  <-  Biostrings::translate( DNAStringSet(vec_codon), genetic.code=mycode, no.init.codon=FALSE, if.fuzzy.codon="error" ) )

  if (f_aaa) {
    vec_aa <- conv_aa2aaa(vec_aa, f_append_aa=f_append_aa)
  } 
  if (f_append_codon) {
    paste(vec_aa, vec_codon, sep="_")
  } else {
    vec_aa
  }

} # conv_codon2aa


# conv_aa2aaa
conv_aa2aaa <- function(vec_aa, f_append_aa=FALSE) {

  vec_aa <- as.character(vec_aa)
  # convet amino-acid one-letter code into the three-letter one.
  suppressWarnings( vec_aaa <- aaa(vec_aa) )

  # 21st and 22nd amino acids, alternative stop codons
  idx <- match("U", vec_aa)
  if (length(idx) > 0) vec_aaa[idx] <- "Sec"
  idx <- match("O", vec_aa)
  if (length(idx) > 0) vec_aaa[idx] <- "Pyl"

  if (f_append_aa) {
    paste(vec_aaa, vec_aa, sep="")
  } else {
    vec_aaa
  }

} # conv_aa2aaa



# verb
verb <- function(...) cat(sprintf(...), sep='', file=stdout())






