suppressMessages(library(Biostrings))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(limma))
suppressMessages(library(cqn)) ## need install
suppressMessages(library(edgeR))
suppressMessages(library(methods))
suppressMessages(library(utils))
suppressMessages(library(stats))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(Hmisc))
suppressMessages(library(Matrix)) # rank
suppressMessages(library(gplots))
suppressMessages(library(cowplot))
suppressMessages(library(scales)) # muted
suppressMessages(library(GGally)) # scatmat .. need install
suppressMessages(library(data.table))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(fastmatch))

source("./differential_expression_analysis.TCC.R")
verb <- function(...) cat(sprintf(...), sep='', file=stdout())
cmat <- Exp


simple_limma_voom_for_mrna  <-
  function(       project ,
                  individual = NULL ,
                  conditions = NULL ,
                  nuisance = NULL ,
                  strategy = NULL ,
                  mtx_count = NULL ,
                  #normMethod = "TMM" ,
                  doVoom = FALSE ,
                  min.expr.counts = 5 ,
                  min.expr.log2cpm = -Inf ,
                  min.expr.num.samples = 2 ,
                  ebayes.trend = TRUE ,
                  ebayes.robust = TRUE ,
                  weights = TRUE ,
                  triples = FALSE ,
                  MHadjust = "BH" ,
                  FDR.thresh = 0.05 ,
                  th_log2fc = 0.26 ,
                  refColumn = NULL ,
                  rundate_appendix = "" ,
                  outfbase = 'TEST' ,
                  gtf = NULL ,
                  verbose = FALSE,
                  target_species = 'mm10' ) {
    
  
    if (! is.null(project)){
      prj <- read.table(project, header=T)
      myData.project  <-  process_project_input( project =prj)
    }
    print(myData.project)
    
    rownames(myData.project)  <-  myData.project$runid
    
    myData.project  <-  myData.project[ order(myData.project$order , myData.project$runid) , ,drop=FALSE]
    
    if (is.null(individual)) {
      individual  <-  myData.project$individual[1]
    }
    
    if (!is.null(strategy)) {
      myData.project  <-  myData.project[ myData.project$strategy == strategy , ,drop=FALSE]
    } # strategy
    
    if (!is.null(conditions)) {
      myData.project  <-  subset( myData.project , condition %in% conditions )
    } # conditions
    
    
    # basic info
    myData.species  <-  myData.project$species[1]
    myData.strategy  <-  myData.project$strategy[1]
    myData.DNAlevel  <-  get_DNA_level_for_strategy( strategy = myData.strategy )
    
    # check
    stopifnot( length(unique(myData.project$strategy)) == 1 )
    stopifnot( length(unique(myData.project$species)) == 1 )
    stopifnot( length(unique(myData.project$individual)) == 1 )
    

    ############################# counts
    #print(rownames(mtx_count))
    
    ### convert to gene names
    f_meta_tags <- grepl("^__|rDNA_promoter", rownames(mtx_count))
    if (any(f_meta_tags)) {
      mtx_meta  <-  mtx_count[f_meta_tags,,drop=F]
      mtx_count  <-  mtx_count[!f_meta_tags,]
    }
    fname  <-  sprintf("%s.counts.raw.gene_id.txt", outfbase )
    write.table( mtx_count , file = fname , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
    
    #mtx_count  <-  convert_ensembl_gene_id_to_gene_name_and_aggregate_matrix( mat = mtx_count , species = myData.species , fun = sum , gtf = myData.trangtf )
  
    if (any(f_meta_tags)) {
      mtx_count <- rbind(mtx_meta, mtx_count)
    }
    ## Modification by Nitish
    mtx_count <- log2(mtx_count+1)
    mtx_count <- mtx_count[!grepl("^NA$", rownames(mtx_count)),]
    #mtx_count <- mtx_count[!is.na(rownames(mtx_count)),]
    fname  <-  sprintf("%s.counts.raw.txt", outfbase )
    write.table( mtx_count , file = fname , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
    
    ########################### prepare
    verb("\tprepare.\n")
    verb("count matrix dimension:\n")
    show(dim(mtx_count))
    
    ################### design
    verb("\n\ndesign.\n")
    design  <-  make_design_matrix( info = myData.project , condition = "condition" , nuisance = nuisance , verbose=verbose )
    
    show(design)
    write.table( design , file = sprintf("%s.design.txt",outfbase) , sep = "\t" , quote = FALSE , row.names = TRUE , col.names = NA )
    
    ###################### contrast
    verb("\n\ncontrast.\n")
    contrast.matrix  <-  make_contrast_matrix( info = myData.project , condition = "condition" , design = design )
    
    write.table( contrast.matrix , file = sprintf("%s.contrast.txt",outfbase) , sep="\t" ,  row.names=TRUE , col.names=NA , quote = FALSE  )
    
    ########################## count preparation
    verb("count preparation.\n")
    # do not remove __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique, since # of total reads need to be transfered to DGEList, but the meta tag rows from the DGEList object should be removed before downstream analyses. https://www.biostars.org/p/379690/,https://support.bioconductor.org/p/104972/
    # https://rdrr.io/bioc/edgeR/man/DGEList.html
    # default: lib.size = colSums(counts): numeric vector giving the total count (sequence depth) for each library.
    
    message('[Condition]')
    print(myData.project$condition)
    myData.counts.dge <- DGEList( counts = mtx_count  ,  group = myData.project$condition  ,  genes = row.names(mtx_count) )
    
    ### filter
    verb("\tfiltering.\n")
    message("[filtering][min.expr.counts] ", min.expr.counts)
    message("[filtering][min.expr.num.samples] ", min.expr.num.samples)
    message("[filtering][min.expr.log2cpm] ", min.expr.log2cpm)
    
    logi.good.count.expr  <-  rowSums(getCounts(myData.counts.dge) >= min.expr.counts)  >=  min.expr.num.samples
    logi.good.log2cpm.expr  <-  rowSums(cpm(myData.counts.dge , log = TRUE) >= min.expr.log2cpm) >= min.expr.num.samples
    logi.is.expr  <-  logi.good.log2cpm.expr  &  logi.good.count.expr
    
    
    # begin of addition by H. Kim
    f_meta_tags <- grepl("^__|rDNA_promoter", rownames(myData.counts.dge))
    logi.is.expr <- logi.is.expr & !f_meta_tags
    logi.is.nonexpr <- !logi.is.expr & !f_meta_tags
    ngene <- length(logi.is.expr) - length(which(f_meta_tags))
    # end of addition
    
    
    verb("\t\t\tnon-expressed genes:  %d  /  %d   =  %2.2f%%\n" , sum(logi.is.nonexpr) , ngene , sum(logi.is.nonexpr) / ngene * 100 )
    
    ### save non-expressed
    #non.expr.names  <-  myData.counts.dge$genes$gene[!logi.is.expr]
    non.expr.names  <-  myData.counts.dge$genes$gene[logi.is.nonexpr]
    # end of modification
    
    myData.nonexpr.mat  <-  mtx_count[ non.expr.names ,]
    
    fname  <-  sprintf("%s.nonexpr.txt" , outfbase)
    write.table( myData.nonexpr.mat , file = fname , sep="\t" , row.names = TRUE , col.names = NA , quote = FALSE )
    
    #### expressed genes only
    # the meta tag rows from the DGEList object should be removed before downstream analyses. https://www.biostars.org/p/379690/,https://support.bioconductor.org/p/104972/
    myData.counts.filt <- myData.counts.dge[logi.is.expr,, keep.lib.sizes=FALSE]
    
    
    ################################# voom stabilization
    pdf( sprintf("%s.limma.plots.pdf",outfbase))
    if (doVoom) {
      verb("\n\nvoom stabilization.\n")
      verb("\tchecking voom.\n")
      the.Elist.voom  <-  voom( counts = myData.counts.filt.normalized , design = design , plot = TRUE )
      
      verb("\tvoomWithQualityWeights.\n")
      the.Elist.voom.qual  <-  voomWithQualityWeights( counts = myData.counts.filt.normalized , design = design , plot = TRUE , normalize.method = "none" , maxiter = 100 )
      
      if (weights) {
        the.Elist  <-  the.Elist.voom.qual
      } else {
        verb("\t\tproceeding WITHOUT weights.\n")
        the.Elist  <-  the.Elist.voom
      } # weights
      
    } else {
      the.Elist  <-  list( E = as.matrix(mtx_count) , genes = data.frame(genes = rownames(mtx_count), stringsAsFactors = FALSE) )
      the.Elist  <-  new("EList" , the.Elist)
    }
    
    ########################################### linear model analysis
    verb("linear model analysis.\n")
    
    
    ### fit
    verb("\tfit\n")
    fit <- lmFit( object = the.Elist , design = design)
    
    
    
    ### contrast fit
    verb("\tcontrast fit.\n")
    cfit <- contrasts.fit( fit = fit , contrasts = contrast.matrix )
    
    
    ### ebayes
    verb("\teBayes\n")
    cfit <- eBayes(fit = cfit , trend = ebayes.trend , robust = ebayes.robust )
    
    ########################### save to file
    verb("\tsave\n")
    write.fit( fit , file = sprintf("%s.fit" , outfbase) , digits = 8 , adjust = MHadjust , method = "separate" , sep = "\t"  )
    write.fit( cfit , file = sprintf("%s.cfit" , outfbase) , digits = 8 , adjust = MHadjust , method = "separate" , sep = "\t"  )
    
    write.table( fit$coefficients ,  file = sprintf("%s.log2cpm.txt",outfbase)   , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )
    write.table( fit$stdev.unscaled * fit$sigma ,  file = sprintf("%s.dev.txt",outfbase)   , row.names = TRUE , col.names = NA ,  sep="\t"  , quote = FALSE  )
    
    ################## write DE tables
    verb("\n\nwrite DE tables.\n")
    
    tu.gene.map  <-  NULL
    write_limma_de_tables(  cfit.fname = sprintf("%s.cfit" , outfbase) , fit = fit , info = myData.project , condition = "condition" , triples = triples , MHadjust = MHadjust ,
                            FDR.thresh = FDR.thresh , th_log2fc = th_log2fc , outfbase = outfbase , map = tu.gene.map ) # th_log2fc argument was added by H. Kim
    
    system(sprintf("gzip  -f  %s.*.txt", outfbase))
    system(sprintf("gzip  -f  %s*.cfit", outfbase))
    system(sprintf("gzip  -f  %s*.fit", outfbase))
    
    ################# limma plots
    verb("limma plots.\n")
    
    test.res  <-  limma::decideTests( cfit ,  adjust = MHadjust , method = "separate" , p.value = FDR.thresh )
    
    # expression values from all samples
    limma::plotDensities( the.Elist , log = TRUE  , main = "density of all log2cpm gene expression from all samples" )
    
    
    # MA,MD
    limma::plotMA(fit , main = "MA plot of fit")
    limma::plotMD(fit , main = "MD plot of fit")
    
    # MDS
    limma::plotMDS( the.Elist , main = "MDS plot")
    
    # heat
    limma::heatDiagram( test.res , coef = cfit$coefficients  )
    
    # venn
    uniqConds  <-  unique(myData.project$condition)
    if (length(uniqConds)  <= 3) {
      verb("\t\tvenn.\n")
      for  (direc  in  c("both","up","down"))  {
        limma::vennDiagram( test.res , include = direc )
      } # direc
    } # venn
    
    dummyPrint  <-  dev.off()
    verb("\n\n\nDone with differential expression analysis.\n\n\n")
  }


simple_limma_voom_for_mrna('./data/my_databases/Jake_Proteomics.txt',mtx_count=cmat, outfbase = 'Proteomics_data', th_log2fc = 0.26, target_species = 'mm39')
### Change by Nitish
## Change RData file name
save.image("WholeProteome.RData")
