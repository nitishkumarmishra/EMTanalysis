


suppressMessages(suppressWarnings(source("./differential_expression_analysis.R")))
suppressMessages(suppressWarnings(source("./functions_homer.R")))


### functions_transcriptome




get_transcriptome_feature_fasta_filename  <-  function( species , feat ) {

    items <- strsplit(feat, "_")[[1]]
    for (i in length(items):1) {

        file  <-  "./out/transcriptome-gtf/m38/Mus_musculus.GRCm38.97.rdna_rn18s.transcriptome.gtf.gz"
        basef  <-  gsub( "\\.gtf\\.gz$", "" , file)
        #featf  <-  paste( basef , feat , "fa.gz" , sep=".")
        featf  <-  paste( basef , paste(items[1:i], collapse='_') , "fa.gz" , sep=".")

        if (file.exists(featf)) break
    }

    return(featf)

} # get_transcriptome_feature_fasta_filename






### GTF

get_gene_biotypes_and_protein_coding  <-  function( gtf ) {
# returns a list with fields "biotype" and "proteinCoding"
#       "biotype"  =  character vector of unique biotypes from gtf
#       "proteinCoding" = character vector of transcript ids that are protein coding

        biot  <-  get_gtf_attribute_field( gtf = gtf , field = "gene_biotype" )
        myData.geneBiotype  <-  data.frame( gene_id = gtf$gene_id , biotype = biot )
        myData.geneBiotype  <-  unique(myData.geneBiotype)

        myData.proteinCoding  <-  unique(myData.geneBiotype$gene_id[myData.geneBiotype$biotype == "protein_coding"])

        return(list(biotype = myData.geneBiotype , proteinCoding = myData.proteinCoding))

} # get_gene_biotypes_and_protein_coding



get_canonical_transcript_per_gene_name  <-  function( gtf ) {
## DESCRIPTION
#       This function will convert "gene names/symbol" to a canonical transript id.
# The problem is that, outrageously, gene name to gene id (ensembl) is one-to-many for some gene names.
# Thus, while gene id to canonical gene transcript is one-to-one, to get from gene name to canonical gene
# transcript is nontrivial.
#
# INPUT
#       gtf             = transcriptome gtf (data frame)
#       canon           = a data frame containing canonical transcript_id per gene_id.
#                        if NULL, then will be inferred from the gtf.
#
# OUTPUT
#       dataframe with columns:
#                gene_id , transcript_id , gene , CDS.length , five_prime_utr.length , three_prime_utr.length


        gtf  <-  prepare_transcriptome_gtf( gtf )
        # 'seqname' 'source' 'feature' 'start' 'end' 'score' 'strand' 'frame' 'attribute' 'transcript_id' 'gene_name' 'gene_id'
        canon  <-  get_canonical_transcripts(gtf)
        utab  <-  unique(gtf[,c("gene_name","gene_id","transcript_id") ,drop=FALSE])

        gdf  <-  merge( utab , canon , by = c("gene_id","transcript_id") )

        gdf  <-  gdf[ order(gdf$CDS.length,gdf$five_prime_utr.length,gdf$three_prime_utr.length) , ,drop=FALSE]
        gdf  <-  gdf[ !duplicated(gdf$gene_name) , ,drop=FALSE]

        return(gdf)

} # get_canonical_transcript_per_gene_name




get_canonical_transcripts  <-  function( gtf ) {
# returns a data frame with 3 columns:  "gene_id", "transcript_id", "cds.length".
# each gene appears only once.

        # order by longest CDS, 5'UTR, 3'UTR (in that order)
        gtf$feat.length  <-  gtf$end - gtf$start + 1
        gtf  <-  gtf[ gtf$feature %in% c("five_prime_utr","three_prime_utr","CDS") , ,drop=FALSE]
        gtf  <-  gtf[ order(gtf$feature, -gtf$feat.length) , ,drop=FALSE]

        # canonincal set
        uniqdf  <-  gtf[ !duplicated(gtf$gene_id) , c("gene_id","transcript_id") ,drop=FALSE]

        # feat length table
        tfeat  <-  reshape2::dcast( gtf , transcript_id ~ feature , value.var = "feat.length" , fill = 0 )
        loc.o  <-  which(colnames(tfeat)  %in%  c("five_prime_utr","three_prime_utr","CDS"))
        colnames(tfeat)[loc.o]  <-  paste( colnames(tfeat)[loc.o] , "length" , sep=".")

        # put feature lengths onto canonical df
        uniqdf  <-  merge( uniqdf , tfeat , by = "transcript_id" , all.x = TRUE )

        return(uniqdf)

} # get_canonical_transcripts




get_all_transciptome_classification_info  <-  function( gtf ) {

	utab  <-  unique(gtf[,c("gene_name","gene_id","transcript_id") ,drop=FALSE])
        biot  <-  get_gene_biotypes_and_protein_coding( gtf = gtf )
        geneBiotype  <-  biot$biotype
        # consider only protein coding genes...
        proteinCodingGenes  <-  biot$proteinCoding

        canonTrans  <-  get_canonical_transcripts( gtf = gtf )
        protCod.canon  <-  canonTrans[ canonTrans$gene_id %in% proteinCodingGenes , ,drop=FALSE]

        canonTransName  <-  get_canonical_transcript_per_gene_name( gtf = gtf )

        retval  <-  list( utab = utab,
			geneBiotype = geneBiotype ,
                        proteinCodingGenes = proteinCodingGenes ,
                        canonicalTranscripts = canonTrans ,
                        canonicalTranscriptsGeneName = canonTransName ,
                        proteinCodingCanonTrans = protCod.canon )

        return(retval)

} # get_all_transciptome_classification_info


# build_transcriptome_info
# input:
#   species:
# output:
#   list$myData.transcriptome.gtf: data.frame of columns ("seqname", "source", "feature", "start", "end", "frame", "attribute", "transcript_id", "gene_name", "gene_id")
#      seqname: transcript
#   list$myData.transcriptome.info
#      utab: data.frame of columns ("gene_name", "gene_id", "transcript_id")
#      geneBiotype: data.frame of columns ("gene_id", "biotype")
#      proteinCodingGenes: vector of gene_id
#      canonicalTranscripts: data.frame of columns ("transcript_id", "gene_id", "CDS.length", "five_prime_utr.length", "three_prime_utr.length")
#      canonicalTranscriptsGeneName: data.frame of columns ("gene_id", "transcript_id", "gene_name", "CDS.length", "five_prime_utr.length", "three_primer_utr.length")
#      proteinCodingCanonTrans: data.frame of columns ("transcript_id", "gene_id", "CDS.length", "five_prime_utr.length", "three_primer_utr.length")
# comment:
# consider using load(sprintf("%s/rdata/jupyter_common_with_gtf.rdata", dir_jupyter))
# see make_jupyter_common.ipynb
# usage:
# list_out <- build_transcriptome_info(species="m38")
build_transcriptome_info <- function(species = species_notebook) {

  myData.species  <-  species
  myData.transcriptome.gtf  <-  load_transcriptome_gtf( species = myData.species )
  myData.transcriptome.gtf  <-  prepare_transcriptome_gtf( gtf = myData.transcriptome.gtf )

  # see functions_transcriptome.R
  #gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "transcript_length" , num = TRUE )
  #gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "five_prime_utr_length" , num = TRUE )
  #gtf  <-  extract_and_append_attribute_field_to_gtf( gtf = gtf , field = "three_prime_utr_length" , num = TRUE )

  myData.transcriptome.info  <-  get_all_transciptome_classification_info( gtf = myData.transcriptome.gtf )
  #myData.transcriptome.feats.fa <- load_all_transcriptome_gtf_feature_fastas(species = myData.species, gtf = myData.transcriptome.gtf)
  #myData.genome.gtf  <-  load_genome_gtf( species = myData.species)

  list(myData.transcriptome.gtf=myData.transcriptome.gtf, myData.transcriptome.info=myData.transcriptome.info)

} # build_transcriptome_info






### gene sets


# load_project_DE_test_tables
# input:
#   file_names: vector
#   project:
#   individual:
# output:
#   data.frame:
#     gene_name	condition1	condition2	expression1	expression2	log2FC	FDR	project	individual
# usage:
# myData.deTest <- load_project_DE_test_tables(c('/home/jupyter/hkim/stjude-nmumg_riboprof/data/2016-11-01--14.45.30.new.mrna.de/emt.161021/emt.161021.limma-voom.mrna.limma.gene.unt48--vs--tgfb48.all.txt','/home/jupyter/hkim/stjude-nmumg_riboprof/data/2016-11-01--14.45.30.new.mrna.de/emt.161021/emt.161021.limma-voom.mrna.limma.gene.tgfb48--vs--tgfbCX5461100nm.all.txt'), project='161021', individual='NMuMG')
#
load_project_DE_test_tables <- function(file_names, project, individual) {

  df_de <- NULL
  for (filename in file_names) {
      df1 <- read.table(file=filename,
                header=TRUE, sep="\t", row.names=1,
                quote="", comment.char="#", stringsAsFactors=F)
      # unt48	tgfb48	log2FC	FDR
      cols <- colnames(df1)
      condition1 <- cols[1]
      condition2 <- cols[2]
      # convert to simpler condition names.
      condition2 <- gsub('tgfbCX5461100nm','tgfbCX5461',condition2)
      colnames(df1) <- c('expression1','expression2','log2FC','FDR','p.value')
      df1$project <- project
      df1$individual <- individual
      df0 <- data.frame(gene_name=rownames(df1), condition1=condition1, condition2=condition2)
      df1 <- cbind(df0,df1)
      rownames(df1) <- NULL
      df_de <- rbind(df_de, df1)
  }

  return(df_de)

}

# build_gene_set
# input:
#   gene_set: 
# output:
#   gene_set$sets: list of list
#     gene_set$sets[['unt48.tgfb48']]$up
#     gene_set$sets[['unt48.tgfb48']]$down
#     gene_set$sets[['unt48.tgfb48']]$notSig
#   gene_set$factors: list of factor: up, down, notSig
#   gene_set$tables: list of data.frame
#     gene_name	condition1	condition2	expression1	expression2	log2FC	FDR	project	individual
#     
# usage:
# gene_set <- build_gene_set(de = myData.deTest, project = "161021", conditions = c("unt48", "tgfb48"), gname="unt48.tgfb48", verbose=T)
# gene_set <- build_gene_set(de = myData.deTest, project = "161021", conditions = c("unt48", "tgfb48", "tgfbCX5461"), gene_set=gene_set, gname="reversible.CX", verbose=T)
#
build_gene_set  <-  function( de , project , conditions , reversible = length(conditions) > 2 , include.notSig = !reversible , FDR.thresh = 0.05 , element = "gene_name" , gene_set=NULL, gname='default', verbose = FALSE ) {

        if (is.null(gene_set)) {
          gene_set <- list()
          gene_set$sets <- list()
          gene_set$factors <- list()
          gene_set$tables <- list()
        } 

        if (is.data.frame(project)) { project  <-  unique(project$project) }

        gset  <-  list()

        if (!reversible) {
                if (verbose) { verb("\tpairwise.\n") }

                stopifnot(length(conditions) == 2)

                ### DE
                gset$up  <-  de[[element]][ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC > 0  &  de$project %in% project ]
                gset$down  <-  de[[element]][ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC < 0  &  de$project %in% project ]

                ### factor and notsig
                if (include.notSig) {
                        gset$notSig  <-  de[[element]][ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR > FDR.thresh  &  de$project %in% project ]
                        gfac  <-  factor(names(gset) , levels = c("notSig","down","up"))
                } else {
                        gfac  <-  factor(names(gset) , levels = c("down","up"))
                } # n.s.

                ### table
                gtab  <-  de[ de$condition1 == conditions[1]  &  de$condition2 == conditions[2]   &  de$project %in% project ,,drop=FALSE]

        } else {
                if (verbose) { verb("\treversible.\n") }

                stopifnot(length(conditions) == 3)
                stopifnot(!include.notSig)
            
                ### rev
                up1  <-  de[[element]][(de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC > 0  &  de$project %in% project)]
                down1  <-  de[[element]][(de$condition1 == conditions[1]  &  de$condition2 == conditions[2]  &  de$FDR <= FDR.thresh  &  de$log2FC < 0  &  de$project %in% project)]

                up2  <-  de[[element]][(de$condition1 == conditions[2]  &  de$condition2 == conditions[3]  &  de$FDR <= FDR.thresh  &  de$log2FC > 0  &  de$project %in% project)]
                down2  <-  de[[element]][(de$condition1 == conditions[2]  &  de$condition2 == conditions[3]  &  de$FDR <= FDR.thresh  &  de$log2FC < 0  &  de$project %in% project)]

                gset$upDown  <-  intersect( up1 , down2 )
                gset$downUp  <-  intersect( down1 , up2 )

                ### factor
                gfac  <-  factor(names(gset) , levels = c("downUp","upDown"))

                ### table
                gtab  <-  de[ ( (de$condition1 == conditions[1]  &  de$condition2 == conditions[2])  | (de$condition1 == conditions[2]  &  de$condition2 == conditions[3]) )  & de$project %in% project ,,drop=FALSE]
        } # rev

        gene_set$sets[[gname]] <- gset
        gene_set$factors[[gname]] <- gfac
        gene_set$tables[[gname]] <- gtab
        return(gene_set)

} # build_gene_set



# compare_gene_set
# input:
# output:
# usage:
# gname  <-  "translationONLY.unt48.tgfb48"
# name1  <-  "riboProf.unt48.tgfb48"
# name2  <-  "polyA.unt48.tgfb48"
# res  <-  compare_gene_set( set1 = gene_set$sets[[name1]] , set2 = gene_set$sets[[name2]] , factor1 = gene_set$factors[[name1]] , table1 = gene_set$tables[[name1]] , table2 = gene_set$tables[[name2]] , type = "complement" )
#
compare_gene_set  <-  function( set1 , set2 , fac , table1 , table2 , type ) {


        common.names  <-  intersect( names(set1) , names(set2) )
        common.names  <-  setdiff( common.names , "notSig" )
        stopifnot(length(common.names) > 1)

        gset  <-  list()
        if (type == "agree") {
                for (namex  in  common.names) {
                        gset[[namex]]  <-  intersect( set1[[namex]] , set2[[namex]] )
                } # namex

        } else if (type == "opposite") {
                stopifnot(length(common.names) == 2)
                for (namex  in  common.names) {
                        name2  <-  setdiff( common.names , namex )
                        gset[[namex]]  <-  intersect( set1[[namex]] , set2[[name2]] )
                } # namex
        } else if (type == "complement") {
                for (namex  in  common.names) {
                        gset[[namex]]  <-  setdiff( set1[[namex]] , set2[[namex]] )
                } # namex
        } else {
                verb("\n\n\nERROR!  unrecognized type=[%s] for comparing gene sets!!!\n", type)
                stop()
        } # type

        ### factor
        gfac  <-  as.character(fac)
        gfac  <-  setdiff( gfac , "notSig" )
        gfac  <-  factor( gfac , levels = setdiff(levels(fac) , "notSig" ) )
        stopifnot(!any(is.na(gfac)))

        ### table
        gtab  <-  rbind( table1 , table2 )

        res  <-  list( set = gset , factor = gfac , table = gtab )

        return(res)

} # compare_gene_set





# build_gene_set_after_comparison
# input:
#   type: {'agree','complement','opposite'}
# output:
#   gene_set$sets
#   gene_set$factors
#   gene_set$tables
# usage:
# gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rna, name1="unt48.tgfb48", name2="rna.unt48.tgfb48", type="complement", gname="unt48.tgfb48.ribo.only")
build_gene_set_after_comparison <- function(gene_set, gene_set2, name1, name2, type, gname) {


  res  <-  compare_gene_set( set1 = gene_set$sets[[name1]] , set2 = gene_set2$sets[[name2]] , fac = gene_set$factors[[name1]] , table1 = gene_set$tables[[name1]] , table2 = gene_set2$tables[[name2]] , type = type )

  gene_set$sets[[gname]] <- res$set
  gene_set$factors[[gname]] <- res$factor
  gene_set$tables[[gname]] <- res$table

  return(gene_set)

}



# build_gene_set_with_rnaseq_riboseq
# input:
#   pattern_rnaseq_gene: (e.g. n-R5s)
#   pattern_riboseq_gene: (e.g. n-R5s)
# output:
#   gene_set: list
#   gene_set$sets: list of list
#     gene_set$sets[['unt48.tgfb48']]$up
#     gene_set$sets[['unt48.tgfb48']]$down
#     gene_set$sets[['unt48.tgfb48']]$notSig
#   gene_set$factors: list of factor: up, down, notSig
#   gene_set$tables: list of data.frame
#     gene_name condition1      condition2      expression1     expression2     log2FC  FDR     project individual
# usage:
# gene_set <- build_gene_set_with_rnaseq_riboseq()
# gene_set <- build_gene_set_with_rnaseq_riboseq(pattern_rnaseq_gene="n-R5s", pattern_riboseq_gene="n-R5s")
build_gene_set_with_rnaseq_riboseq <- function(project_rnaseq="170224", project_riboseq="161021", rundata_appendix=".rdna_rn18s", individual="NMuMG", level="htseq_gene", f_include_tgfb_vs_cx5461=FALSE, pattern_rnaseq_gene=NULL, pattern_riboseq_gene=NULL) {


  list_out <- de_analysis_type_to_directory_and_tag(type = "mrna", project_rnaseq)
  bdir <- list_out$dir
  btag <- list_out$tag

  dir_de <- sprintf("./%s/%s%s/%s", bdir, project_rnaseq, rundata_appendix, individual)
  fname1 <- sprintf("%s/%s.%s.%s.unt48--vs--tgfb48.all.txt.gz", dir_de, project_rnaseq, individual, btag)
  fname2 <- sprintf("%s/%s.%s.%s.tgfb48--vs--tgfbCX5461100nm.all.txt.gz", dir_de, project_rnaseq, individual, btag)
  rnaseq.deTest <- load_project_DE_test_tables(c(fname1, fname2), project = project_rnaseq, individual = individual)

  if (!is.null(pattern_rnaseq_gene)) {
    f <- grepl(pattern_rnaseq_gene, rnaseq.deTest$gene_name)
    rnaseq.deTest <- rnaseq.deTest[f,]
  }

  dir_de <- sprintf("./%s/%s%s/%s", bdir, project_riboseq, rundata_appendix, individual)
  fname1 <- sprintf("%s/%s.%s.%s.unt48--vs--tgfb48.all.txt.gz", dir_de, project_riboseq, individual, btag)
  fname2 <- sprintf("%s/%s.%s.%s.tgfb48--vs--tgfbCX5461100nm.all.txt.gz", dir_de, project_riboseq, individual, btag)
  riboseq.deTest <- load_project_DE_test_tables(c(fname1, fname2), project = project_riboseq, individual = individual)

  if (!is.null(pattern_riboseq_gene)) {
    f <- grepl(pattern_riboseq_gene, riboseq.deTest$gene_name)
    riboseq.deTest <- riboseq.deTest[f,]
  }

  # gene_set_rnaseq
  gene_set_rnaseq <- NULL
  gene_set_rnaseq <- build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("unt48", "tgfb48"), gene_set = gene_set_rnaseq, gname = "unt48.tgfb48.DEtranscription", verbose = T)
  if (f_include_tgfb_vs_cx5461) {
    gene_set_rnaseq <- build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "tgfb48.tgfbCX5461.DEtranscription", verbose = T)
  }
  gene_set_rnaseq <- build_gene_set(de = rnaseq.deTest, project = project_rnaseq, conditions = c("unt48", "tgfb48", "tgfbCX5461"), gene_set = gene_set_rnaseq, gname = "reversible.transcription.CX", verbose = T)


  # gene_set
  gene_set <- NULL
  gene_set <- build_gene_set(de = riboseq.deTest, project = project_riboseq, conditions = c("unt48", "tgfb48"), gene_set = gene_set, gname = "unt48.tgfb48.DEtranslation", verbose = T)
  if (f_include_tgfb_vs_cx5461) {
    gene_set <- build_gene_set(de = riboseq.deTest, project = project_riboseq, conditions = c("tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "tgfb48.tgfbCX5461.DEtranslation", verbose = T)
  }
  gene_set <- build_gene_set(de = riboseq.deTest, project = project_riboseq, conditions = c("unt48", "tgfb48", "tgfbCX5461"), gene_set = gene_set, gname = "reversible.translation.CX", verbose = T)

  # translationONLY
  gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rnaseq, name1="unt48.tgfb48.DEtranslation", name2="unt48.tgfb48.DEtranscription", type="complement", gname="unt48.tgfb48.DEtranslationONLY")
  if (f_include_tgfb_vs_cx5461) {
    gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rnaseq, name1="tgfb48.tgfbCX5461.DEtranslation", name2="tgfb48.tgfbCX5461.DEtranscription", type="complement", gname="tgfb48.tgfbCX5461.DEtranslationONLY")
  }
  gene_set <- build_gene_set_after_comparison(gene_set, gene_set_rnaseq, name1="reversible.translation.CX", name2="reversible.transcription.CX", type="complement", gname="reversible.translationONLY.CX")

  gene_set

} # build_gene_set_with_rnaseq_riboseq





hypergeo_test_all_feat  <-  function( data , feat , sample , good.val , bad.val = NULL ) {
#       data            data.frame with columns feat and group.
#                               The feat col gives a factor, for eac element.
#                                       e.g. the start codon.
#                               The sample col specifies which element are in the drawn sample.
#                                       e.g. "up","down"  would all be in the same sample.
#                                               NA indicates not in the sample.
#       feat            a column in the data df which contains the group labels (a factor).
#       sample          a coulmn in the data df which denotes which elements were drawn
#                               as a sample.
#       good.val        one of the values in the sample col which indicates which elements
#                               were drawn in the sample.
#       bad.val         if non-NULL, then the population size =  num good.val + num bad.val,
#                               otherwise, population size size = nrow(data)
#
#
# feat may be, e.g. start_codon  ("ATG","CGC",etc)
#               or  transcription factor


  ####### check
  if (!all(c(feat,sample) %in% colnames(data))) {
        verb("feat=[%s],sample=[%s] not in data colnames:\n", feat,sample)
        show(colnames(data))
        stop()
  } # check

  if (!(good.val  %in%  data[[sample]])) {
        verb("good.val=[%s] not in sample=[%s]\n", good.val , sample)
        stop()
  } # check



  logi.na.samp  <-  is.na(data[[sample]])
  logi.na.feat  <-  is.na(data[[feat]])


  ### effective population
  if (is.null(bad.val)) {
        logi.pop  <-  !logical(length = nrow(data))
  } else {
        logi.pop  <-  data[[sample]] %in% c(good.val,bad.val)  &  !logi.na.samp
  }

  data  <-  data[ logi.pop ,]
  N.pop.size  <-  nrow(data)


  # make factor
  if (!is.factor(data[[feat]])) {
        data[[feat]]  <-  factor(data[[feat]])
  }

  # non-trivial
  if (length(levels(data[[feat]][!logi.na.feat])) == 0) { return() }

  logi.samp  <-  data[[sample]] == good.val  &  !logi.na.samp
  k.samp.size  <-  sum(logi.samp)

  res  <-  list()

  uniq.levs  <-  unique(as.vector(data[[feat]][!logi.na.feat]))
  if (is.null(uniq.levs)  ||  length(uniq.levs) == 0) { return() }

  for (levx  in  uniq.levs) {
        logi.good  <-  data[[feat]] == levx  &  !is.na(data[[feat]])
        m.pop.good  <-  sum(logi.good)
        n.pop.bad  <-  N.pop.size - m.pop.good
        q.samp.good  <-  sum( logi.good  &  logi.samp )

        p.low  <-  phyper( q = q.samp.good , m = m.pop.good , n = n.pop.bad , k = k.samp.size , lower.tail = TRUE )
        p.upp  <-  phyper( q = q.samp.good , m = m.pop.good , n = n.pop.bad , k = k.samp.size , lower.tail = FALSE )

        expected.perc  <-  m.pop.good / (m.pop.good + n.pop.bad) * 100

        res[[levx]]  <-  data.frame(    level = levx ,
                                        p.depleted = p.low , p.enriched = p.upp ,
                                        num.success.sample = q.samp.good , num.success.pop = m.pop.good ,
                                        sample.size = k.samp.size , pop.size = N.pop.size ,
                                        perc.success.sample = q.samp.good / k.samp.size * 100 ,
                                        perc.success.expected = m.pop.good / N.pop.size * 100 )
        res[[levx]]$log2FC.perc.success  <-  log2(res[[levx]]$perc.success.sample) - log2(res[[levx]]$perc.success.expected)
  } # levx

  res  <-  do.call( rbind , res )
  res$group  <-  good.val

  # smallest p-value
  small.p  <-  pmin( res$p.depleted , res$p.enriched )
  res  <-  res[ order(small.p) , ]

  ### FDR
  if (length(uniq.levs) > 20 ) {
        res$FDR  <-  p.adjust( p = res$small.p , method = "BH" )
        res  <-  res[ order(res$FDR) ,]
  } # FDR

  return(res)

} # hypergeo_test_all_feat





hypergeo_test_all_groups  <-  function( data , feat , sample ) {

        res.bg  <-  list()

        the.samps  <-  as.vector(data[[sample]])
        uniq.groups  <-  unique(the.samps[!is.na(the.samps)])
        for (groupx  in  uniq.groups) {
                res.bg[[groupx]]  <-  hypergeo_test_all_feat( data = data , feat = feat , sample = sample , good.val = groupx )
        } # groupx

        res.bg  <-  do.call( rbind , res.bg )

        return(res.bg)

} # hypergeo_test_all_groups



# load_protein_motif_found
load_protein_motif_found <- function( file.hits , motif_name, fname_fasta ) {


   indf  <-  tryCatch(read.table( file = file.hits , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , fill = TRUE ), error=function(e) NULL)
   if (is.null(indf)) return(indf)

   colnames(indf)  <-  c("transcript_id","motif_start","sequence","name","strand","score")
   indf$motif_end <- indf$motif_start + nchar(indf$sequence) - 1

   idx <- match(indf$transcript_id, myData.transcriptome.info$utab[,"transcript_id"])
   indf$gene_name <- myData.transcriptome.info$utab[idx, "gene_name"]

   # add a column of total length of protein
   aa <- readAAStringSet(fname_fasta)
   idx <- match(indf$transcript_id, names(aa))
   indf[,"seq_len"] <- width(aa)[idx]

   indf <- indf[, c("transcript_id", "sequence", "name", "motif_start", "motif_end", "gene_name", "seq_len")]

}


# read_scanprosite_motif_found
read_scanprosite_motif_found <- function(fname_found, motif_name, fname_fasta) {

             data <- readLines(con <- file(fname_found)); close(con)
             df_tmp <- data.frame(transcript_id=data, sequence="", name=motif_name, motif_start=NA, motif_end=NA, stringsAsFactors=F)
             transcript_id <- ''
             for (i in 1:length(data)) {
                vec <- strsplit(data[i],"[>:]")[[1]]
                vec <- trimws(vec)
                if (grepl("ENS", vec[2])) {
                        transcript_id <- vec[2]
                } else {
                        df_tmp[i,1] <- transcript_id
                        vec <- strsplit(data[i],"[ ]+")[[1]]
                        df_tmp[i,c(2,4,5)] <- vec[c(5,2,4)]
                }
             }
             cudf <- df_tmp[!grepl("^>", df_tmp[,1]),,drop=F]
	     suppressWarnings(	cudf$motif_start <- as.numeric(cudf$motif_start) )
	     suppressWarnings( cudf$motif_end <- as.numeric(cudf$motif_end) )
	     idx <- match(cudf$transcript_id, myData.transcriptome.info$utab[,"transcript_id"])
             cudf$gene_name <- myData.transcriptome.info$utab[idx, "gene_name"]

             # add a column of total length of protein
             aa <- readAAStringSet(fname_fasta)
	     idx <- match(cudf$transcript_id, names(aa))
	     cudf[,"seq_len"] <- width(aa)[idx]

	     cudf
}


# write_fasta_with_boundary
write_fasta_with_boundary <- function( fname_fasta, list_boundary, featx, outpref, alphabet) {

  f_search_codon_pattern_in_cds <- FALSE
  if (grepl("CDS", fname_fasta) && alphabet=="AA") {
    f_search_codon_pattern_in_cds <- TRUE
    alphabet <- "RNA"
    if (any(names(list_boundary) == featx)) {
      # update list_boundary
      s <- list_boundary[[featx]][1]
      e <- list_boundary[[featx]][2]
      if (s < 0) s <- s*3        else s <- (s-1)*3+1
      if (e < 0) e <- (e+1)*3-1  else e <- e*3
      list_boundary[[featx]] <- c(s,e)
    }
  }

  switch(alphabet,
    "DNA"={ ss <- readDNAStringSet(fname_fasta) },
    "RNA"={
       #ss <- readRNAStringSet(fname_fasta)
       ss <- readDNAStringSet(fname_fasta)
    },
    "AA"={ ss <- readAAStringSet(fname_fasta) },
    {}
  )

  if (f_search_codon_pattern_in_cds) {
     mycode <- GENETIC_CODE
     #Selenocysteine  Sec/U   UGA     https://en.wikipedia.org/wiki/Selenocysteine
     mycode[["TGA"]] <- "U"
     #Pyrrolysine     Pyl/O   UAG     https://en.wikipedia.org/wiki/Pyrrolysine
     mycode[["TAG"]] <- "O"

     # trim untranslated tail sequences.
     n_nucleotides_of_tail <- width(ss) %% 3
     ss <- subseq(ss, 1, -1-n_nucleotides_of_tail) 
     # translation
     suppressWarnings( aa  <-  Biostrings::translate( ss, genetic.code=mycode, no.init.codon=FALSE, if.fuzzy.codon="error" ) )

     # remove the last stop codon for interproscan.sh.
     # see make_transcriptome_feature_protein_fasta.R
     f <- as.character(subseq(aa, -1,-1)) %in% c("*","U","O")
     ss[f] <- subseq(ss[f], 1, -4)
     aa[f] <- subseq(aa[f], 1, -2)

     # remove polymorphic pseudogenes that contain stop codons inside aa sequences.
     # see make_transcriptome_feature_protein_fasta.R
     # ENSMUST00000126664.7
     # ENSMUST00000134773.3
     # ENSMUST00000135748.1
     f <- grepl("\\*", as.character(aa))
     ss <- ss[!f]
  }

  if (any(names(list_boundary) == featx)) {
    n <- width(ss)
    s <- list_boundary[[featx]][1]
    e <- list_boundary[[featx]][2]
    if (s < 0) s <- pmax(-n, s) else s <- pmin(s, n)
    if (e < 0) e <- pmax(-n, e) else e <- pmin(e, n)

    f.na <- (sign(s*e) < 0) & (abs(s)+abs(e) > n)
    ss[f.na] <- subseq(ss[f.na], 1, width=0)
    ss[!f.na] <- subseq(ss[!f.na], s[!f.na], e[!f.na])
  }

  # write a temporary fasta file
  fname_out <- sprintf("%s.fa", outpref)
  writeXStringSet(ss, filepath=fname_out, compress=FALSE, format="fasta")

  fname_out

}


# performt_tests
perform_tests <- function(input, feats.to.test, input_type, gene_set, f_protein_coding_only=TRUE, df_gene_info=NULL, max_gset=Inf, dir_out, fname_prefix='', list_boundary=list("protein_nterminus"=c(1,200), "protein_middle"=c(201,-201), "protein_cterminus"=c(-200,-1), "CDS_nterminus"=c(1,200), "CDS_middle"=c(201,-201), "CDS_cterminus"=c(-200,-1)), verbose=TRUE) {

  list_cudf <- list()
  alltestdf  <-  data.frame(); n_feat <- 0
  for (featx  in  feats.to.test) {
        if (verbose) verb("%s\n", featx)

        outLocalDir <- sprintf("%s/%s" , dir_out, featx )
        dir.create(outLocalDir, recursive = TRUE , showWarnings = FALSE)
        outLocalPrefix  <-  sprintf("%s/%s.%s", outLocalDir, fname_prefix, featx)

        ### sumdf, pcdf
        n_feat <- n_feat+1
        switch(input_type,

          # Ybx1
          "rna-binding_protein_mrna_motif"={
	     # /home/hkim5/riboprof/out/transcriptome-gtf/m38/Mus_musculus.GRCm38.97.rdna_rn18s.transcriptome.five_prime_utr.fa.gz
	     # >ENSMUST00000000001
             tran.fname  <-  get_transcriptome_feature_fasta_filename( species = myData.species , feat = featx )
             fname_fasta <- write_fasta_with_boundary( fname_fasta=tran.fname, list_boundary, featx, outLocalPrefix, alphabet="RNA")
             if (verbose) verb("\tfind_motifs_in_fasta: %s\n", tran.fname)
             Sys.setenv(HOMER="/home/hkim5/packages/homer_v4.9/bin")
	     find_motif_in_fasta( motif=input[['motif']], name=input[['motif_name']], fasta = fname_fasta , outpref = outLocalPrefix  , mismatch = 0 , verbose = verbose )
             if (verbose) verb("\tload_homer_motif_annotation\n")
             cudf  <-  load_homer_motif_annotation( hits = sprintf("%s.found" , outLocalPrefix ) , bed = TRUE )
             if (is.null(cudf)) {
                system(sprintf("rm -f %s", fname_fasta))
		next
	     }
             cudf  <-  subset( cudf , strand == "+" )
	     suppressMessages(suppressWarnings( sumdf  <-  cudf  %>%  group_by(chr,name)  %>%  summarise(num.hits = n()) ))
             sumdf  <-  as.data.frame(sumdf)
	     if (f_protein_coding_only) {
               pcdf <- merge(myData.transcriptome.info$proteinCodingCanonTrans, sumdf, by.x = "transcript_id", by.y = "chr", all.x = TRUE)
	     } else {
               pcdf <- merge(myData.transcriptome.info$canonicalTranscripts, sumdf, by.x = "transcript_id", by.y = "chr", all.x = TRUE)
	     }
             pcdf$num.hits[is.na(pcdf$num.hits)] <- 0
             # hasCUbox, noCUbox
             #pcdf$feature  <-  ifelse(pcdf$num.hits > 0 , input[['motif_name']] , paste0("no_",input[['motif_name']]))
             pcdf$feature  <-  ifelse( pcdf$num.hits > 0 , featx , paste0("no_",featx))
             pcdf <- merge(pcdf, myData.transcriptome.info$canonicalTranscriptsGeneName[, c("transcript_id", "gene_name"), drop = FALSE], by = "transcript_id", all.x = TRUE)
             # pcdf: data.frame
             # gene_name	gene_id	transcript_id	CDS.length	five_prime_utr.length	three_prime_utr.length	name	num.hits	feature
             stopifnot(nrow(pcdf) == length(unique(pcdf$transcript_id)) )
          },

          # Larp1
          "table_mrna_5utr_cds_3utr"={
             if (verbose) verb("\tload %s.\n", input[n_feat])
             df <- read.table(input[n_feat], header=T, row.names=1, sep='\t', quote='')
             strcol <- make.names(tail(strsplit(featx, '_')[[1]],1))
             idx <- which(df[,strcol]==1)
             sym_human <- df[idx,'Gene.Name']
             sym_human <- sym_human[!is.na(sym_human)]
             sym_mouse <- get_orthologous_gene_names( from = "human" , to = "mouse" , conv = orthodf , genes = sym_human )
             cudf <- data.frame(gene_name=sym_mouse, name=featx, num.hits=1)
	     sumdf <- cudf
             pcdf  <-  merge( myData.transcriptome.info$canonicalTranscriptsGeneName , sumdf , by = "gene_name", all.x = TRUE)
             pcdf$num.hits[is.na(pcdf$num.hits)]  <-  0
             # mTOR_active_Larp1_5'UTR, no_mTOR_active_Larp1_5'UTR
             pcdf$feature  <-  ifelse( pcdf$num.hits > 0 , featx , paste0("no_",featx))
             # pcdf: data.frame
             # gene_name	gene_id	transcript_id	CDS.length	five_prime_utr.length	three_prime_utr.length	name	num.hits	feature
             stopifnot(nrow(pcdf) == length(unique(pcdf$transcript_id)) )
          },

          "table_gene_1st_column"={
             if (verbose) verb("\tload %s.\n", input[n_feat]) 
             sumdf <- read.table(input[n_feat], header=T, row.names=1, sep='\t', quote='')
             sym_human <- sumdf[,1]
             sym_mouse <- get_orthologous_gene_names( from = "human" , to = "mouse" , conv = orthodf , genes = sym_human )
             cudf <- data.frame(gene_name=sym_mouse, name=featx, num.hits=1)
             sumdf <- cudf
             pcdf  <-  merge( myData.transcriptome.info$canonicalTranscriptsGeneName , sumdf , by = "gene_name", all.x = TRUE)
             pcdf$num.hits[is.na(pcdf$num.hits)]  <-  0
             # mTOR_inactive_TR, no_mTOR_inactive_TR
             pcdf$feature  <-  ifelse( pcdf$num.hits > 0 , featx , paste0("no_",featx))
             # pcdf: data.frame
             # gene_name	gene_id	transcript_id	CDS.length	five_prime_utr.length	three_prime_utr.length	name	num.hits	feature
             stopifnot(nrow(pcdf) == length(unique(pcdf$transcript_id)) )
          },
          {}
       )

	# store pcdf
	list_cudf[[featx]] <- cudf

        # gene_set
        switch(class(gene_set),
          "enrichResult"={
             df_enricher <- as.data.frame(gene_set)
             df_enricher <- df_enricher[order(df_enricher$p.adjust),]
             gene_set_names <- gsub(' ','_',df_enricher$Description)
          },
          { gene_set_names <- names(gene_set$sets) }
        ) 

        alldf  <-  data.frame()

        n_gset <- 0
        for (gsetx  in  gene_set_names) {
                n_gset <- n_gset + 1
                if (n_gset > max_gset) break
                if (verbose) verb("\t\t%s\n", gsetx)

                outLocalDir <-  sprintf("%s/%s/%s", dir_out , featx, gsetx)
                dir.create(outLocalDir ,  recursive = TRUE , showWarnings = FALSE)
                outLocalPrefix  <-  sprintf("%s/%s.%s.%s" , outLocalDir, fname_prefix, featx, gsetx )

                # gtab: define universe (protein coding genes detected in a platform)
                # mdf: define genes selected, colnames=c('gene_name','direc')
                switch(class(gene_set),
                   "enrichResult"={
                      if (is.null(df_gene_info)) {
                         gtab <- data.frame(gene_name=pcdf$gene_name, biotype='protein_coding')
                      } else {
                         gtab <- df_gene_info
                      }
                      syms <- strsplit(df_enricher[n_gset,'geneID'],'/')[[1]]
                      mdf <- data.frame(gene_name=syms, direc=df_enricher[n_gset,'Description'])
                   },
                   {
                      gtab  <-  gene_set$tables[[gsetx]]
                      mdf  <-  reshape2::melt( gene_set$sets[[gsetx]] )
                   }
                )
                gtab  <-  merge( gtab , pcdf , by = "gene_name" , all.y = TRUE )
                # begin addition by H. Kim
                f <- !is.na(gtab$gene_name)
                gtab <- gtab[f, , drop = F]
                # end addition
                colnames(mdf)  <-  c("gene_name","direc")
                mdf  <-  mdf[ !is.na(mdf$direc)  &  mdf$direc != "notSig" , ,drop=FALSE]

                # protein coding, with 3' UTR
                mdf  <-  subset( mdf , gene_name  %in%  gtab$gene_name )
                groups  <-  mdf
                colnames(groups)  <-  c("element","group")

                # begin modification by H. Kim
                #tdf  <-  test_hypergeometric_enrichment( universe = gtab , groups = groups , prop = "has.CUbox" , element = "gene_name" , strict = TRUE )
                classes <- gtab[, c("gene_name", "feature")]
                colnames(classes) <- c("element", "group")            
                tdf <- test_hypergeometric_enrichment(universe = gtab$gene_name, classes = classes, groups = groups)            
                # end modification
                tdf$geneset  <-  gsetx
                tdf$feature  <-  featx
                alldf  <-  rbind( alldf , tdf )

                annotdf  <-  merge( gtab ,  mdf , by = "gene_name" , all = TRUE )
                fname  <-  sprintf("%s.annotation.txt" , outLocalPrefix)
                write.table( annotdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

        } # gsetx


        alldf  <-  subset( alldf , class == featx )

        # begin modification by H. Kim
        #alldf  <-  alldf[ order(alldf$p.over.represented) , ,drop=FALSE]
        alldf  <-  alldf[ order(alldf$p.overrepresented) , ,drop=FALSE]
        # end modification 
        alltestdf  <-  rbind( alltestdf , alldf )

        # save
        outLocalDir <-  sprintf("%s/%s", dir_out, featx )
        dir.create(outLocalDir ,  recursive = TRUE , showWarnings = FALSE)
        outLocalPrefix <- sprintf("%s/%s.%s", outLocalDir, fname_prefix, featx)

        fname  <-  sprintf("%s.tests.txt", outLocalPrefix)
        write.table( alldf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

        if (exists("fname_fasta")) {
          # delete temporaty files
          system(sprintf("rm -f %s", fname_fasta))
        }

  } # featx

  ### save
  outLocalDir <-  dir_out
  dir.create(outLocalDir ,  recursive = TRUE , showWarnings = FALSE)
  outLocalPrefix  <-  sprintf("%s/%s", outLocalDir, fname_prefix)

  fname  <-  sprintf("%s.tests.txt", outLocalPrefix)
  write.table( alltestdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
  # begin of addition by H. Kim
  if (nrow(alltestdf) == 0 || (!"geneset" %in% colnames(alltestdf))) {
	 return(list(tdf=NULL, tmat=NULL, list_cudf=list_cudf))
  }
  # end of addition

  ### plot
  tdf  <-  alltestdf
  tdf$minp  <-  pmin(tdf$p.underrepresented , tdf$p.overrepresented)
  tdf$log10p  <-  log10(tdf$minp)
  # begin modification by H. Kim
  #tdf  <-  tdf[ tdf$minp <= 1e-2  &&  tdf$class == "hasCUbox" ,, drop=FALSE]
  #tdf  <-  tdf[ tdf$minp <= 1e-2 ,, drop=FALSE]
  # include all lines instead of selecting enriched lines for drawing
  #tdf$direc  <-  ifelse( tdf$p.underrepresented < 1e-2 , "underrepresented" , "overrepresented" )
  #tdf$dlog10p  <-  ifelse( tdf$p.overrepresented < 1e-2 , abs(tdf$log10p) , tdf$log10p )
  tdf$direc  <-  ifelse( tdf$p.overrepresented > tdf$p.underrepresented, "underrepresented" , "overrepresented" )
  tdf$dlog10p  <-  ifelse( tdf$p.overrepresented < tdf$p.underrepresented, abs(tdf$log10p), tdf$log10p )
  # when success.population is too small, set nonSig 
  f <- tdf$success.population < 10
  tdf[f, "dlog10p"] <- 0
  # end modification 

  tmat  <-  acast( data = tdf , formula = geneset + group ~ feature , value.var = "dlog10p" , fill = 0 )
  # begin of addition by H. Kim
  rownames(tmat) <- gsub('_', ' ', rownames(tmat))
  f <- feats.to.test %in% colnames(tmat)
  feats.to.test <- feats.to.test[f]
  # end of addition
  tmat  <-  tmat[, feats.to.test,drop=FALSE]
  tmat  <-  tmat[ !str_detect(rownames(tmat) , "no_") ,, drop=FALSE]
  f <- tdf$group == gsub('_',' ',tdf$geneset)
  if (all(f)) {
    rownames(tmat) <- sapply(rownames(tmat), function(x) {
      vec <- str_split(x,'_')[[1]]
      tail(vec,1) } )
  }

  thresh  <-  15
  tmat[tmat >  thresh]  <-  thresh
  tmat[tmat < -thresh]  <-  -thresh

  # save
  fname  <-  sprintf("%s.tests.plot.tdf.txt", outLocalPrefix)
  write.table( tdf , file = fname , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
  fname  <-  sprintf("%s.tests.plot.tmat.txt", outLocalPrefix)
  write.table( tmat , file = fname , quote = FALSE , sep = "\t" , row.names = TRUE , col.names = TRUE )

  return(list(tdf=tdf, tmat=tmat, list_cudf=list_cudf))


} # performt_tests





# process_orf_resdf
# input:
#   resdf: generated by test_ORF_classification_read_count_between_conditions()
#   gene_set: 
#   good.comps: vector
#   classx: {"uORF", ... }
# ouput:
#   list$tmat
# usage:
# resdf <- test_ORF_classification_read_count_between_conditions( gtf = NULL , rpf = NULL , species = "m38" , groups = groups, project = "161021" , classes = NULL , outfbase = NULL , verbose = FALSE )
# list_tests <- process_orf_resdf(resdf, gene_set, good.comps=c("unt48 tgfb48" , "tgfb48 tgfbCX5461100nm"), classx="uORF")
process_orf_resdf <- function(resdf, gene_set, good.comps, classx) {

  gsdf <- reshape2::melt(gene_set$sets)
  colnames(gsdf) <- c("gene_name", "direc", "geneset")
  gsdf$group <- paste(gsdf$geneset, gsdf$direc, sep = " ")
  gsdf <- gsdf[gsdf$direc != "notSig", , drop = FALSE]
  groups <- gsdf

  rdf  <-  merge( resdf , unique(gsdf[,c("geneset","group","direc"),drop=FALSE]) , by = "group" , all.x = TRUE )
  rdf$comparison  <-  paste( rdf$condition1 , rdf$condition2 , sep = " " )

  rdf  <-  rdf[rdf$comparison  %in% good.comps ,,drop=FALSE]
  rdf$dlog10p  <-  ifelse( rdf$log2FC > 0 , -log10(rdf$p.value) , log10(rdf$p.value) )
  rdf$dlog10p[rdf$p.value >= 0.05]  <-  0
  rdf  <-  rdf[order(rdf$geneset,rdf$direc),,drop=FALSE]


  subdf <- rdf[rdf$classification == classx & !is.na(rdf$classification), , drop = FALSE]
  amat <- acast(data = subdf, formula = group ~ comparison, value.var = "dlog10p", fill = 0)
  if (nrow(amat) == 0 || ncol(amat) == 0) {
    return(NULL)
  }
  amat <- amat[intersect(unique(rdf$group), rownames(amat)), good.comps, drop = FALSE]

  thresh <- 10
  tmat <- amat
  tmat[tmat > thresh] <- thresh
  tmat[tmat < -thresh] <- -thresh

  heat.max.break <- max(abs(tmat)) + 1
  heat.max.break <- max(heat.max.break, 4)
  heat.by.break <- 2 * heat.max.break/100

  list(tmat=tmat)

} # process_orf_resdf




# verb
verb <- function(...) cat(sprintf(...), sep='', file=stdout())



