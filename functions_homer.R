



##########  packages

library(plyr)
library(dplyr)
library(reshape2)
library(methods)
library(stringr)
library(Biostrings)

source("./functions_enrichment_test.R")






######################################### verb

### use this function for stderr output
verb <- function(...) cat(sprintf(...), sep='', file=stderr())














######################################################
######################################################
######################################################
######################################################
########################  load #######################
######################################################
######################################################
######################################################
######################################################





# This loads the file resulting from "findMotifs.pl" with the "-find" option.
#' Load motif hits from searching for a motif within target sequence.
#'
#' After searching for motifs in a target fasta, e.g. with "find_motif_in_fasta",
#' this function will load the results.
#' Note that the "sequence" given by a HOMER hit always corresponds to the forward strand, even if the motif was found on the reverse strand.
#'
#' @param  hits		filename giving the results from HOMER motif search, e.g. from function "find_motif_in_fasta"
#' @param  bed		if TRUE, will convert the results to BED format (no loss of information).
#'
#' @return  if bed = TRUE, then will return a dataframe in BED format.
#'		if bed = FALSE, will return a dataframe in HOMER results format.
load_homer_motif_annotation  <-  function( hits , bed = FALSE ) {

	# begin of modification by H. Kim
	#indf  <-  read.table( file = hits , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , fill = TRUE )
	indf  <-  tryCatch( read.table( file = hits , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , fill = TRUE ), error=function(e) NULL )
        if (is.null(indf)) return(indf)
        # end of modification

	colnames(indf)  <-  c("chr","offset","sequence","motif","strand","score")

	if (bed) {
		indf  <-  homer_motif_annotation_to_bed( annot = indf )
	} # bed

	return(indf)

} # load_homer_motif_annotation








# Internal
# This loads the file created by "seq2profile.pl" function
load_homer_motif  <-  function( motif ) {

	header  <-  readLines( con = motif , n = 1 )

	### matrix
	inmat  <-  read.table( file = motif , header = FALSE , row.names = NULL , sep="\t" , stringsAsFactors = FALSE , quote = "\"" , skip = 1 )
	inmat  <-  as.matrix(inmat)

	inmotif  <-  list( header = header , matrix = inmat )

	return(inmotif)

} # load_homer_motif









# INTERNAL
get_acceptable_iupac_letters  <-  function() {
	let  <-  c("A","C","G","T","R","Y","S","W","K","M","B","D","H","V","N")
	return(let)
} # get_acceptable_iupac_letters




















######################################################
######################################################
######################################################
######################################################
#################  motif annotate  ###################
######################################################
######################################################
######################################################
######################################################



#' Interface to find a known motif within a set of sequences using HOMER
#'
#' Given a known motif and a set of sequences, this script will use the HOMER package to 
#' annotate/identify instances of the motif within the target fasta.
#'
#' @param  motif	character string giving the motif, in IUPAC format (use "T" not "U").  e.g. "SCTBYC" (CU box motif)
#' @param  name		character string giving the name of this motif.
#' @param  fasta	Biostrings DNAStringSet object or filename of a fasta file containing target sequences in which the motif will be searched for.
#' @param  outpref	character string giving output filename prefix.
#' @param  mismatch	integer giving number of allowable mismatches when searching for this motif.
#' @param  discard.target	remove target fasta temporary file when done.
#'
#' @return   Nothing. This function formulates the command line call and uses a system call to execute it. The output is written to a file.
#'		To read the results, use the function "load_homer_motif_annotation"
#'
#' @export
find_motif_in_fasta  <-  function( motif , name , fasta , outpref , mismatch = 0 , discard.target = TRUE , verbose = TRUE ) {

	### make Homer motif file
	if (verbose) verb("\tmake Homer motif file.\n")

        # begin of addition by H. Kim
        if (grepl('\\.', motif)) {
          mo.fname <- motif
	  inmotif  <-  load_homer_motif( motif = mo.fname )
        } else {
        # end of addition
	  # check
	  iupac  <-  get_acceptable_iupac_letters()
	  stopifnot(all(sapply( str_split(motif , "" ) , FUN=function(x)  x %in% iupac )))

	  mo.fname  <-  sprintf("%s.motif" , outpref )

	  cmd  <-  sprintf("$HOMER/seq2profile.pl  %s  %d  %s > %s" , motif , mismatch , name , mo.fname )
	  system(cmd)

	  # check
	  inmotif  <-  load_homer_motif( motif = mo.fname )
	  stopifnot(nrow(inmotif$matrix) == nchar(motif) )
        # begin of addition by H. Kim
        }
	# end of addition


	local.fa.fname  <-  sprintf("%s.target.fa" , outpref )

	if (class(fasta) == "DNAStringSet") {
		writeXStringSet( fasta , local.fa.fname )
	} else if (is.character(fasta)  &&  length(fasta) == 1) {
		system(sprintf("zcat  -f  %s  >  %s" , fasta , local.fa.fname))
	} else {
		verb("\n\n\nunrecognized fasta!!!!\n")
		show(class(fasta))
		stop()
	} # fasta


	### find motifs
	if (verbose) verb("\tfind motifs in fasta.\n")

	ffound  <-  sprintf("%s.found" , outpref)
	#find.dir  <-  sprintf("%s.find" , outpref )
	#cmd  <-  sprintf("$HOMER/findMotifs.pl  %s  fasta  %s  -fasta %s  -find %s  >  %s.found.bed" , fasta , find.dir , fasta , mo.fname , outpref )
	cmd  <-  sprintf("$HOMER/homer2  find  -i %s  -m %s  -o %s  -offset 1 " , local.fa.fname , mo.fname , ffound )
	# with offset = 1, then if the motif starts at the first position in the FASTA and is on the forward strand, then in the results file it will have position = 1

	system(cmd)


	if (discard.target) {
		file.remove(local.fa.fname)
	} # discard.target
	
} # find_motif_in_fasta








# Internal
homer_motif_annotation_to_bed  <-  function( annot ) { 

	annot$start  <-  ifelse( annot$strand == "+" , annot$offset - 1 , annot$offset - nchar(annot$sequence) )
	annot$end  <-  annot$start + nchar(annot$sequence)

	abed  <-  data.frame( chr = annot$chr , start = annot$start , end = annot$end , name = annot$motif , score = annot$score , strand = annot$strand , sequence = annot$sequence )
	return(abed)

} # homer_motif_annotation_to_bed
















######################################################
######################################################
######################################################
######################################################
#######################  test ########################
######################################################
######################################################
######################################################
######################################################



#' Find a given motif in transcripts and test for enrichment in gene sets
#'
#' This function will take a known motif and find its occurrences among sequences in a fasta file.
#' It will then test if transcripts containing the given motif are over-represented in each gene set.
#'
#' @param  motif	character string giving motif in IUPAC coding.
#' @param  name		name of motif
#' @param  fasta	Biostrings DNAStringSet or filename of fasta (can be .gz) to be tested (e.g. transcript 5' UTRs). Names of sequences should be transcript_ids.
#' @param  outpref	character string giving output filename prefix
#' @param  sets		data frame containing gene sets of interest.  column names must include "group" and "gene_name" or "transcript_id" or "gene_id"
#' @param  universe     character vector defining the universe. THe values of the character vector should correspond to the "level" argument.
#'                      i.e. if level == "gene_name", then "universe" should consist of gene names.
#'                      If NULL, the rownames of "quant" will be used.
#'                      If subset is not NULL, then the universe will be subsetted according to 'subset' argument.
#' @param  trantable	result from function "get_gene_transcript_table"
#' @param  canon                list returned by "get_canonical_transcript_tables"
#' @param  subset               If subset == "protein_coding", then the universe will be intersected with canonical transcripts for protein coding genes.
#'                              If subset == "non_protein_coding", then the universe will be intersected with canonical transcripts for non-protein coding genes.
#'                              If subset is NULL, then the universe will not be filtered by protein-coding or non-protein coding genes, and instead
#' @param  level	"gene_name" or "gene_id" or "transcript_id", specifying which kind of IDs are provided in 'universe' and 'sets'.
#' @param  gtf          transcriptome gtf, as obtained from "load_transcriptome_gtf"
#' @param  species	character string giving the species being analyzed.  Only necessary if quant is NULL and gtf, info, or fasta is NULL.
#'
#' @export
test_motif_enrichment_in_transcripts  <-  function( motif , name , fasta , outpref , sets , universe , trantable = NULL , canon = NULL , subset = "protein_coding" , level , gtf = NULL , species = NULL , verbose = FALSE ) { 

	stopifnot(level %in% c("gene_name","gene_id"))
	stopifnot(subset %in% c("protein_coding" , "non_protein_coding") | is.null(subset))


	############# find motifs in fasta
	if (verbose) { verb("\t\tfind motifs in fasta.\n") }

	find_motif_in_fasta( motif = motif , name = name , fasta = fasta , outpref = outpref , discard.target = FALSE )

	targfname  <-  sprintf("%s.target.fa", outpref)
	targfa  <-  readDNAStringSet(targfname)
	tnames  <-  names(targfa)

	# remove
	file.remove(targfname)
	remove(targfa)



	############### load motif
	if (verbose) { verb("\t\t\tload motif.\n") }

	hits.fname  <-  sprintf("%s.found", outpref)
	modf  <-  load_homer_motif_annotation( hits = hits.fname , bed = TRUE )



	############# process motif hits
	if (verbose) { verb("\t\tprocess motif hits.\n") }

	### must be forward strand for transcripts
	modf  <-  modf[ modf$strand == "+" ,,drop=FALSE]

	### num hits
	modf  <-  modf  %>%  group_by(chr,name)  %>%  mutate(num.hits = n())
	modf  <-  as.data.frame(modf)
	modf$has.motif  <-  modf$num.hits > 0




	

	#### gtf
	if (is.null(gtf)) {
		if (verbose) { verb("\t\tload gtf.\n") }
		gtf  <-  load_transcriptome_gtf( species = species )
	} # proteinCodingTranscripts
	gtf  <-  prepare_transcriptome_gtf( gtf = gtf )


	# transcript and gene name tables.
	if (verbose) { verb("\t\ttables.\n") }

	if (is.null(trantable)) {
		trantable  <-  get_gene_transcript_table( gtf = gtf , verbose = verbose )
	} # trantable

	if (is.null(canon)) {
		if (verbose) { verb("\t\tcanon.\n") }
		canon  <-  get_canonical_transcript_tables( gtf = gtf , verbose = verbose )
	}


	### universe
	# subset and universe to transcript id
	if (verbose) { verb("\t\tinferring universe.\n") }

	universe  <-  subset_and_get_canonical_transcript( universe = universe , subset = subset , level = level , canon = canon , verbose = verbose )
	show(head(universe))
	show(length(universe))
	stopifnot(!any(is.na(universe)))


	### restrict universe
	if (verbose) { verb("\t\trestrict universe.\n") }

	universe  <-  intersect( universe , tnames )


	### convert gene set to transcript id and restrict
	if (verbose) { verb("\t\trestrict sets.\n") }

	sets  <-  convert_sets_to_transcript_id_and_restrict_to_universe( sets = sets , level = level , universe = universe , transcripts = trantable , verbose = verbose )


show(head(sets))

        ############ hypergeo
        if (verbose) { verb("\t\tmake classes.\n") }

	hasdf  <-  modf[ modf$num.hits > 0 ,,drop=FALSE]
	classes  <-  unique(data.frame( element = hasdf$chr , group = hasdf$name ))
	classes  <-  classes[ classes$element %in% universe ,,drop=FALSE]
	

	### test
        if (verbose) { verb("\t\thypergeo.\n") }

	groups  <-  sets
	groups$element  <-  sets$transcript_id
	geodf  <-  test_hypergeometric_enrichment( universe = universe , classes = classes , groups = groups , verbose = verbose )



	######### annotate
	if (verbose) { verb("\t\tannotate.\n") }

	andf  <-  hasdf
	andf  <-  merge( sets , modf , by.x = "transcript_id" , by.y = "chr" )
	

	### return
	resl  <-  list(enrich = geodf , annotate = andf )
	return(resl)

} # test_motif_enrichment






























