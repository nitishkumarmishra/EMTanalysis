



##########  packages
library(methods)
library(stats)











#' Test enrichment and depletion via hypergeometric distribution
#'
#' @param  universe		character vector of element names for the entire universe
#' @param  classes		Data frame. Master table containing all gene <-> class annotations.
#'					Must have columns "element","group", where element is e.g. gene name,
#'					and where group is e.g. the regulating transcript factor, etc.
#'					The universe may have the same element in multiple different classifications 
#'					(e.g. the same gene regulated by multiple different transcriptiion factors.)
#'					It may also have group = NA. These will not be tested as a class.
#' @param  groups		data frame: must have columns "element" and "group". element is e.g. gene name, group is e.g. "up", "down", etc.
#'					The same element may appear in multiple different groups.
#' @return dataframe of p-value info
#'
#' @export
test_hypergeometric_enrichment  <-  function( universe , classes , groups , verbose = FALSE ) {

	##### check
	if (verbose) { verb("\tchecking.\n") }
	
	stopifnot(is.data.frame(groups))
	stopifnot(is.data.frame(classes))
	stopifnot(all(complete.cases(groups)))
	stopifnot(all(complete.cases(classes)))
	stopifnot(all(c("element","group") %in% colnames(classes)))
	stopifnot(all(c("element","group") %in% colnames(groups)))
	stopifnot(length(setdiff(groups$element , universe)) == 0)
	stopifnot(length(setdiff(classes$element , universe)) == 0)


	##### prepare
	if (verbose) { verb("\tpreparing.\n") }

	# unique
	classes  <-  unique(classes[,c("group","element"),drop=FALSE])
	groups  <-  unique(groups[,c("group","element"),drop=FALSE])

	# prepare
	pop.size  <-  length(unique(universe))

	classes$group  <-  factor( as.character(classes$group) , levels = unique(as.character(classes$group)) )

	ctab  <-  table(classes$group)
	classdf  <-  data.frame( group = names(ctab) , pop.success = as.numeric(ctab) )


	##### check again
	if (verbose) { verb("\tchecking again.\n") }

	stopifnot(!any(is.na(groups)))
	stopifnot(!any(is.na(classes)))
	stopifnot(!any(is.na(classdf)))


	resdf  <-  data.frame()
	for (groupx  in  unique(groups$group)) {
		if (verbose) { verb("\t\t%s\n" , groupx ) }

		genes  <-  groups$element[groups$group == groupx]

		subclass  <-  classes[ classes$element %in% genes  ,,drop=FALSE]
		stab  <-  table(subclass$group)

		scdf  <-  data.frame( group = names(stab) , sample.success = as.numeric(stab) )
		overdf  <-  merge( classdf , scdf , by = "group" , all = TRUE)
		stopifnot(!any(is.na(overdf)))

		q  <-  overdf$sample.success
		m  <-  overdf$pop.success
		n  <-  pop.size - overdf$pop.success
		k  <-  length(genes)

		# Recall that:
		#	lower.tail = TRUE  -->  P[X <= x]
		#	lower.tail = FALSE -->  P[X  > x]
		
		utail  <-  phyper( q = q , m = m , n = n , k = k , lower.tail = FALSE , log.p = FALSE )
		dens.x  <-  dhyper( x = q , m = m , n = n , k = k , log = FALSE )
		p.over  <-  utail + dens.x
		p.under  <-  1.0 - utail
		
		rdf  <-  data.frame( class = overdf$group , group = groupx ,
					success.sample = q , sample.size = k , success.population = m , population.size = pop.size ,
					frac.success.sample = q / k , frac.success.population = m / pop.size ,
					p.overrepresented = p.over , p.underrepresented = p.under )
		resdf  <-  rbind( resdf , rdf )
	} # groupx
	
	resdf  <-  resdf[order(resdf$group , resdf$p.overrepresented , resdf$success.population , resdf$sample.size) ,,drop=FALSE]
	
	return(resdf)
					
} # test_hypergeometric_enrichment







