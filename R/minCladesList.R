#' @title minimum length set of clades with tip labels
#' 
#' Returns a smallest set of non-intersecting, non-sister clades of an input tree, and which together contain all of the tip labels supplied in the vector tiplabs.
#' REQUIRES ggtree
#' 
#' @param tree A tree object of class phylo or treedata.
#' @param tiplabs A character vector with a subset of the tip labels of the input tree
#' @param quiet Logical indicating whether to suppress some messages printed to screen. Default TRUE. This does not affect warnings.
#' @param warn Logical indicating if a warning should be generated when one or more tiplabs do not occur in the input tree. Default TRUE.
#' @return A list of character vectors, each containing a subset of the names supplied to tiplabs, and each characterizing a clade that is non-intersecting and non-sister to any other clade in the list.
#' @export minCladesList
minCladesList <- function(tree,tiplabs,quiet=T,warn=T) {
	## Assigning individuals to independent clades
	tib     <- tibble::as_tibble(tree)
	numtips <- length(tree$tip.label)
	if(!all(tiplabs %in% unlist(tib[,'label']))){
		if(warn){
			missingLabs <- tiplabs[grepl(F,tiplabs %in% unlist(tib[,'label']))]
			warning(sprintf("%s tiplabs not in tree: %s ...",length(missingLabs),paste(head(missingLabs),collapse=',')))
		}
	}
	samplesRemaining <- intersect(tiplabs, unlist(tib[,'label']))
	if(!quiet){
		print(sprintf("Assigning %s tips with shared factor to fewest number of clades",length(samplesRemaining)))
	}
	# name a list
	cladelist <- list()
	# initial sample node
	node.temp <- c(na.omit(unlist(tib[tib[,'label']==samplesRemaining[1],'node'],use.names=F)))
	while (length(samplesRemaining) > 0 ){
		# parent node of node.temp
		parent.temp <- c(na.omit(unlist(tib[tib[,'node']==node.temp,'parent'],use.names=F)))
		# all descendents of parent
		grp     <- tibble::as_tibble(tidytree::groupClade(tree,parent.temp))
		members <- c(na.omit(unlist(grp[c(grp[,'group']==1),'label'],use.names=F)))
		# if all members are in samplesRemaining, set node.temp to parent.temp
		if(all(members %in% samplesRemaining)){
			node.temp <- parent.temp
		} else {
			if(node.temp<=numtips){
				members.off <- tib$label[unlist(tib$node==node.temp)]
			} else {
				grp.off     <- tibble::as_tibble(tidytree::groupClade(tree,node.temp))
				members.off <- c(na.omit(unlist(grp.off[c(grp.off[,'group']==1),'label'],use.names=F)))
			}
			clade.temp  <- list(members.off)
			names(clade.temp) <- paste0("clade",node.temp)
			cladelist   <- c(cladelist,clade.temp)
			samplesRemaining <- setdiff(samplesRemaining,members.off)
			node.temp <- c(na.omit(unlist(tib[tib[,'label']==samplesRemaining[1],'node'],use.names=F)))
			if(!quiet){
				print(sprintf("Clade with %s tips added to list of clades. %s unassigned tips remaining.",length(members.off),length(samplesRemaining)))
			}
		}
	}
	cladelist
}
