#' @title substitute a clade with another clade
#' 
#' replaces a specified crown clade in a tree with a specified crown clade from a second tree.
#' 
#' @param tree1 phylo object containing the clade that will be replaced
#' @param tree2 phylo object containing the replacement clade; if tips2 is NULL, then tree2 must be rooted
#' @param tips1 vector with tip names in tree1 with MRCA defining the clade to be replaced
#' @param tips2 NULL (default) or a vector with tip names in tree2 with MRCA defining the replacement clade; if NULL, then tree2 is the replacement clade.
#' @param use.stem One of "tree1", "tree2", or "diff" (Default "tree1"). Determines the stem branch length of the substitute clade on the output tree. Setting to "tree1" or "tree2" uses the stem branch length of the clade to be replaced on tree1 or the replacement clade on tree2, respectively.
#' @return phylo object
#' @export cladesub
cladesub <- function(tree1,tree2,tips1,tips2=NULL,use.stem="tree1"){
	if(is.null(tips2)){
		tips2 <- tree2$tip.label
	}
	tips1 <- tree1$tip.label[unlist(phangorn::Descendants(tree1,treeio::MRCA(tree1,tips1)))]
	tips2 <- tree2$tip.label[unlist(phangorn::Descendants(tree2,treeio::MRCA(tree2,tips2)))]
	if(any(tips2 %in% setdiff(tree1$tip.label,tips1))){
		dups <- tips2[tips2 %in% setdiff(tree1$tip.label,tips1)]
		#cops <- paste0(dups,"_copy")
		tree2$tip.label[tree2$tip.label %in% dups] <- paste0(tree2$tip.label[tree2$tip.label %in% dups],"_copy")
		tips2[tips2 %in% dups] <- paste0(tips2[tips2 %in% dups],"_copy")
	}
	#return(list(tree2,tips2))
	tr1    <- ape::rtree(2,tip.label=c("TEMPTIP1","TEMPTIP2"))
	clade  <- treeio::drop.tip(tree2,setdiff(tree2$tip.label,tree2$tip.label[unlist(phangorn::Descendants(tree2,treeio::MRCA(tree2,tips2)))]),collapse.singles=F)
	clade$root.edge <- 1
	tr1$root.edge   <- 1
	tree3a <- ape::bind.tree(tree1,tr1,where=treeio::MRCA(tree1,tips1),position=0)
	tree4a <- treeio::drop.tip(tree3a,tip=tips1)
	tree3  <- ape::bind.tree(tree4a,clade,where=MRCA(tree4a,tr1$tip.label),position=0)
	tree4  <- treeio::drop.tip(tree3,tip=tr1$tip.label)
	tib4   <- tibble::as_tibble(tree4)
	if(use.stem=="tree1"){
		tib4[treeio::MRCA(tree4,tips2),'branch.length'] <- unlist(tibble::as_tibble(tree1)[treeio::MRCA(tree1,tips1),'branch.length'])
	}
	if(use.stem=="tree2"){
		tib4[treeio::MRCA(tree4,tips2),'branch.length'] <- unlist(tibble::as_tibble(tree2)[treeio::MRCA(tree2,tips2),'branch.length'])
	}
	if(use.stem=="diff"){
		age.parent.mrca1 <- max(nodeHeights(tree1)) - nodeheight(tree1,parent(tree1,treeio::MRCA(tree1,tips1)))
		#age.mrca1        <- max(nodeHeights(tree1)) - nodeheight(tree1,treeio::MRCA(tree1,tips1))
		age.mrca2        <- max(nodeHeights(tree2)) - nodeheight(tree2,treeio::MRCA(tree2,tips2))
		newlength <- age.parent.mrca1-age.mrca2
		if(newlength<0){
			stop("crown age of replacement clade is older than stem age of clade to be replaced. set 'use.stem' to value other than 'diff'")
		}
		tib4[treeio::MRCA(tree4,tips2),'branch.length'] <- newlength
	}
	tree5  <- treeio::as.phylo(tib4)
	tree5
}

