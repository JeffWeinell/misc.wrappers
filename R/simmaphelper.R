#' @title Generate transition models
#' 
#' Makes a list of all possible model schemes for transitions between a set of states. The returned object can be passed to the function fitmodels
#' 
#' @param states A character vector with the names of possible character states.
#' @param include.asymmetric Whether or not to include assymetric rate models in the list of returned models. Default TRUE; (FALSE not yet implemented)
#' @param constraint A numeric matrix with values 1 and/or 0, indicating whether or not specific transition types should be allowed or restricted, respectively. The default value of NULL allows all model types in the returned list of models. Providing a constraint matrix filters the set of models to only those that satisfy the constraint.
#' @return A list of transition matrices, for all possible partitions of pairwise transitions into parameters describing transitions.
#' @export transitionModels
transitionModels <- function(states,include.asymmetric=TRUE,constraint=NULL) {
	max.parameters=(length(states)^2)-(length(states))
	xmat     <- expand.grid(replicate(max.parameters, 1:max.parameters, simplify=FALSE))
	xmat3    <- unique(do.call(rbind,apply(xmat,1,function(x){(rank(unique(unlist(x)))[match(x, unique(unlist(x)) )])},simplify=F)))
	filter1  <- apply(xmat3,1,function(x){rank(unlist(unique(x)))}, simplify=F)
	filter2  <- lapply(lengths(filter1),function(x){c(1:x)})
	xmat4    <- xmat3[sapply(1:length(filter1),function(x){all(filter1[[x]]==filter2[[x]])}),]
	zeromats <- as.matrix(expand.grid(replicate(max.parameters, 0:1, simplify=FALSE)))
	xmat5    <- unique(do.call(rbind,lapply(2:nrow(zeromats),function(x){ unique(xmat4 %*% diag(zeromats[x,])) })))
	filter3  <- apply(xmat5,1,function(x){rank(unlist(unique(x)))}, simplify=F)
	filter4  <- lapply(lengths(filter3),function(x){c(1:x)})
	xmat6    <- xmat5[sapply(1:length(filter3),function(x){all(filter3[[x]]==filter4[[x]])}),]

	allmodels.matrix <- cbind(rep(0,nrow(xmat6)),xmat6[,(1:length(states))],rep(0,nrow(xmat6)),xmat6[,(1:length(states)+length(states))],rep(0,nrow(xmat6)))
	models.list <- lapply(1:nrow(allmodels.matrix),function(x){tempQ=matrix(as.numeric(allmodels.matrix[x,]),nrow=length(states),ncol=length(states),dimnames=list(c(states),c(states)))})
	if(!include.asymmetric){
		models.list <- models.list[sapply(models.list,isSymmetric)]
	}
	if(!is.null(constraint)){
		models.list.constrained <- lapply(models.list,function(x){x*constraint})
		models.pass <- sapply(1:length(models.list),function(x){all(models.list[[x]]==models.list.constrained[[x]])})
		models.list <- models.list[models.pass]
	}
	models.list
}

#' @title Fit data and a list of Mk models with simmap
#' 
#' Makes a list of all possible model schemes for transitions between a set of states.
#' 
#' @param tree newick character string with phylogenetic tree
#' @param priors.tips matrix with character state priors for the tips of the tree
#' @param model.transitions a square matrix or list of square matrices for the transition rate parameter schemes; a list with all possible parameter schemes can be generated with the function transitionModels
#' @param priors.root a numeric vector with the priors for the character states at the root node. Default is a uniform distribution.
#' @return A list of transition matrices, for all possible partitions of pairwise transitions into parameters describing transitions.
#' @export fitmodels
fitmodels <- function(tree,priors.tips,model.transitions,priors.root=NULL){
	priors.tips <- as.matrix(priors.tips)
	nstates <- ncol(priors.tips)
	if(!is(model.transitions,"list")){
		if(is(model.transitions,"data.frame")){
			model.transitions <- list(as.matrix(model.transitions))
		} else{
			if(is(model.transitions,"array")){
				model.transitions <- list(model.transitions)
			}
		}
	}
	if(is.null(priors.root)){
		priors.root <- (1:nstates)/nstates
	}
	tt.list <- lapply(model.transitions,function(x){tryCatch(phytools::fitMk(tree,priors.tips, model=x,pi=priors.root), error=function(e) e, warning=function(w) w)})
	#tt <- tryCatch(phytools::fitMk(tree,priors.tips, model=model.transitions), error=function(e) e, warning=function(w) w)
	models.pass <- sapply(tt.list,is,'fitMk')
	if(all(!models.pass)){
		stop("model fitting failed for all models")
	} else {
		tt.list2 <- tt.list[models.pass]
		model.transitions2 <- model.transitions[models.pass]
		#tt.list.pass <- 
		#tempfit <- list()
		#AIClist <- list()
		#length(tempfit) <- length(AIClist) <- length(tt.list2)
		tempfit <- lapply(model.transitions2,function(x){phytools::fitMk(tree,priors.tips,model=x)})
		AIClist <- lapply(tempfit,stats::AIC)
	}
	result <- list(models=model.transitions2,fitMk.objects=tempfit,AIC=unlist(AIClist))
	result
}
#' @examples
#' # load tree
#' tr <- ape::read.tree("~/Trachylepis-tree_WithAllOutgroups_NewNames2.new")
#' # load matrix of tip names and tip priors
#' tp <- as.matrix(read.csv("~/Trachylepis_parity_20Aug2018.csv", row.names=1))
#' # names of character states
#' st <- colnames(tp)
#' # generate all possible models
#' mod <- transitionModels(st)
#' ## fit the first 20 transition rates model scheme to the data, and use a uniform prior for the root state.
#' test1 <- fitmodels(tree=tr,priors.tips=as.matrix(tp),model.transitions=mod[1:20])
#' bestModel  <- test1$models[[which(test1$AIC==min(test1$AIC))]]
#' make.bestModel <- phytools::make.simmap(tr,tp,model=bestModel, nsim=1000)
#' describe.bestModel <- phytools::describe.simmap(make.bestModel,plot=F)
#' plot(describe.bestModel)
#' ## fit the first 10 transition rates model scheme to the data, and set the root state prior to c(0,1,0) for the three character states.
#' test2 <- fitmodels(tree=tr,priors.tips=as.matrix(tp),model.transitions=mod[1:10],priors.root=c(0,1,0))
#' bestModel2  <- test2$models[[which(test2$AIC==min(test2$AIC))]]
#' make.bestModel2 <- phytools::make.simmap(tr,tp,model=bestModel2, nsim=1000)
#' describe.bestModel2 <- phytools::describe.simmap(make.bestModel2,plot=F)
#' plot(describe.bestModel2)






