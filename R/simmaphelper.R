#' @title Generate transition models
#' 
#' Makes a list of all possible model schemes for transitions between a set of states. The returned object can be passed to the function fitmodels
#' 
#' @param states A vector with values representing character states.
#' @param include.asymmetric Whether or not to include assymetric rate models in the list of returned models. Default TRUE; (FALSE not yet implemented)
#' @param constraint A numeric matrix with values 1 and/or 0, indicating whether or not specific transition types should be allowed or restricted, respectively. The default value of NULL allows all model types in the returned list of models. Providing a constraint matrix filters the set of models to only those that satisfy the constraint.
#' @return A list of transition matrices, for all possible partitions of pairwise transitions into parameters describing transitions.
#' @export transitionModels
transitionModels <- function(states,include.asymmetric=TRUE,constraint=NULL) {
	if(length(states)<2){
		stop("'states' must be a vector with length > 1 ")
	}
	ns   <- length(states)
	mp   <- (ns^2)-ns
	PM   <- gtools::permutations(n=(mp+1),r=mp,v=(0:mp),repeats.allowed=T)
	RPM0 <- unique(t(apply(PM,1,misc.wrappers::reassign)))
	RPM  <- RPM0[rowSums(RPM0)!=0,]
	qmat <- ((diag(ns)-1)^2)
	models.list <- lapply(1:nrow(RPM),function(x){tempQ=qmat; tempQ[!diag(T,nrow=ns,ncol=ns)] <- RPM[x,] ; dimnames(tempQ)=list(c(states),c(states)); tempQ})

	### this section could be implemented before ennumerating permutations, which would greatly speed things up by reducing the maximum number of parameters
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

#' @title membership to equivalency matrix
#' 
#' For a vector with entries indicating partition membership, thus function returns an equivalency matrix describing vector elements belong to the same partition
#' 
#' @param v a vector
#' @return numerical equivalency relationship matrix, with 1s use indicating equivalency and zeros otherwise
#' @export equiv
equiv <- function(v){
	PM=as.matrix(expand.grid(seq(length(v)),seq(length(v))))
	EM=matrix(0,nrow=length(v),ncol=length(v))
	EM[apply(PM,1,function(x){ v[x[1]] == v[x[2]]})] <- 1
	dimnames(EM) <- replicate(2,dimnames(v))
	EM
}

#' @title reassign
#' 
#' Transforms a vector or matrix (mode character or numerical) into an integer object with the same dimension as the input object, with values on the range from 1:length(unique(x)).
#' For an input vector x1, an integer vector x2 is returned with the following properties:
#'   (1) x2 equivalency relations of x2 are identical to those in x1
#'   (2) sum(x2) <= sum(x1), for any x1 == abs(x1).
#'   (3) unique(x2[x1>0]) == 1:length(unique(x1[x1>0])); if ignore0=FALSE, then unique(x2) == 1:length(unique(x1))
#'
#' If x is a matrix then x is coerced to a vector and then coerced back to a matrix before returning.
#' 
#' 
#' @param x vector or matrix
#' @param ignore0 logical indicating if zeros should be ignored. Default true.
#' @return integer vector if x is a vector; integer matrix if x is a matrix
#' @export reassign
reassign <- function(x,ignore0=T){
	if(ignore0 && all(x==0)){
		res2 <- x
	} else {
		x0 <- x
		if(ignore0 && any(x==0)){
			x  <- x0[x0!=0]
		}
		if(is.null(dim(x))){
			# result <- rank(unique(x))[match(x,unique(x))]
			res1 <- apply(unique(misc.wrappers::equiv(x)),2,function(i){grep(1,i)})
		} else {
			#result <- matrix(rank(unique(c(x)))[match(c(x),unique(c(x)))],nrow=nrow(x),ncol=ncol(x))
			res1    <- matrix(apply(unique(misc.wrappers::equiv(x)),2,function(i){grep(1,i)}),,nrow=nrow(x),ncol=ncol(x))
		}
		if(ignore0 && any(x0==0)){
			res2 <- x0
			res2[x0!=0] <- res1
		} else {
			res2 <- res1
		}
		dimnames(res2) <- dimnames(x0)
	}
	res2
}

# These functions are very close but not quite there.
# Accurate for values of mp 1–5
# mp  <- max.parameters <- 4
equivalencyMats <- function(mp){
	# Zero matrix
	# ZM <- matrix(0,mp,mp)
	# List of starting matrices
	ZML <- lapply(1:mp,function(x){m=diag(mp); m[1:x,1:x] <- 1; m})
	## For each matrix in ZML, reorder rows and columns on the permutation set of 1:mp.
#	PM  <- gtools::permutations(n=mp,r=mp)
	# generates all possible reorderings of rows; each row reorder occurs with an identical column reorder.
	PM  <- gtools::permutations(n=mp,r=mp)
	LML <- list()
	## obtain the list of equivalence relations matrices
	for (i in 1:length(ZML)){
		ZMLi  <- lapply(1:nrow(PM),function(x){ZML[[i]][PM[x,],PM[x,]]})
		LMLi  <- unique(do.call(rbind,lapply(ZMLi,c)))
		#LML  <- rbind(LML,LMLi)
		LML   <- c(LML,list(LMLi))
	}
	LMLf <- do.call(rbind,LML)
	setparts1   <- apply(LMLf,1,function(x){matrix(x,nrow=mp,ncol=mp)},simplify=F)
	p1          <- setparts1
	numSchemes1 <- length(setparts1)
	pairsMat0   <- as.matrix(expand.grid(c(1:numSchemes1),c(1:numSchemes1)))
	pairsMat    <- pairsMat0[(pairsMat0[,1] < pairsMat0[,2]),]

	# pairwise intersections, pairwise unions, pairwise clipping, pairwise products
	i1    <- apply(pairsMat0, 1, function(x){intMat(setparts1[[x[1]]],setparts1[[x[2]]])}, simplify=F)
	u1    <- apply(pairsMat0, 1, function(x){unionMat(setparts1[[x[1]]],setparts1[[x[2]]])}, simplify=F)
	c1    <- apply(pairsMat0, 1, function(x){clipMat(setparts1[[x[1]]],setparts1[[x[2]]])}, simplify=F)
	
	
	# All expected partitions generated in the set of pairwise union sets; this is probably the fastest method.
	u1test1 <- sapply(u1,transitive)
	u1test2 <- sapply(u1,symmetric)
	u1test3 <- sapply(u1,reflexive)
	u2      <- u1[u1test1 & u1test2 & u1test3]
	u2f     <- unique(do.call(rbind,lapply(u2,c)))
	#result  <- apply(u2f,1,function(x){matrix(x,nrow=mp,ncol=mp)},simplify=F)
	
	
	# sets of matrices combined and then filtered to those unique ones
	#s1      <- c(p1,i1,u1,c1,pr1)
	s1      <- c(p1,i1,u1,c1)
	s1test1 <- sapply(s1,transitive)
	s1test2 <- sapply(s1,symmetric)
	s1test3 <- sapply(s1,reflexive)
	s2      <- s1[s1test1 & s1test2 & s1test3]
	s2f     <- unique(do.call(rbind,lapply(s2,c)))
	s2r     <- apply(s2f,1,function(x){matrix(x,nrow=mp,ncol=mp)},simplify=F)
	if(mp<6){
		result <- s2r
	} else {
		pm2    <- as.matrix(expand.grid(c(1:length(s2r)),c(1:length(s2r))))
		pr1    <- apply(pm2[pm2[,1]<pm2[,2],], 1, function(x){(s2r[[x[1]]] %*% s2r[[x[2]]])}, simplify=F)
		pr2    <- pr1[sapply(pr1,function(x){transitive(x) & symmetric(x) & reflexive(x)})]
		pr2r   <- unique(do.call(rbind,lapply(pr2,c)))
		result <- apply(pr2r,1,function(x){matrix(x,nrow=mp,ncol=mp)},simplify=F)
	}

	if(FALSE){
		# All expected partitions generated in the set of pairwise product sets
		pr1test1 <- sapply(pr1,transitive)
		pr1test2 <- sapply(pr1,symmetric)
		pr1test3 <- sapply(pr1,reflexive)
		pr2      <- pr1[pr1test1 & pr1test2 & pr1test3]
		pr2f     <- unique(do.call(rbind,lapply(pr2,c)))
	
		# nothing new from pairwise intersections
		i1test1 <- sapply(i1,transitive)
		i1test2 <- sapply(i1,symmetric)
		i1test3 <- sapply(i1,reflexive)
		i2      <- i1[i1test1 & i1test2 & i1test3]
		i2f     <- unique(do.call(rbind,lapply(i2,c)))
		
		# All expected partitions generated in the set of pairwise union sets; this is probably the fastest method.
		# Should be slower than using union sets
		c1test1 <- sapply(c1,transitive)
		c1test2 <- sapply(c1,symmetric)
		c1test3 <- sapply(c1,reflexive)
		c2      <- c1[c1test1 & c1test2 & c1test3]
		c2f     <- unique(do.call(rbind,lapply(c2,c)))
	}
	result
}

partitions <- function(mp){
	setparts <- equivalencyMats(mp)
	setparts.unique <- lapply(setparts,unique)
	setIDs <- list()
	for(i in 1:length(setparts.unique)){
		setsi0 <- apply(setparts.unique[[i]],1,function(x){which(x==1)},simplify=F)
		setsi  <- setsi0[order(sapply(setsi0,min))]
		setIDi <- sapply(1:mp,function(x){which(sapply(setsi,function(s){x %in% s}))})
		setIDs <- c(setIDs,list(setIDi))
	}
	sets.mat <- do.call(rbind,setIDs)
	sets.mat
}


unionMat <- function(m1,m2){
	m1[m1>0] <- 1
	m2[m2>0] <- 1
	nm <- nrow(m1)
	um <- (m1 + m2)
	um[um>0] <- 1
	um
}

intMat <- function(m1,m2){
	m1[m1>0] <- 1
	m2[m2>0] <- 1
	nm <- nrow(m1)
	um <- (m1 + m2)
	um[um<2] <- 0
	um[um>0] <- 1
	um
}


clipMat <- function(m1,m2){
	m1[m1>0] <- 1
	m2[m2>0] <- 1
	nm <- nrow(m1)
	cm <- (m1 + m2)
	cm[!cm==1] <- 0
	cm[diag(T,dim(cm)[1])] <- 1
	cm
}

compSet <- function(m){
	cm  <- IM(m)
	cm[m == 0] <- 1
	cm
}

# tests if a matrix is transitive
transitive <- function(m){
	m2 <- m %*% m
	#mdiff <- m2-m
	#all(mdiff<=1)
	all(unionMat(m,m2)==m)
}

# tests if a matrix is symmetric
symmetric  <- function(m){
	all(m == t(m))
}

# tests if a matrix is reflexive
reflexive  <- function(m){
	im  <- diag(dim(m)[1])
	umi <- unionMat(m,im)
	all(m==umi)
}

# tests if two matrices are disjoint
disjoint  <- function(m1,m2){
	all(ntMat(m1,m2) == 0)
}

IM <- function(m){
	diag(dim(m)[1])
}

ZM <- function(m){
	m * 0
}

bsum <- function(m2){
	apply(m2,1,max)
}








