#' @title genind2diffs function
#' 
#' Generates a diffs matrix from a genind object
#' 
#### This function supercedes bed2diffs_v1 and bed2diffs_v2
#' @param genind.obj Input genind object (see adegenet package for description of genind class). Non-biallelic sites are dropped if they exist.
#' @param ploidy Integer indicating the ploidy of individuals
#' @param include.indv.names If the output matrix should include the names of the individuals as rownames and column names. Default FALSE
#' @param output.file Optional character string to save the result to a file.
#' @return A list of length three   with (1) a numeric matrix corresponding to the Diff matrix used for EEMS, (2) the number of individuals, (3) number of biallelic sites in genind.obj
#' @export genind2diffs
genind2diffs <- function(genind.obj,ploidy=2,include.indv.names=F,output.file=NULL){
	gen        <- genind.obj
	adegenet::ploidy(gen) <- ploidy
	stopifnot(identical(gen@type, 'codom'))
	Geno       <- gen@tab
	## Names of loci that are not biallelic
	multi.loci <- names(which(gen@loc.n.all != 2))
	## Columns corresponding to the non-biallelic loci.
	multi.cols <- gsub("\\..+","",colnames(Geno)) %in% multi.loci
	### Removing columns for loci that are not biallelic
	if(length(which(multi.cols))){
		Geno <- Geno[, -which(multi.cols)]
	}
	# Number of individuals
	nIndiv <- nrow(Geno)
	# Number of sites
	nSites <- ncol(Geno)
	# Missing data sites
	Miss   <- is.na(Geno)
	### Impute NAs with the column means of Geno
	# This is a numeric matrix with dimensions the same as Geno; all values in a column are the same and equal to the mean of non-NA values in the corresponding Geno column (= twice the allele frequencies)
	Mean   <- matrix(colMeans(Geno, na.rm = TRUE),nrow = nIndiv, ncol = nSites, byrow = TRUE)
	### Set the means in Mean to zero for observed genotypes, because these don't need to be imputed.
	# Mean[Miss == FALSE] <- 0
	Mean[!Miss] <- 0
	### Set the missing genotypes to 0 (used to be NA).
	Geno[Miss]  <- 0
	### Missing genotypes will take on the value from Mean and values at observed genotypes do not change.
	Geno     <- Geno + Mean
	## Compute similarities
	Sim      <- (Geno %*% t(Geno)) / nSites
	## Self-similarities
	SelfSim  <- diag(Sim)
	## vector of 1s
	vector1s <- rep(1, nIndiv)
	## Uses some linear algebra to generate a `diffs` matrix (version 1)
	Diffs1 <- (SelfSim %*% t(vector1s)) + (vector1s %*% t(SelfSim)) - (2 * Sim)
	rownames(Diffs1) <- colnames(Diffs1)
	######
	# Diffs2 <- matrix(0, nIndiv, nIndiv)
	#for (i in seq(nIndiv - 1)) {
	#	for (j in seq(i + 1, nIndiv)) {
	#		x <- Geno[i, ]
	#		y <- Geno[j, ]
	#		Diffs2[i, j] <- mean((x - y)^2, na.rm = TRUE)
	#		Diffs2[j, i] <- Diffs2[i, j]
	#	}
	#}
	# colnames(Diffs2) <- colnames(Diffs1)
	# rownames(Diffs2) <- rownames(Diffs1)
	### Round diffs to six digits
	Diffs1 <- round(Diffs1, digits = 6)
#	 Diffs2 <- round(Diffs2, digits = 6)
	### Check which Diffs matrix to use by taking the eigenvalue
#	eigvals.Diffs_v1 <- sort(round(eigen(Diffs1)$values, digits = 2))
#	 eigvals.Diffs_v2 <- sort(round(eigen(Diffs2)$values, digits = 2))
	# if(length(which(eigvals.Diffs_v1>0))==1){
	diffs <- Diffs1
	# } else {
	# 	if(length(which(eigvals.Diffs_v2>0))==1){
#	 	diffs <- Diffs2
	# 	} else {
	# 	stop("No Diffs matrix with one positive eigenvalue")
	# 	}
	# }
	if(!include.indv.names){
		rownames(diffs) <- NULL
		colnames(diffs) <- NULL
	}
	if(!is.null(output.file)){
		write.table(diffs,output.file,quote=F,sep=" ",col.names=F,row.names=F)
	}
	result <- list(diffs,nIndiv,nSites)
	names(result) <- c("diffs","nIndiv","nSites")
	result
} ### End genind2diffs function



