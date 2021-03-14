#' @title genind2diffs function
#' 
#' Generates a diffs matrix from a genind object
#' 
#' 
#### This function supercedes bed2diffs_v1 and bed2diffs_v2
#' @param genind.obj Input genind object (see adegenet package for description of genind class). Non-biallelic sites are dropped if they exist.
#' @param ploidy Integer indicating the ploidy of individuals
#' @param include.indv.names If the output matrix should include the names of the individuals as rownames and column names. Default FALSE
#' @param output.file Optional character string to save the result to a file.
#' @return A list of length three   with (1) a numeric matrix corresponding to the Diff matrix used for EEMS, (2) the number of individuals, (3) number of biallelic sites in genind.obj
#' @export genind2diffs
genind2diffs <- function(genind.obj,ploidy=2,include.indv.names=F,output.file=NULL){
	gen         <- genind.obj
	ploidy(gen) <- ploidy
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
	## Impute NAs with the column means (= twice the allele frequencies); a row of means; a matrix with nIndiv identical rows of means
	Mean   <- matrix(colMeans(Geno, na.rm = TRUE),nrow = nIndiv, ncol = nSites, byrow = TRUE)
	## Set the means that correspond to observed genotypes to 0
	Mean[Miss == 0] <- 0
	## Set the missing genotypes to 0 (used to be NA) 
	Geno[Miss == 1] <- 0
	Geno     <- Geno + Mean
	## Compute similarities
	Sim      <- Geno %*% t(Geno) / nSites
	## self-similarities
	SelfSim  <- diag(Sim)
	## vector of 1s
	vector1s <- rep(1, nIndiv)
	## This chunk generates a `diffs` matrix (version 1)
	Diffs1 <- SelfSim %*% t(vector1s) + vector1s %*% t(SelfSim) - 2 * Sim
	rownames(Diffs1) <- colnames(Diffs1)
	######
	Diffs2 <- matrix(0, nIndiv, nIndiv)
	for (i in seq(nIndiv - 1)) {
		for (j in seq(i + 1, nIndiv)) {
			x <- Geno[i, ]
			y <- Geno[j, ]
			Diffs2[i, j] <- mean((x - y)^2, na.rm = TRUE)
			Diffs2[j, i] <- Diffs2[i, j]
		}
	}
	colnames(Diffs2) <- colnames(Diffs1)
	rownames(Diffs2) <- rownames(Diffs1)
	### Round diffs to six digits
	Diffs1 <- round(Diffs1, digits = 6)
	Diffs2 <- round(Diffs2, digits = 6)
	### Check which Diffs matrix to use by taking the eigenvalue
	eigvals.Diffs_v1 <- sort(round(eigen(Diffs1)$values, digits = 2))
	eigvals.Diffs_v2 <- sort(round(eigen(Diffs2)$values, digits = 2))
	if(length(which(eigvals.Diffs_v1>0))==1){
		diffs <- Diffs1
	} else {
		if(length(which(eigvals.Diffs_v2>0))==1){
		diffs <- Diffs2
		} else {
		stop("No Diffs matrix with one positive eigenvalue")
		}
	}
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


#' @title runEEMs function
#' 
#' 
#' Wrapper function for running EEMS. Optionally generates the '*.outer' file from the '*.coords' file automatically.
#' 
#' 
#' @param eems.exe.path Character string with path to eems executable.
#' @param input.data An object of class genind or vcfR, or a character string with path to a diffs file.
#' @param coord Character string with path to coordinates file, which has two columns with longitude and latitude coordinates (decimal degree format) of individuals in the genind object. Columns should be space separated.
#' @param outer Either NULL (the default) or a character string with path to the '*.outer' file, which has two columns with longitude and latitude coordinates (decimal degree format) defining the perimeter of a polygon covering the region to model with EEMS. See EEMS documentation.
#' If NULL, this function will generate a "*.outer" file automatically; the user will be prompted to accept the auto-generated *.outer file prior to running eems.
#' If you don't want to use the automatically generated *.outer file, you can create a custom one here: https://www.keene.edu/campus/maps/tool/
#' @param ask.use.outer When outer = NULL, this logical determines if the auto-generated outer coordinates should be plotted, at which point the user must confirm if they want to use the coordinates. Default is TRUE.
#' @param data.type On of the following character strings: "diffs" or "vcf". Ignored if input.data is a character string. The default value NULL is coerced to "diffs".
#' @param output.dirpath Character string with path to directory where EEMS input and output files should be saved. This directory must not exist prior to running this function, but its parent directory must exist.
#' Files and directories generated by this function:
#'    /output.dirpath/
#'    /output.dirpath/data
#'    /output.dirpath/data/data.coord
#'    /output.dirpath/data/data.diffs
#'    /output.dirpath/data/data.outer
#'    /output.dirpath/mcmc    ### This will contain as many mcmc subdirectories as the value of nchains.
#'    /output.dirpath/params  ### This will contain as many parameter files as the value of nchains.
#' @param n.sites The number of sites. Determined automatically if input data is a genind or vcf object, but must be supplied if input.data is a path to a diffs matrix.
#' @param ploidy Ploidy of individuals. Default is 2.
#' @param nDemes Number of demes. Default is 300. This can be a single number or a numeric vector of length nchains with value to use for each chain.
#' @param numMCMCIter Number of MCMC iterations (i.e., chain length).Default is 10000000. This can be a single number or a numeric vector of length nchains.
#' @param numBurnIter Number of Burnin iterations (i.e., chain length). Default is 1000000. This can be a single number or a numeric vector of length nchains.
#' @param numThinIter Number of iterations to ignore before sampling the next MCMC iteration. Default is 9999. This can be a single number or a numeric vector of length nchains.
#' @param nchains Number of chains to run. Default is 3.
#' @param setup.only If the function should setup data, parameters, and output directories but not run EEMS. Default is TRUE
#' @return Nothing is returned.
#' @export runEEMs
runEEMs <- function(eems.exe.path, input.data, coord, outer=NULL, ask.use.outer=T, data.type="diffs", output.dirpath, n.sites, ploidy=2, nDemes=300, numMCMCIter = 10000000, numBurnIter = 1000000, numThinIter = 9999, nchains=3, setup.only=T){
 # data.dirname <- paste0(sample(c(letters,LETTERS,0:9),size=10,replace=T),collapse="")
 # data.dirpath <- paste0(tempdir(),"/",data.dirname)
 # dir.create(data.dirpath)
	input.data <- data

	if(dir.exists(output.dirpath)){
		stop(paste("Output directory:",output.dirpath,"already exists. Use a new output.dirpath"))
	} else {
		dir.create(output.dirpath)
	}
	### Check input arguments and stop if wrong format.
	if(length(nDemes)!=1 | length(nDemes)!= length(nchains)){
		stop("nDemes should have length 1 or length equal to nchains")
	} else {
		if(length(nDemes)==1){
			nDemes = rep(nDemes,nchains)
		}
	}
	if(length(numMCMCIter)!=1 | length(numMCMCIter)!= length(nchains)){
		stop("numMCMCIter should have length 1 or length equal to nchains")
	} else {
		if(length(numMCMCIter)==1){
			numMCMCIter = rep(numMCMCIter,nchains)
		}
	}
	if(length(numBurnIter)!=1 | length(numBurnIter)!= length(nchains)){
		stop("numBurnIter should have length 1 or length equal to nchains")
	} else {
		if(length(numBurnIter)==1){
			numBurnIter = rep(numBurnIter,nchains)
		}
	}
	if(length(numThinIter)!=1 | length(numThinIter)!= length(nchains)){
		stop("numThinIter should have length 1 or length equal to nchains")
	} else {
		if(length(numThinIter)==1){
			numThinIter = rep(numThinIter,nchains)
		}
	}
	### Create directories in output.dirpath to hold a copy of the input (data) files, parameter files, and mcmc output
	# Create input directory containing copy of input files
	input.dirpath <- paste0(output.dirpath,"/data")
	dir.create(input.dirpath)
	# Create output mcmc directory and subdirectories
	mcmcpath         <- paste0(output.dirpath,"/mcmc")
	mcmcpath.subdirs <- paste0(mcmcpath,"/",paste0(nDemes,"Demes-chain",1:nchains))
	dir.create(mcmcpath)
	sapply(X=mcmcpath.subdirs,FUN=dir.create)
	# Create directory to hold params files
	params.path <- paste0(output.dirpath,"/params")
	dir.create(params.path)
	# Copy *.coord file to input.dirpath and rename as "data.coord"
	system(paste("cp",coord,paste0(input.dirpath,"/data.coord")))
	# If path to the outer file is NULL, generate one automatically, otherwise copy it to input.dirpath and rename as "data.outer"
	if(is.null(outer)){
		data_outer <- misc.wrappers::create.outer(coords,method=1,buffer.adj=0,coords.radius=0.01,max.fractal.dimension=1.1,plot.outer=ask.use.outer,ask.use=ask.use.outer,output.path=paste("cp",outer,paste0(input.dirpath,"/data.outer")))
	} else {
		system(paste("cp",outer,paste0(input.dirpath,"/data.outer")))
	}
	if(is(input.data,"vcfR")){
		genind     <- vcfR::vcfR2genind(input.data)
		input.data <- genind
	}
	if(is(input.data,"genind")){
		data.diffs <- misc.wrappers::genind2diffs(genind.obj=input.data,output.file=paste0(input.dirpath,"/data.diffs"))
		diffs      <- data.diffs["diffs"]
		nIndiv     <- c(data.diffs["nIndiv"])
		nSites     <- data.diffs["nSites"]
	} else {
		# stop("input.data must be an object of class 'vcfR' or 'genind'")
		if(is(input.data,"character")){
			if(data.type=="vcf"){
				vcf.data   <- vcfR::read.vcfR(input.data)
				genind     <- vcfR::vcfR2genind(vcf.data)
				data.diffs <- misc.wrappers::genind2diffs(genind.obj=genind,output.file=paste0(input.dirpath,"/data.diffs"))
				diffs      <- data.diffs["diffs"]
				nIndiv     <- c(data.diffs["nIndiv"])
				nSites     <- data.diffs["nSites"]
			}
			if(data.type=="diffs"){
				data.diffs <- read.table(input.data,sep=" ",header=F)
				nIndiv     <- ncol(data.diffs)
				nSites     <- n.sites
				system(paste("cp",input.data,paste0(input.dirpath,"/data.diffs")))
			}
		}
	}
	### Generate diffs file, saving to input.dirpath with name "data.diffs"
	# data.diffs <- genind2diffs(genind.obj=genind,output.file=paste0(input.dirpath,"/data.diffs"))
	# diffs      <- data.diffs["diffs"]
	# nIndiv     <- c(data.diffs["nIndiv"])
	# nSites     <- data.diffs["nSites"]
	params.files.path <- paste0(params.path,"/params-chain",1:nchains,"_",nDemes,".ini")
	for(i in 1:nchains){
		L1 <- paste0("datapath = ",paste0(input.dirpath,"/data"))
		L2 <- paste0("mcmcpath = ",mcmcpath.subdirs[i])
		L3 <- paste0("nIndiv = ",as.integer(nIndiv))
		L4 <- paste0("nSites = ",as.integer(nSites))
		L5 <- paste0("nDemes = ",as.integer(nDemes[i]))
		L6 <- paste0("diploid = ",TRUE)
		L7 <- paste0("numMCMCIter = ",as.integer(numMCMCIter[i]))
		L8 <- paste0("numBurnIter = ",as.integer(numBurnIter[i]))
		L9 <- paste0("numThinIter = ",as.integer(numThinIter[i]))
		param.lines <- c(L1,L2,L3,L4,L5,L6,L7,L8,L9)
		writeLines(param.lines,params.files.path[i])
	}
	command.write <- c("#!/bin/bash",paste(eems.exe.path,"--params",params.files.path))
	writeLines(command.write,paste0(output.dirpath,"/runeems_snps.sh"))
	if(!setup.only){
		command.exe   <- gsub(" & $","",paste(command.write,collapse=" & "))
		system(command.exe)
	}
} ### End runEEMs function



