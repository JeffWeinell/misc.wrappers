#' @title Run run_structure from SNP data in a VCF file and and plot results.
#'
#' This function is a wrapper that enables running STRUCTURE on SNP data in either a VCF file or vcfR object
#' 
#' Notes:
#' Requires installation of STRUCTURE software
#' 
#' @param x 'vcfR' object (see package::vcfR) or a character string with path to a SNPs dataset formatted according to the 'format' argument. Currently VCF or 'structure' (a type of STRUCTURE format) can be used.
#' @param format Character string indicating the format of the data. Currently only "VCF" or "structure" allowed. Other types may be added. Ignored if x is a vcfR object.
#' @param coords Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.
#' @param mainparams.path Character string with path to the mainparams file. Default is NULL, in which case the mainparams file is generated from values supplied to arguments of this function.
#' @param extraparams.path Character string with path to the extraparams file. Default is NULL, in which case the extraparams file is generated from values supplied to arguments of this function.
#' @param kmax Numerical vector with set of values to use for K. Default 40.
#' @param burnin Integer with how many initial MCMC samples to ignore. Default is 1000. NEED TO CHECK IF THIS IS REASONABLE. If this argument is NULL, then BURNIN must be defined in the file mainparams file and 'mainparams.path' must not be NULL.
#' @param numreps Chain length. Default 10000.
#' @param runs Number of times to repeat the mcmc analysis . Default 5.
#' @param ploidy Integer ≥ 1 indicating ploidy, or NULL (the default), in which case ploidy is determined automatically from the input data (only works if 'format' = "VCF" or "vcfR").
#' @param missing Integer used to code missing alleles, or NULL (the default), in which case missing data is identified automatically from the input file (only works if 'format' = "VCF" or "vcfR").
#' @param onerowperind Logical indicating if the input data file codes individuals on a single or multiple rows. Default is NULL (only works if 'format' = "VCF" or "vcfR"), in which case a temporary structure file is created and onerowperind is coerced to TRUE.
##' @param save.as Where to save the output PDF. Default is NULL.
#' @param save.in Character string with path to directory where output files should be saved. The directory will be created and should not already exist. Default is NULL, in which case output is saved to a new folder (name randomly generated) in the current directory.
##' @param tolerance Tolerance for convergence, i.e., the change in marginal likelihood required to continue.
##' @param prior Type of prior to use. Default "simple"
##' @param full Whether or not to generate output files holding variation of Q, P, and marginal likelihood, in addition to the files holding means. Default FALSE.
##' @param seed Value to use as a seed for reproducing results. Default NULL.
#' @param structure.path Character string with path to folder containing the structure executable called 'structure.py'
#' @param samplenames NULL. Not yet implemented.
#' @param cleanup Whether or not the original output files should be deleted/replaced with one, simple table holding all of the information usually spread across multiple files and tables. Default TRUE.
#' @param include.out Character vector indicating which type of files should be included as output in addition to the usual structure output. Default is c(".pdf","popfiles"). ".pdf" generates EvannoPlots and admixture barplots, and "popfiles" generates an easySFS-format popfile with individual assignments to populations for each K.
#' @param debug Logical indicating whether or not to print messages indicating the internal step of the function. Default FALSE. Typically only used for development.
#' @param ... Additional arguments passed to STRUCTURE. Not yet implemented in the future may include 'LABEL', 'POPDATA', 'POPFLAG', 'LOCDATA', 'PHENOTYPE', 'EXTRACOLS', 'MARKERNULLMES', 'RECESSIVEALLELES', 'MAPDISTANCES', 'PHASED', 'PHASEINFO', 'MARKOVPHASE', and 'NOTAMBIGUOUS'
#' @param overwrite Whether or not to overwrite previous results. Default FALSE.
#' @return List of plots
#' @export run_structure
run_structure <- function(x, format="VCF", coords=NULL, mainparams.path=NULL, extraparams.path=NULL, burnin=1000, kmax=10, numreps=10000, runs=5, ploidy=NULL, missing=NULL, onerowperind=NULL, save.in=NULL, structure.path=NULL, samplenames=NULL, cleanup=TRUE, include.out=c(".pdf","popfiles"), debug=FALSE, ...,overwrite=FALSE){
	# list with user-specified arguments
	argslist <- list(...)
	if(length(argslist)>0){
		names(argslist) <- toupper(names(argslist))
	}
	# mainparams agruments that can be determined from input data automatically when 'format' = "VCF" or "vcfR"
	if(FALSE){
		# calculated
		NUMINDS = NULL
		# calculated
		NUMLOCI = NULL
		# derived from x
		INFILE = NULL
		# derived from save.as
		OUTFILE
	}
	if(is.null(save.in)){
		save.in <- file.path(getwd(),paste(sample(c(letters,LETTERS,rep(0:9,3)),10,replace=T),collapse=""))
	}
	if(dir.exists(save.in)){
		if(!overwrite){
			stop("Output directory should not already exist. Use a different name for 'save.in' argument.")
		}
	}
	Krange=1:kmax
	if(format=="VCF" | is(x,"vcfR")){
		if(is(x,"vcfR")){
			vcf.obj <- vcf <- x
		} else {
			vcf <- x
			vcf.obj     <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F,convertNA=FALSE)
		}
		gt.mat      <- gsub(":.+","",vcf.obj@gt[,-1])
		# Detect ploidy from genotype matrix of vcf
		test.sample <- unlist(gt.mat)[!is.na(unlist(gt.mat))][1]
		ploidy      <- length(unlist(strsplit(gt.mat[1],split="[/,|]")))
		if(is.null(samplenames)){
			samplenames <- colnames(vcf.obj@gt)[-1]
		}
		### Generate structure file from vcfR object
		# vcf2structure function may not be ready
		outstr          <- file.path(getwd(),paste0(paste(sample(c(letters,LETTERS,rep(0:9,3)),30,replace=T),collapse=""),".str"))
		strconvert      <- vcf2structure(vcf=vcf.obj,out=outstr,OneRowPerIndv=FALSE)
		str.path0       <- strconvert[[1]]
		str.path        <- tools::file_path_sans_ext(str.path0)
		mainparams.df0  <- strconvert[[2]]
	} else {
		if(format=="structure"){
			str.path0       <- x
			str.path        <- tools::file_path_sans_ext(str.path0)
			mainparams.df0  <- data.frame(MAXPOPS = NA, BURNIN = NA, NUMREPS = NA, INFILE = str.path0, OUTFILE = NA, NUMINDS = NA, NUMLOCI = NA, PLOIDY = NA, MISSING = NA, ONEROWPERIND = NA, LABEL = NA, POPDATA = NA, POPFLAG = NA, LOCDATA = NA, PHENOTYPE = NA, EXTRACOLS = NA, MARKERNAMES = NA, RECESSIVEALLELES = NA, MAPDISTANCES = NA, PHASED = NA, PHASEINFO = NA, MARKOVPHASE = NA, NOTAMBIGUOUS = NA)
		} else {
			stop("Data must be supplied in VCF or structure format")
		}
	}
	## Reorganize the mainparams and extraparams data frames so that other values can be easily added, and to make it easy to convert rows to strings
	mainparams.df1           <- data.frame(definitions=rep("# define",nrow(mainparams.df0)),variables=rownames(mainparams.df0),values=mainparams.df0[,"mainparams.known"])
	rownames(mainparams.df1) <- rownames(mainparams.df0)
	## Load a tempate extraparams data frame
	# Not sure this should be "extraparams.df <- data(extraparams.template.df,package="misc.wrappers")
	extraparams.df0 <- extraparams.template.df
	### Number of individuals
	numind        <- length(samplenames)
	### text size to use for labels (individual names) in ggplots of admixture plots
	label.size    <- min((288/numind),7)
	if(debug) message("step 0")
	if(FALSE){
		### Checking that python2 exists and is executable
		if(is.null(python.path)){
			python.testpath <- find.exe.path("python")
			if(python.testpath!=1){
				python.path <- python.testpath
			} else {
				stop("No valid python identified. Set 'python.path' to location of python.")
			}
		}
	}
	### Checking that the structure program exists and is executable
	if(is.null(structure.path)){
		structure.testpath <- find.exe.path("structure.py")
		if(structure.testpath!=1){
			structure.path <- structure.testpath
		} else {
			stop("No valid path to structure identified. Set 'structure.path' to location of 'structure.py' executable file")
		}
	}
	if(!file.exists(structure.path)){
		stop(paste0(structure.path," does not exist"))
	} else {
		if(file.info(structure.path)$isdir){
			structure.path <- file.path(structure.path,"structure")
		}
	}
	if(!file.exists(structure.path)){
		stop(paste0(structure.path," does not exist"))
	} else {
		if(check.if.executable(structure.path)!=0){
			stop("structure 'structure.py' not executable")
		}
	}
	### Loading coordinates of individuals if coordinates were supplied.
	if(!is.null(coords)){
		if(is(coords,"array") | is(coords,"data.frame")){
			coords <-  coords[,c(1:2)]
		} else {
			if(file.exists(coords)){
				coords   <- read.table(coords)[,c(1:2)]
			}
		}
		colnames(coords) <- c("Lon","Lat")
		if(!is.null(rownames(coords))){
			### Check that all individuals with coords are in the vcf file, and vice versa.
			if(!all(samplenames %in% rownames(coords) & rownames(coords) %in% samplenames)){
				stop("All individuals in coords file must be in vcf")
			}
		}
		
		maxK <- min(nrow(unique(coords)),(numind-1))
	} else {
		maxK <- (numind-1)
	}
	if(kmax > maxK){
		Krange <- 1:maxK
	}
	kmax <- max(Krange)
	##### Define/create output directory and define output files.
	if(debug) message("step 2")
	outdir.temp <- save.in
	dir.create(outdir.temp)
	outfile.temp  <- file.path(outdir.temp,"structure.log")
	### Move the structure file into outdir.temp
	infile.temp0 <- file.path(getwd(),mainparams.df1["INFILE","values"])
	infile.temp1 <- file.path(outdir.temp,mainparams.df1["INFILE","values"])
	file.rename(infile.temp0,infile.temp1)
	### Move the file with IDs and samplenames into outdir.temp
	IDs.mat.temp0 <- gsub(".str$","_sampleIDs.txt",infile.temp0)
	IDs.mat.temp1 <- gsub(".str$","_sampleIDs.txt",infile.temp1)
	file.rename(IDs.mat.temp0,IDs.mat.temp1)
	###### Defining other settings for structure
	if(debug) message("step 1")
	##### Updating the mainparams file
	mainparams.df1["MAXPOPS","values"]  <- kmax
	mainparams.df1["INFILE","values"]   <- infile.temp1
	mainparams.df1["OUTFILE","values"]  <- outfile.temp
	mainparams.df1["BURNIN","values"]   <- burnin
	mainparams.df1["NUMREPS","values"]  <- numreps
	mainparams.df1["NUMINDS","values"]  <- numind
	#mainparams.df["NUMLOCI","values"] <- 
	mainparams.df1["PLOIDY","values"]   <- ploidy
	#mainparams.df1["MISSING","values"]      <- 
	#mainparams.df1["ONEROWPERIND","values"] <- 
	#mainparams.df1["LABEL","values"]        <- 
	#mainparams.df1["POPDATA","values"]      <- 
	#mainparams.df1["POPFLAG","values"]      <- 
	#mainparams.df1["LOCDATA","values"]      <- 
	#mainparams.df1["PHENOTYPE","values"]    <- 
	#mainparams.df1["EXTRACOLS","values"]    <- 
	#mainparams.df1["MARKERNAMES","values"]  <- 
	# Logical
	if(exists("recessivealleles")){
		if(recessivealleles){
			mainparams.df1["RECESSIVEALLELES","values"] <- 1
		} else {
			mainparams.df1["RECESSIVEALLELES","values"] <- 0
		}
	} else {
		mainparams.df1["RECESSIVEALLELES","values"] <- 0
	}
	if(mainparams.df1["RECESSIVEALLELES","values"]==1){
		# Logical
		if(exists("mapdistances")){
			if(mapdistances){
				mainparams.df1["MAPDISTANCES","values"] <- 1
			} else {
				mainparams.df1["MAPDISTANCES","values"] <- 0
			}
		} else {
			mainparams.df1["MAPDISTANCES","values"] <- 0
		}
		# Logical 
		if(exists("phased")){
			if(phased){
				mainparams.df1["PHASED","values"] <- 1
			} else {
				mainparams.df1["PHASED","values"] <- 0
			}
		} else {
			mainparams.df1["PHASED","values"] <- 0
		}
		# Logical 
		if(exists("phaseinfo")){
			if(phaseinfo){
				mainparams.df1["PHASEINFO","values"] <- 1
			} else {
				mainparams.df1["PHASEINFO","values"] <- 0
			}
		} else {
			mainparams.df1["PHASEINFO","values"] <- 0
		}
		# Logical
		if(exists("markovphase")){
			if(markovphase){
				mainparams.df1["MARKOVPHASE","values"] <- 1
			} else {
				mainparams.df1["MARKOVPHASE","values"] <- 0
			}
		} else {
			mainparams.df1["MARKOVPHASE","values"] <- 0
		}
		# Integer
		if(exists("notambiguous")){
			mainparams.df1["NOTAMBIGUOUS","values"] <- notambiguous
		} else {
			mainparams.df1["NOTAMBIGUOUS","values"] <- NA
		}
	} else {
		mainparams.df1["MAPDISTANCES","values"]     <- 0
		mainparams.df1["PHASED","values"]           <- 0
		mainparams.df1["PHASEINFO","values"]        <- 0
		mainparams.df1["MARKOVPHASE","values"]      <- 0
		mainparams.df1["NOTAMBIGUOUS","values"]     <- NA
	}
	# Drop variables with NA values.
	if(any(is.na(mainparams.df1[,"values"]))){
		mainparams.df2 <- mainparams.df1[!is.na(mainparams.df1[,"values"]),]
	} else {
		mainparams.df2 <- mainparams.df1
	}
	
	### This section reproduces the template extraparams dataframe that is returned if when calling the 'extraparams.template.df' data
	if(FALSE){
		logicalparams   <- c(linkage=0, noadmix=0, usepopinfo=0, locprior=0, onefst=0, inferalpha=1, popalphas=0, inferlambda=0, popspecificlambda=0)
		# freqscorr = FPRIORMEAN*FPRIORSD, but is set to NA to indicate that it should not (usually) be specified in the extraparams file or command line.
		numericalparams <- c(freqscorr=NA, alpha=0.5, lambda=1)
		# "fpriormean" and "fpriorsd" only applicable when freqscorr is TRUE
		# "log10rmin","log10rmax","log10propsd","log10rstart" only applicable when 'linkage' is TRUE
		# The following priors are only applicable when "usepopinfo" is TRUE: "gensback" (integer), "migrprior" (numeric in [0,1]; suggested 0.001 to 0.1), "pfrompopflagonly" (logical)
		# The following priors are only applicable when "locprior" is TRUE: "locispop" (logical), "locpriorinit" (numeric), "maxlocprior" (numeric; suggested value 20).
		logicalpriors        <- c(unifprioralpha=0, pfrompopflagonly=0, locispop=0)
		integerpriors        <- c(gensback=2)
		# Priors set to NA should not (usually) be specified in the extraparams file or command line, and STRUCTURE's default values will be used. The defaults are not reported in the STRUCTURE documentation.
		numericalpriors      <- c(fpriormean=NA, fpriorsd=NA, alphamax=NA, alphapriora=NA, alphapriorb=NA, log10rmin=NA, log10rmax=NA, log10propsd=NA, log10rstart=NA, migrprior=0.1, locpriorinit=NA, maxlocprior=NA)
		# Output options.
		# 'sitebysite' is applicable only when 'linkage' is TRUE
		logicaloutputoptions   <- c(printnet=1, printlambda=0, printqsum=0, sitebysite=0, printqhat=1, printlikes=0, echodata=0, ancestdist=1)
		# You may want to set 'intermedsave' to a value 5-10, to get intermediate results before numreps is reached.
		integeroutputoptions   <- c(updatefreq=0, intermedsave=5)
		numericaloutputoptions <- c(ancestpint=0.9, numboxes=0.01)
		### Miscellaneous options
		logicalmisc   <- c(computeprob=1, startatpopinfo=0, randomize=1, reporthitrate=1)
		# admburnin applicable only when 'recombine' = 1; metrofreq only applicable when noadmix=0
		integermisc   <- c(admburnin=500, seed=NA, metrofreq=500)
		# 'alphapropsd' affects mixing
		numericalmisc <- c(alphapropsd=0.5)
		# vector and data frame with all of the logicals
		all.extra.logicals    <- c(logicalparams,logicalpriors,logicaloutputoptions,logicalmisc)
		all.extra.logicals.df <- data.frame(definitions="#define",variables=toupper(names(all.extra.logicals)),values=unname(all.extra.logicals),type="boolean")
		# data frame with all of the numericals (doubles)
		all.extra.numericals  <- c(numericalparams,numericalpriors,numericaloutputoptions,numericalmisc)
		all.extra.numericals.df <- data.frame(definitions="#define",variables=toupper(names(all.extra.numericals)),values=unname(all.extra.numericals),type="double")
		# data frame with all of the integers
		all.extra.integers  <- c(integerpriors,integeroutputoptions,integermisc)
		all.extra.integers.df  <- data.frame(definitions="#define",variables=toupper(names(all.extra.integers)),values=unname(all.extra.integers),type="integer")
		# data frame with all possible parameters that could be included in the extraparams file
		extraparams.template.df <- rbind(all.extra.logicals.df,all.extra.numericals.df,all.extra.integers.df)
	}
	#### Updating extraparams.df
	if(length(argslist)>0){
		if(any(extraparams.df0[,"variables"] %in% names(argslist))){
			args2update <- intersect(extraparams.df0[,"variables"] %in% names(argslist))
			for(i in 1:length(args2update)){
				extraparams.df0[args2update[i],"variables"] <- argslist[args2update[i]]
			}
		}
	}
	#### Drop variables with values NA
	if(any(is.na(extraparams.df0[,"values"]))){
		extraparams.df1 <- extraparams.df0[!is.na(extraparams.df0[,"values"]),]
	} else {
		extraparams.df1 <- extraparams.df1
	}
	#### Change the mode of the values column to 'character', because a mixture of value types is used.
	mode(extraparams.df1[,"values"]) <- "character"
	if(FALSE){
		### Drop variables if they would be ignored given the values of the other variables. STRUCTURE should do this already, so skipping for now.
		### Ancestry Models
		#		No admixture
		#		Admixture
		#		Linkage
		#		Prior Pop Information
		#			LOCPRIOR with PopData or LocData
		#			USEPOPINFO with MIGRPRIOR
		#			USEPOPINFO with PopFlag=1 for learning samples and MIGRPRIOR>0
	}
	# Drop the "type" column from extraparams data frame
	extraparams.df <- extraparams.df1[,c("definitions","variables","values")]
	mainparams.df  <- mainparams.df2
	if(is.null(mainparams.path)){
		mainparams.path <- file.path(outdir.temp,"mainparams")
		# write mainparams to mainparams.path
		#mainparams.strings <- sapply(X=1:nrow(mainparams),FUN=function(x) paste(x,collapse="\t") ), paste(,collapse="\t")
		#writeLines(mainparams.df,mainparams.path,colnames=FALSE,rownames=TRUE,sep="/t",quote=FALSE)
		write.table(mainparams.df,file=mainparams.path,col.names=F,row.names=F,sep="\t",quote=FALSE)
	} else {
		# copies mainparams into the output directory
		file.copy(from=mainparams.path,to=file.path(outdir.temp,"mainparams"))
	}
	if(is.null(extraparams.path)){
		extraparams.path <- file.path(outdir.temp,"extraparams")
		write.table(extraparams.df,file=extraparams.path,col.names=F,row.names=F,sep="\t",quote=FALSE)
	} else {
		# copies extraparams into the output directory
		file.copy(from=extraparams.path,to=file.path(outdir.temp,"extraparams"))
	}
	#for(i in 1:runs){
	#	for(K in Krange){
	#		outfile.K <- paste0(tools::file_path_sans_ext(outfile.temp),"_K",K,"_run",i,".log")
	#		command1     <- paste0(structure.path," -K ",K," -m ",mainparams.path," -e ",extraparams.path," -o ",outfile.K)
	#		run.command1 <- system(command1)
	#	}
	#}
	for(i in 1:runs){
		for(K in Krange){
			outfile.K <- paste0(tools::file_path_sans_ext(outfile.temp),"_K",K,"r",i,".log")
	#		outfile.K <- paste0(tools::file_path_sans_ext(outfile.temp),"_K",K,".log")
			command1     <- paste0(structure.path," -K ",K," -m ",mainparams.path," -e ",extraparams.path," -o ",outfile.K)
			run.command1 <- system(command1)
		}
	}
	if(".pdf" %in% include.out){
		#### Make Evanno Plots. Once the other functions for plotting are made I'll change save.as to NULL.
		ggEvanno <- EvannoPlots(input.dir=outdir.temp, save.as=file.path(outdir.temp, "EvannoPlots.pdf"))
		#### Admixture barplots
		admixture.barplot <- admixturePlots(xdir=outdir.temp, userun=c(1:runs))
		#### Geographic-iterpolation of admixture coefficients

	}
	if("popfiles" %in% include.out){
		#### For each K, generate a popfile that can be used by easysfs
		for(i in 1:kmax){
			create_popfile(xdir=outdir.temp, K=i)
		}
	}
	## Remove 0-byte files in outdir.temp, if they exist
	outfiles       <- list.files(outdir.temp,full.names=T)
	outfileinfo.df <- file.info(outfiles,extra_cols = FALSE)
	if(any(outfileinfo.df[,"size"]==0)){
		doremoval <- file.remove(outfiles[which(outfileinfo.df[,"size"]==0)])
	}
	#### Close any open graphics windows, if they exist
	if(!is.null(dev.list())){
		graphics.off()
	}
	return(NULL)
	stop("function not ready for implementation")
	### Make Assignment plots

	### Make interpolated-assignment maps if coords is not null

	#qfiles         <- list.files(outdir.temp, full.names=T, pattern="log_f$")
	#admixture.plot <- admixturePlots(x=qfiles, labels=samplenames)
	
	# Barplots; probably better to use the ggplot method
	# qplots <- pophelper::plotQ(qlist=aqlist, exportpath=outdir.temp)
	##### CODE BELOW HERE NOT UPDATED #####
	#### Load the marginal likelihood scores (the metric of fit for each K)
	logpaths <- list.files(outdir.temp,pattern="^.+\\.log$",full.names=TRUE)
	loglines <- lapply(logpaths,readLines)
	logfiles <- basename(logpaths)
	### Value of K associated with each output logfile. This info is extracted from the filename and corresponds to the number immediately before the ".log" extension.
	#Kvals.out.temp <- do.call(rbind,strsplit(logfiles,split=".",fixed=TRUE))
	#Kvals.out      <- as.numeric(Kvals.out.temp[,rev(1:ncol(Kvals.out.temp))][,2])
	KLogs          <- as.numeric(sapply(logfiles,function(x){unlist(strsplit(x,split=".",fixed=TRUE))[2]}))
	### replicate of K associated with each log file
	repsLogs       <- as.numeric(gsub("^rep","",gsub("\\..+$","",logfiles)))
	### Matrix with marginal likelihood from each rep and each K
	margL.mat <- do.call(rbind,lapply(Krange,FUN=function(x){A=unlist(loglines[which(KLogs==x)]);B=A[grep("^Marginal Likelihood = .+",A)]; as.numeric(gsub("Marginal Likelihood = ","",B,fixed=T))}))
	rownames(margL.mat) <- paste0("K",Krange)
	colnames(margL.mat) <- paste0("rep",1:numreps)
	mean.margL      <- apply(margL.mat,MARGIN=1,FUN=mean,na.rm=TRUE)
	range.margL.mat <- do.call(rbind,lapply(X=1:nrow(margL.mat),FUN=function(x){range(margL.mat[x,],na.rm=TRUE)}))
	if(debug) message("step 3")
	margL.df      <- data.frame(marginalLikelihood=unname(unlist(c(margL.mat))),Kval=rep(Krange, numreps))
	margL.df$Kval <- factor(margL.df$Kval, levels=c(1:nrow(margL.df)))
	margLPlot     <- ggplot2::ggplot(margL.df, ggplot2::aes(x=Kval, y=marginalLikelihood)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("marginal likelihood (",numreps," replicates) vs. number of ancestral populations (K)"), x="Number of ancestral populations", y = "marginalLikelihood") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
	if(debug) message("step 4")
	#### Load the Q matrices
	qpaths          <- list.files(outdir.temp,pattern="^.+\\.meanQ$",full.names=TRUE)
	qmats.list      <- lapply(qpaths,read.table)
	qfiles   <- basename(qpaths)
	### K value associated with each q matrix
	KQmats          <- sapply(qmats.list,ncol)
	### replicate of K associated with each q matrix
	repsQmats       <- as.numeric(gsub("^rep","",gsub("\\..+$","",qfiles)))
	### for a particular K, the replicate with the highest marginal likelihood
	bestReps        <- unname(sapply(1:max(Krange),FUN=function(x){which(margL.mat[x,]==max(margL.mat[x,]))[1]}))
	### for a particular K, the qmatrix of the replicate with the highest marginal likelihood
	qmats.list.best <- lapply(1:max(Krange),FUN=function(x){qmats.list[[which(KQmats==x)[bestReps[(x)]]]]})
	#### List holding population assignment probabilities for each K
	# Could also try using the mean among qmatrices, but that could be problematic if cluster1 is treated as cluster2 in another iteration.
	if(debug) message("step 5")
	slist          <- qmats.list.best
	Krange.plot    <- setdiff(Krange,1)
	admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
	assignmentPlot <- list(); length(assignmentPlot)  <- length(Krange.plot)
	if(debug) message("step 6")
	if(!is.null(coords)){
		mapplot        <- list(); length(mapplot)     <- length(Krange.plot)
		x.min <- min((coords[,1]-0.5))
		x.max <- max((coords[,1]+0.5))
		y.min <- min((coords[,2]-0.5))
		y.max <- max((coords[,2]+0.5))
		world_sf      <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
		world_sp      <- rnaturalearth::ne_countries(scale=10,returnclass="sp")
		#current_sf    <- sf::st_crop(world_sf,xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max)
		#current.gg.sf <- ggplot2::geom_sf(data=current_sf,colour = "black", fill = NA)
	} else {
		mapplot <- NULL
	}
	if(debug) message("step 7")
	for(K in Krange.plot){
		if(debug) message(paste0("K=",K," step 7.1"))
		i=(K-1)
		if(K <= 15){
			myCols          <- goodcolors2(n=K)
		}
		if(K>15){
			myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
		}
		if(debug) message(cat("\r",paste0("K=",K," step 7.2")))
		q.matrix  <- slist[[K]]
		rownames(q.matrix) <- samplenames
		colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
		indv.pop           <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))})
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
		if(debug) message(cat("\r",paste0("K=",K," step 7.3")))
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Admixture Proportion",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
		admixturePlot[[i]]  <- posterior.gg
		if(debug) message(cat("\r",paste0("K=",K," step 7.4")))
		indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels             <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignment.K         <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
		assignmentPlot[[i]]  <- assignment.K
		if(debug) message(cat("\r",paste0("K=",K," step 7.5")))
		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols, 9)
			if(debug) message(cat("\r",paste0("K=",K," step 7.6")))
			mapplot.i       <- suppressWarnings(tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE, map.polygon=world_sp))
			if(debug) message(cat("\r",paste0("K=",K," step 7.7")))
			mapplot2.i      <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA) + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
			if(debug) message(cat("\r",paste0("K=",K," step 7.8")))
			mapplot[[i]]    <- mapplot2.i + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			#mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
		} else {
			mapplot[[i]] <- NULL
		}
	}
	if(debug) message("step 8")
	result <- c(list(margLPlot),admixturePlot,assignmentPlot,mapplot)
	#### Save the PDF
	if(debug) message("step 9")
	if(".pdf" %in% include.out){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	if(debug) message("step 10")
	if(format=="VCF"){
		# Delete the temporary structure file
		if(file.exists(str.path0)){
			remove.str <- file.remove(str.path0)
		}
	}
	
	if(debug) message("step 11")
	#### Generate compiled output files so that the other output files can be deleted.
	if(".margLlog" %in% include.out){
		### Organizing a single data frame to hold all data from all "*.log" files generated by structure. This data frame will be written to a single file to replace the many '*.log' files.
		log.df <- NULL
		for(i in 1:length(KLogs)){
			loglines.i   <- loglines[[i]]
			loglines2.i  <- loglines.i[-c(grep("^Iteration",loglines.i))]
			loglines3.i  <- loglines2.i[1:(length(loglines2.i)-3)]
			init.i.lines <- loglines3.i[grep("^M",loglines3.i)]
			init.margL.i <- sapply(1:length(init.i.lines),function(x){as.numeric(rev(unlist(strsplit(init.i.lines[x],split=" ")))[1])})
			init.i.mat   <- unname(cbind((1:length(init.margL.i)*-1),init.margL.i,rep("--",length(init.margL.i)),rep("--",length(init.margL.i))))
			iter.i.lines <- loglines3.i[-grep("^M",loglines3.i)]
			iter.i.mat   <- do.call(rbind,lapply(1:length(iter.i.lines),function(x){unlist(strsplit(iter.i.lines[x],split=" "))}))
			logmat.i     <- rbind(init.i.mat,iter.i.mat)
			logmat2.i    <- cbind(logmat.i,rep(KLogs[i],nrow(logmat.i)),rep(repsLogs[i],nrow(logmat.i)))
			colnames(logmat2.i) <- c("Iteration", "Marginal_Likelihood", "delta_Marginal_Likelihood", "Iteration_Time_secs","K","replicate")
			log.i.df <- data.frame(K_Replicate_Iteration=paste(logmat2.i[,"K"],logmat2.i[,"replicate"],logmat2.i[,"Iteration"],sep=":"), MarginalLikelihood=logmat2.i[,"Marginal_Likelihood"],deltaMarginalLikelihood=logmat2.i[,"delta_Marginal_Likelihood"], IterationTime.seconds = logmat2.i[,"Iteration_Time_secs"])
			log.df <- rbind(log.df,log.i.df)
		}
		write.table(x=log.df,file=paste0(tools::file_path_sans_ext(save.as),".margLlog"),row.names=F,col.names=T,quote=F,sep="\t")
	}
	if(debug) message("step 12")
	if(".Qlog" %in% include.out){
		### Organizing a single data frame to hold all q matrices for all replicates of every K. This data frame will be written to a single file so that the many '*.meanQ' files produced by fastStrucure can be deleted.
		q.df <- NULL
		for(i in 1:length(KQmats)){
			qmatrix.i    <- qmats.list[[i]]
			rownames(qmatrix.i) <- samplenames
			colnames(qmatrix.i) <- paste0("cluster",1:ncol(qmatrix.i))
			q.i.df <- data.frame(individual=rep(rownames(qmatrix.i),ncol(qmatrix.i)), cluster=as.numeric(gsub("^cluster","",rep(colnames(qmatrix.i),each=nrow(qmatrix.i)))), assignment=c(unlist(unname(qmatrix.i))),K=rep(KQmats[i],(KQmats[i]*numind)),replicate=rep(repsQmats[i],(KQmats[i]*numind)))
			# posterior.df <- data.frame(indv=rep(rownames(qmatrix.i),ncol(qmatrix.i)), pop=rep(colnames(qmatrix.i),each=nrow(qmatrix.i)), assignment=c(unlist(unname(qmatrix.i))),replicate=rep(repsQmats[i],(KQmats[i]*numind)))
			q.df <- rbind(q.df,q.i.df)
		}
		write.table(x=q.df,file=paste0(tools::file_path_sans_ext(save.as),".Qlog"),row.names=F,col.names=T,quote=F,sep="\t")
	}
	if(debug) message("step 13")
	if(".Plog" %in% include.out){
		### Organizing a single data frame to hold all data from all "*.meanP" files generated by structure. This data frame will be written to a single file to replace the many '*.meanP' files.
		ppaths <- list.files(outdir.temp,pattern="^.+\\.meanP$",full.names=TRUE)
		if(length(ppaths)>0){
			pmats.list <- lapply(ppaths,read.table)
			pfiles     <- basename(ppaths)
			### Value of K associated with each output logfile. This info is extracted from the filename and corresponds to the number immediately before the ".log" extension.
			KPmats          <- as.numeric(sapply(pfiles,function(x){unlist(strsplit(x,split=".",fixed=TRUE))[2]}))
			### replicate of K associated with each log file
			repsPmats       <- as.numeric(gsub("^rep","",gsub("\\..+$","",pfiles)))
			p.df <- NULL
			for(i in 1:length(KPmats)){
				pmatrix.i    <- pmats.list[[i]]
				# rownames(qmatrix.i) <- samplenames
				colnames(pmatrix.i) <- paste0("cluster",1:ncol(pmatrix.i))
				p.i.df <- data.frame(cluster=as.numeric(gsub("^cluster","",rep(colnames(pmatrix.i),each=nrow(pmatrix.i)))), value=c(unlist(unname(pmatrix.i))),K=rep(KPmats[i],(KPmats[i]*nrow(pmatrix.i))),replicate=rep(repsPmats[i],(KPmats[i]*nrow(pmatrix.i))))
				# posterior.df <- data.frame(indv=rep(rownames(qmatrix.i),ncol(qmatrix.i)), pop=rep(colnames(qmatrix.i),each=nrow(qmatrix.i)), assignment=c(unlist(unname(qmatrix.i))),replicate=rep(repsQmats[i],(KQmats[i]*numind)))
				p.df <- rbind(p.df,p.i.df)
			}
			write.table(x=p.df,file=paste0(tools::file_path_sans_ext(save.as),".Plog"),row.names=F,col.names=T,quote=F,sep="\t")
		}
	}
	# Delete the folder with .log, .meanQ, and .meanP files; the info from these are compiled in a single file.
	if(debug) message("step 14")
	if(cleanup){
		system(paste0("rm -R ",outdir.temp))
	}
	result
}
#' @examples
#' library(misc.wrappers)
#' ## Example 1:
#' # Path to VCF with SNPs
#' vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
#' # Run structure 30 times each for K=1-10
#' run_structure(x=vcf.path,kmax=10,reps=30,save.as="fs_simK4.pdf",include.out=c(".pdf"))
#' 
#' ## Example 2: Same SNP dataset as example 1, but here we also provide Lon/Lat coordinates of individuals to geographically interpolate admixture coefficients.
#' vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4.vcf.gz")
#' coords.path <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4_coords.txt")
#' run_structure(x=vcf.path, coords=coords.path, kmax=10, reps=30, save.as="fs_simK4_withCoords.pdf", include.out=c(".pdf"))


#' @title Convert VCF file/object to STRUCTURE file (for SNP data)
#' 
#' Converts/writes a VCF file or vcfR object containing SNPs into STRUCTURE format, which is written to a file.
#' Note: Does not yet handle 'recessive alleles' or phase information (see STRUCTURE manual). This information is optional for STRUCTURE but may be of interest to some users.
#' 
#' @param vcf Character string with path to input VCF, or an object of class vcfR.
#' @param MarkerNames Either a logical (TRUE or FALSE), or a character string with marker names. If TRUE (the default), the marker names are constructed from the CHROM and POS columns of the VCF file: 'CHROM_POS' (e.g. "102_4").
#' @param InterMarkerDists Either a logical indicating if intermarker distances should be included from CHROM and POS columns of the VCF, or a vector of integers with inter-marker distances to use. Default is FALSE. A warning is generated if multiple sites per locus are present in the input VCF and 'InterMarkerDists' is FALSE. See STRUCTURE manual for details on supplying inter-marker distances.
#' @param IndvNames Logical indicating if the first column of the STRUCTURE file should contain the names of the individuals, or a character string with names to to use. Default is TRUE, in which case the names will be the names used in the VCF file.
#' @param PopData Either NULL or an integer vector indicating user-defined population assignments of individuals. The default is NULL (population data not included in output file). If non-NULL, PopData is written in the second column of the output file.
#' @param PopFlag Either NULL or a logical vector indicating whether or not STRUCTURE should use the PopData information for the particular individual. Default is NULL, in which case the PopFlag column is not written in the output file. If supplied, PopData must be non-NULL.
#' @param OneRowPerIndv Logical specifying if each genotypes should be written to a single row or multiple rows. Default TRUE, but STRUCTURE works equivalently with either format.
#' @param LocData NULL (the default) or a vector of integers specifying user-defined sampling locality for each individual.
#' @param Phenotype NULL (the default) or a vector of integers specifying the value of a phenotype of interest for each individual.
#' @param OtherData NULL (the default) or either a matrix or data frame with as many rows as individuals and columns specifying any other information of interest associated with individuals.
#' @param MissingData Number to use for missing-data. Default is -9.
#' @param out Path where output structure file should be written (default NULL). The mainparams file is also written using the same name but with the extension '.params'.
#' @return A list with [[1]] Value of 'out', which contains the path to the output file, and [[2]] values that should be used for some of the parameters in the mainparams file.
#' @export vcf2structure
vcf2structure <- function(vcf, IndvNames=TRUE, OneRowPerIndv=TRUE, MarkerNames=TRUE, MissingData=c(-9), out=NULL, InterMarkerDists=FALSE, PopData=NULL, PopFlag=NULL, LocData=NULL, Phenotype=NULL, OtherData=NULL){
	if(is.null(out)){
		out  <- file.path(getwd(),paste0(paste(sample(c(letters,LETTERS,rep(0:9,3)),30,replace=T),collapse=""),".str"))
	}
	# Check class of vcf argument value.
	if(is(vcf,"vcfR")){
		vcf.obj <- vcf
	} else {
		if(is(vcf,"character")){
			# read VCF file into R as a vcfR object
			vcf.obj     <- vcfR::read.vcfR(file=vcf,convertNA=FALSE)
		} else {
			stop("vcf argument must be eith a character string with path to VCF file, or a vcfR object")
		}
	}
	## Data frame with site (SNP) specific information
	fix.mat <- vcf.obj@fix
	# Extract genotypes for each individual and hold genotypes in a character matrix; if polyploid, haplotypes are separated by "/" or "|"
	gt.mat      <- gsub(":.+","",vcf.obj@gt[,-1])
	# Detect ploidy from gt.mat
	test.sample <- unlist(gt.mat)[!is.na(unlist(gt.mat))][1]
	ploidy      <- length(unlist(strsplit(gt.mat[1],split="[/,|]",fixed=F)))
	# Names of individuals
	samplenames <- colnames(vcf.obj@gt)[-1]
	# Replace "." with "_", because structure does not allow the former to be used in samplenames
	samplenames <- gsub(".","_",samplenames,fixed=T)
	# Number of individuals
	numind      <- length(samplenames)
	# Unique ID for each site. If MarkerNames is TRUE then these will be used for MarkerNames.
	SiteNames <- paste(fix.mat[,"CHROM"],fix.mat[,"POS"],sep="_")
	# Number of sites
	nsites <- nrow(fix.mat)
	# locus ID number from CHROM column
	uniqueLoci <- unique(unique(fix.mat[,"CHROM"]))
	# Number of loci (aka chromosomes or linkage groups)
	nloci  <- length(uniqueLoci)
	if(nsites > nloci & !InterMarkerDists){
		warning("Some loci have multiple sites, but InterMarkerDists is FALSE, and therefore STRUCTURE will incorrectly treat linked sites as unlinked")
	}
	# Optional information rows and columns need to be considered for inclusion on output.
	### MarkerNames row
	if(is(MarkerNames,"logical")){
		if(MarkerNames){
			MarkerNames <- paste(fix.mat[,"CHROM"],fix.mat[,"POS"],sep="_")
		} else {
			MarkerNames <- NULL
		}
	} else {
		if(length(MarkerNames)!=nsites){
			stop("If MarkerNames is not a logical it should be a vector with length equal to the number of sites")
		}
	}
	### InterMarkerDists row
	if(is(InterMarkerDists,"logical")){
		if(InterMarkerDists){
			if(nsites==nloci){
				InterMarkerDists <- rep(-1,nsites)
			} else {
				imd.list         <- lapply(X=1:nloci,FUN=function(x){grep(paste0("^",uniqueLoci[x],"_"),SiteNames,value=F)})
				InterMarkerDists <- unlist(lapply(X=imd.list,FUN=function(x){x[1]=-1;x}))
			}
		} else {
			InterMarkerDists <- NULL
		}
	} else {
		if(length(InterMarkerDists)!=nsites){
			stop("InterMarkerDists must be TRUE, FALSE, or an integer vector with length equal to the number of sites")
		}
	}
	### RecessiveAlleles information not implemented
	RecessiveAlleles <- NULL
	if(FALSE){
		stop("Not implemented. This information not required by STRUCTURE and usually not typically used for SNP data.")
	}
	### PhaseInfo information not yet implemented but might add later. This info might be extractable from some VCFs. See vcf.obj@meta
	PhaseInfo <- NULL
	if(FALSE){
		stop("Not yet implemented")
		vcf.meta <- vcf.obj@meta
	}
	### IndvNames; hopefully there isnt a limit on the lengths of names
	if(is(IndvNames,"logical")){
		if(IndvNames){
			IndvNames <- samplenames
		} else {
			IndvNames <- NULL
		}
	} else {
		if(length(IndvNames)!=numind){
			stop("IndvNames must be TRUE, FALSE, or a character or numeric vector with a length equal to the number of individuals")
		}
	}
	### PopData and PopFlag
	if(!is.null(PopData)){
		if(length(PopData)!=numind){
			stop("If PopData is non-NULL it must be a character or numeric vector with a length equal to the number of individuals")
		} else {
			if(!is.null(PopFlag)){
				if(length(PopFlag)!=numind | !all(test %in% c(0,1))){
					stop("If PopFlag is non-NULL it must be a character or numeric vector with a length equal to the number of individuals, with each value either 0 (don't use PopData for individual) or 1 (use PopData for individual)")
				}
			}
		}
	} else {
		# If PopData is NULL, ignore the user-defined value for PopFlag argument and set PopFlag to NULL.
		PopFlag <- NULL
	}
	if(!is.null(LocData)){
		if(length(LocData)!=numind | any(is.na(suppressWarnings(as.integer(LocData))))){
			stop("If LocData is non-NULL it must be a vector of integers with a length equal to the number of individuals")
		}
	}
	if(!is.null(Phenotype)){
		if(length(Phenotype)!=numind | any(is.na(suppressWarnings(as.integer(Phenotype))))){
			stop("If Phenotype is non-NULL it must be a vector of integers with a length equal to the number of individuals")
		}
	}
	if(!is.null(OtherData)){
		if(is(OtherData,c("matrix","array","data.frame"))){
			if(nrow(OtherData)!=numind){
				stop("If OtherData is non-NULL it must either be (1) a vector with its length equal to number of individuals, or (2) a matrix (or data frame coercible to a matrix) with as many rows as genotypes individuals.")
			} else {
				ExtraCols <- ncol(OtherData)
				OtherData <- do.call(rbind,lapply(X=1:nrow(OtherData),FUN=function(x){paste(OtherData[x,],collapse=" ")}))
			}
		} else {
			if(length(OtherData)!=numind){
				stop("If OtherData is non-NULL it must either be (1) a vector with its length equal to number of individuals, or (2) a matrix (or data frame coercible to a matrix) with as many rows as genotypes individuals.")
			} else {
				ExtraCols <- 1
			}
		}
	} else {
		ExtraCols <- 0
	}
	
	### Replace each individual name with an integer ID and save a two-column table that translates from integer ID to full name
	IntIDs     <- 1:length(IndvNames)
	IntIDs.mat <- cbind(IntIDs,IndvNames)
	write.table(IntIDs.mat,paste0(tools::file_path_sans_ext(out),"_sampleIDs.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
	
	if(OneRowPerIndv){
		### Transpose the genotype matrix
		t.gt.mat <- t(gt.mat)
		### Split genotype strings by VCF haplotype separators ("/" or "|") and then reorganize into matrix with n columns per individual for n-ploidy individuals.
		mat.temp0 <- unname(do.call(rbind,lapply(1:nrow(t.gt.mat), FUN=function(x){unlist(strsplit(t.gt.mat[x,],split="[/,|]"))})))
		### Replace VCF missing data value "." with value supplied by MissingData argument.
		mat.temp1 <- gsub(".",MissingData,mat.temp0,fixed=TRUE)
		###
		# MarkerNames <- rep(MarkerNames,each=ploidy)
		# 
	} else {
		### split genotype strings by VCF haplotype separators ("/" or "|") and then reorganize into matrix with n columns per individual for n-ploidy individuals.
		mat.temp0 <- t(unname(do.call(rbind,lapply(1:nrow(gt.mat), FUN=function(x){unlist(strsplit(gt.mat[x,],split="[/,|]"))}))))
		### replace VCF missing data value "." with value supplied by MissingData argument
		mat.temp1 <- gsub(".",MissingData,mat.temp0,fixed=TRUE)
		#### at this point the structure genotype matrix (mat.temp1) is in the correct format but columns and rows are not named. Now the other optional information rows and columns need to be considered for inclusion on output.
		# Repeat each value for the non-genotype columns twice 
		IntIDs      <- rep(IntIDs,each=ploidy)
		PopData     <- rep(PopData,each=ploidy)
		PopFlag     <- rep(PopFlag,each=ploidy)
		LocData     <- rep(LocData,each=ploidy)
		Phenotype   <- rep(Phenotype,each=ploidy)
		OtherData   <- rep(OtherData,each=ploidy)
	}
	
	### Individual/Genotype matrix
	indv.gt.mat      <- cbind(IntIDs,PopData,PopFlag,LocData,Phenotype,OtherData,mat.temp1)
	### Converting each row of the Individual/Genotype matrix into a space-delineated character string
	indv.gt.strings  <- unlist(lapply(1:nrow(indv.gt.mat),FUN=function(x){paste(indv.gt.mat[x,],collapse="\t")}))
	### Collapsing each of MarkerNames, RecessiveAlleles, InterMarkerDists, and PhaseInfo into a space-delineated character string
	MarkerNames         <- paste(MarkerNames,collapse="\t")
	RecessiveAlleles    <- paste(RecessiveAlleles,collapse="\t")
	InterMarkerDists    <- paste(InterMarkerDists,collapse="\t")
	PhaseInfo           <- paste(PhaseInfo,collapse="\t")
	### Combining the other rows (MarkerNames, RecessiveAlleles, InterMarkerDist, and PhaseInformation) with indv.gt.strings
	result.strings   <- c(MarkerNames,RecessiveAlleles,InterMarkerDists,indv.gt.strings,PhaseInfo)
	### Remove empty strings
	result.strings   <- result.strings[nchar(result.strings) > 0]
	### Write output STRUCTURE file
	writeLines(text=result.strings,con=out)
	### Generate the mainparams file
	# For now just return but dont write a list holding the values of the parameters in mainparams.
	mainparams <- matrix(data=NA, nrow=23,ncol=1)
	rownames(mainparams) <- c("MAXPOPS","BURNIN","NUMREPS","INFILE","OUTFILE","NUMINDS","NUMLOCI", "PLOIDY","MISSING","ONEROWPERIND","LABEL","POPDATA","POPFLAG","LOCDATA","PHENOTYPE","EXTRACOLS","MARKERNAMES","RECESSIVEALLELES","MAPDISTANCES","PHASED","PHASEINFO","MARKOVPHASE","NOTAMBIGUOUS")
	mainparams["INFILE",]           <- basename(out)
	mainparams["NUMINDS",]          <- numind
	mainparams["NUMLOCI",]          <- nsites
	mainparams["PLOIDY",]           <- ploidy
	mainparams["MISSING",]          <- MissingData
	mainparams["ONEROWPERIND",]     <- if(OneRowPerIndv) 1 else 0
	mainparams["LABEL",]            <- if(is.null(IndvNames)) 0 else 1
	mainparams["POPDATA",]          <- if(is.null(PopData)) 0 else 1
	mainparams["POPFLAG",]          <- if(is.null(PopFlag)) 0 else 1
	mainparams["LOCDATA",]          <- if(is.null(LocData)) 0 else 1
	mainparams["PHENOTYPE",]        <- if(is.null(Phenotype)) 0 else 1
	mainparams["EXTRACOLS",]        <- ExtraCols
	mainparams["MARKERNAMES",]      <- if(MarkerNames=="") 0 else 1
	mainparams["RECESSIVEALLELES",] <- if(RecessiveAlleles=="") 0 else 1
	mainparams["MAPDISTANCES",]     <- if(InterMarkerDists=="") 0 else 1
	mainparams["PHASEINFO",]        <- if(PhaseInfo=="") 0 else 1
	mp.df <- data.frame(mainparams.known=c(mainparams),row.names=rownames(mainparams))
	return(list(output.path=out,mp.df))
}
#' @examples
#' library(misc.wrappers)

#' @title Generate Evanno Method Plots
#'
#' This function uses functions from pophelper and ggplot2 to generate a pdf with the four Evanno Method Plots.
#' 
#' @param input.dir directory with output files of STRUCTURE analysis. Default is the current directory.
#' @param save.as Where to save the output PDF. Default is "evannoPlots.pdf" in the current directory.
#' @return g table with gpplots arranged in a 2x2 grid
#' @export EvannoPlots
EvannoPlots <- function(input.dir=getwd(),save.as="EvannoPlots.pdf"){
	# Character vector with paths to output structure files
	qfiles <- list.files(input.dir, full.names=T, pattern="log_f$")
	# Read qfiles into R
	qlist  <- pophelper::readQ(files=qfiles)
	# tabulate structure results in a data frame
	tr1   <- pophelper::tabulateQ(qlist)
	sr1   <- pophelper::summariseQ(tr1)
	evStr.df <- suppressWarnings(pophelper::evannoMethodStructure(sr1))
	# Define K as a factor
	evStr.df[,"k"]    <- factor(evStr.df[,"k"], levels=c(1:nrow(evStr.df)))
	# ln probability of the data
	#elpdmean.plot0 <-  ggplot2::ggplot() + ggplot2::geom_boxplot(data=evStr.df, ggplot2::aes(x=k, y=elpdmean), fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("A"), x="K", y = "L(K)") + ggplot2::xlim(levels(evStr.df$k))
	elpdmean.plot0 <-  ggplot2::ggplot() + ggplot2::geom_point(data=evStr.df, ggplot2::aes(x=k, y=elpdmean), fill='black',color='black',shape=21) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("A"), x="K", y = "L(K)") + ggplot2::xlim(levels(evStr.df$k))
	elpdmean.plot  <-  elpdmean.plot0    + ggplot2::geom_segment(data=evStr.df,ggplot2::aes(x=k, y=elpdmin, xend=k, yend=elpdmax),arrow=grid::arrow(angle=90,ends="both",length=grid::unit(0.1, "inches")))
	# first derivative
	#lnk1.plot     <-  ggplot2::ggplot() + ggplot2::geom_boxplot(data=evStr.df[2:nrow(evStr.df),], ggplot2::aes(x=k, y=lnk1), fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("B"), x="K", y = "L'(K)") + ggplot2::xlim(levels(evStr.df$k))
	lnk1.plot0     <-  ggplot2::ggplot() + ggplot2::geom_point(data=evStr.df[2:nrow(evStr.df),], ggplot2::aes(x=k, y=lnk1), fill='black',color='black',shape=21) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("B"), x="K", y = "L'(K)") + ggplot2::xlim(levels(evStr.df$k))
	lnk1.plot      <-  lnk1.plot0 + ggplot2::geom_segment(data=evStr.df[2:nrow(evStr.df),],ggplot2::aes(x=k, y=lnk1min, xend=k, yend=lnk1max),arrow=grid::arrow(angle=90,ends="both",length=grid::unit(0.1, "inches")))
	# second derivative
	#lnk2.plot     <-  ggplot2::ggplot() + ggplot2::geom_boxplot(data=evStr.df[2:(nrow(evStr.df)-1),], ggplot2::aes(x=k, y=lnk2), fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("C"), x="K", y = "|L''(K)|") + ggplot2::xlim(levels(evStr.df$k))
	lnk2.plot0     <-  ggplot2::ggplot() + ggplot2::geom_point(data=evStr.df[2:(nrow(evStr.df)-1),], ggplot2::aes(x=k, y=lnk2), fill='black',color='black',shape=21) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("C"), x="K", y = "|L''(K)|") + ggplot2::xlim(levels(evStr.df$k))
	lnk2.plot      <-  lnk2.plot0 + ggplot2::geom_segment(data=evStr.df[2:(nrow(evStr.df)-1),], ggplot2::aes(x=k, y=lnk2min, xend=k, yend=lnk2max),arrow=grid::arrow(angle=90,ends="both",length=grid::unit(0.1, "inches")))
	# delta K
	if(!all(is.na(evStr.df[,"deltaK"]))){
		# deltaK.plot <-  ggplot2::ggplot(data=evStr.df[!is.na(evStr.df[,"deltaK"]),], ggplot2::aes(x=k, y=deltaK, na.rm = TRUE)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("D"), x="K", y = "delta K") + ggplot2::xlim(levels(evStr.df$k))
		deltaK.plot <-  ggplot2::ggplot(data=evStr.df[!is.na(evStr.df[,"deltaK"]),], ggplot2::aes(x=k, y=deltaK, na.rm = TRUE)) + ggplot2::geom_point(fill='black',color='black',shape=21) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("D"), x="K", y = "delta K") + ggplot2::xlim(levels(evStr.df$k))
	} else {
		blank.df    <- data.frame(x=c(0,0.5,1),y=c(0,0.5,1),text=c("","delta K is NA for all K",""))
		deltaK.plot <- ggplot2::ggplot(blank.df,ggplot2::aes(x=x,y=y)) + ggplot2::geom_text(blank.df,mapping=ggplot2::aes(x=x,y=y,label=text)) + ggplot2::labs(title="", x="", y = "") + ggplot2::theme_void()
	}
	# Define viewport for saving
	vp           <- grid::viewport(height=grid::unit(0.95,"npc"),width=grid::unit(0.95,"npc"))
	# hold plots as a list of arranged grobs
	evannoPlots0 <- lapply(X=lapply(list(elpdmean.plot,lnk1.plot,lnk2.plot,deltaK.plot), FUN=ggplot2::ggplotGrob), FUN=gridExtra::arrangeGrob,vp=vp)
	# arrange the plots in a 2x2 grid
	evannoPlots  <- gridExtra::arrangeGrob(grobs=evannoPlots0,layout_matrix=rbind(c(1,2),c(3,4)),respect=TRUE,top="Evanno Method Plots")
	# grid::grid.newpage(); grid::grid.draw(evannoPlots)
	# Save pdf
	if(!is.null(save.as)){
		pdf(height=6,width=8,file=save.as,onefile=TRUE)
			grid::grid.draw(evannoPlots)
		dev.off()
	}
	evannoPlots
}


#' @title Generate a pdf with admixture plots for each K
#' 
#' Takes as input the output files from strucure runs and generates a pdf with plots of admixture.
#' 
#' @param x Character vector with paths to input files. Ignored if 'xdir' is non-NULL. These are either '*.log_f' files generated by STRUCTURE, or a sinle '*.Qlog' file generated by 'run_DAPC', 'run_fastStructure', or 'run_SNMF' functions.
#' @param xdir Optional path to the directory containing all STRUCTURE output files plus the '.*_sampleIDs.txt' file generated by run_structure. Default NULL. If supplied, this overides save.as, labels, and x arguments.
#' @param labels Character vector with names of individuals. Default NULL. Ignored if 'xdir' non-NULL.
#' @param save.as Character string with path/name to use for the output PDF file. Default is to save the output in the current directory with the name "admixturePlots.pdf".
#' @param userun Number or numerical vector indicating which runs should be used for admixture plots. Default is 1 (the first run). When multiple runs are used, the mean is used across runs after aligning clusters.
#' @return NULL; generates pdf with a barplot of admixture for each K
#' @export admixturePlots
admixturePlots <- function(x, xdir=NULL, labels=NULL, save.as=file.path(getwd(),"admixturePlots.pdf"), userun=1,save=T){
	if(!is.null(xdir)){
		# Character vector with paths to input files
		qfiles         <- c(list.files(xdir, full.names=T, pattern="log_f$"), list.files(xdir, full.names=T, pattern="Qlog$"))
		save.as        <- file.path(xdir,"admixturePlots.pdf")
		labels <- NULL
	} else {
		qfiles <- x
	}
	if(length(grep("log_f$", qfiles))>0){
		# Read qfiles into R
		qlist      <- pophelper::readQ(files=qfiles)
		names.df.path  <- list.files(xdir, full.names=T, pattern="_sampleIDs.txt$")
		samplenames.df <- read.table(names.df.path, header=T)
		samplenames    <- samplenames.df$IndvNames
	} else {
		if(length(grep("Qlog$", qfiles))==1){
			qtab <- read.table(qfiles, header=T,sep="\t")
			if(!"replicate" %in% colnames(qtab)){
				if(nrow(unique(qtab[,c("individual", "cluster","K")])) == nrow(qtab)){
					qtab[,"replicate"] <- 1
				} else {
					if(floor(nrow(qtab)/nrow(unique(qtab[,c("individual", "cluster","K")]))) == ceiling(nrow(qtab)/nrow(unique(qtab[,c("individual", "cluster","K")])))){
						qtab <- qtab[order(test[,"K"], test[,"cluster"], test[,"individual"]),]
						qtab[,"replicate"] <- rep(1:(nrow(qtab)/nrow(unique(qtab[,c("individual", "cluster","K")]))),each=nrow(unique(qtab[,c("individual", "cluster","K")])))
					} else {
						stop("number of replicates differ not equal")
					}
				}
			}
			if(all(is.na(qtab[,"replicate"]))){
				qtab$replicate <- 1
			}
			qid         <- unique(qtab[,c("K","replicate")])
			### List of Q matrices (data frames) from the rectangular Qlog data frame
			qlist        <- lapply(X=1:nrow(qid), FUN=function(x){resB=do.call(cbind,lapply(X=1:qid[x,"K"], FUN=function(z) { A=qtab[which(qtab[,"replicate"]==qid[x,"replicate"] & qtab[,"K"]==qid[x,"K"] ),]; A[A[,"cluster"]==z, "assignment"]})); colnames(resB)=paste0("Cluster",1:qid[x,"K"]); rownames(resB)=1:nrow(resB); resC <- as.data.frame(resB); resC})
			names(qlist) <- sapply(X=1:nrow(qid),FUN=function(x) paste0("structure_K",qid[x,"K"],"r",qid[x,"replicate"],".log_f"))
			samplenames  <- qtab[1:length(unique(qtab[,"individual"])),"individual"]
			#### Need to set the following attributes for each data frame in qlist or else alignK function will fail.
			# ind, k, loci, reps, elpd, mvll, vll # each is an integer except elpd and mvll, which are numeric floats
			# This loop adds arbitrary values for required attributes, allowing us to use the alignK function for qmatrices programs other than STRUCTURE
			for(i in 1:length(qlist)){
				attr(qlist[[i]], 'ind')    <- nrow(qlist[[i]])
				attr(qlist[[i]], 'k')      <- ncol(qlist[[i]])
				attr(qlist[[i]], 'loci')   <- 5000
				attr(qlist[[i]], 'burnin') <- 9999
				attr(qlist[[i]], 'reps')   <- 99999
				attr(qlist[[i]], 'elpd')   <- -99999.9
				attr(qlist[[i]], 'mvll')   <- -99999.9
				attr(qlist[[i]], 'vll')    <- 2000
			}
		} else {
			stop("Q matrices must be supplied in either 'log_f' files produced by STRUCTURE (or 'run_structure' wrapper function) or 'Qlog' files like those produced by 'run_SNMF', 'run_fastStructure', or 'run_DAPC' functions")
		}
	}
	numind     <- nrow(qlist[[1]])
	if(!is.null(labels)){
		samplenames <- labels
	} else {
		if(is.null(xdir)){
			samplenames <- 1:numind
		}
	}
	label.size <- min((288/numind),7)
	# Set rownames of each matrix in qlist as the names of samples
	qlist2 <- lapply(X=1:length(qlist), FUN=function(x) {A=qlist[[x]]; rownames(A)=samplenames; A})
	# Aligned qlist
	aqlist <- pophelper::alignK(qlist)
	# List of aligned q matrices with rownames as samplenames
	aqlist2   <- lapply(X=1:length(aqlist), FUN=function(x) {A=aqlist[[x]]; rownames(A)=samplenames; A})
	aqlist2.K <- sapply(aqlist2, ncol)
	kmax      <- max(aqlist2.K)
	numruns   <- max(table(aqlist2.K))
	### Just use the first run if userun are al greater than the number of runs
	if(any(userun <= numruns)){
		userun <- userun[userun %in% 1:numruns]
	} else {
		userun <- 1
	}
	#for(i in 1:unique(aqlist2.K)){
	#	aqlist.i <- do.call(cbind,aqlist2[aqlist2.K == i])
	#	margin()
	#}
	# List of aligned q matrices for first run of each K. Will edit this to use instead the mean across runs.
	slist      <- aqlist2[match(1:kmax, sapply(aqlist2, ncol))]
	slist.list <- list(); length(slist.list) <- numruns
#	slist      <- qlist2[match(1:kmax,sapply(qlist2,ncol))+(userun-1)]
	#slist2     <- list(); length(slist2) <- length(slist)
	for(i in 1:numruns){
		slist.list[[i]] <- aqlist2[match(1:kmax,sapply(aqlist2, ncol)) + (i-1)]
	}
	slist.list2 <- slist.list[userun]
	if(length(slist.list2)>1){
		### Average across runs
		for(i in 1:kmax){
			X    <- lapply(X=slist.list2, FUN=function(x){ c(t(x[[i]])) })
			Y    <- do.call(rbind, X)
			Z    <- colMeans(Y)
			zmat <- matrix(Z,ncol=i,byrow=T)
			rownames(zmat) <- samplenames
			colnames(zmat) <- colnames(slist[[i]])
			zdf            <- as.data.frame(zmat)
			slist2[[i]]    <- zdf
			#Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
			#apply(X=Y, MARGIN=c(1, 2), FUN=mean, na.rm = TRUE)
		}
	} else {
		slist2    <- slist.list2[[1]]
	}
	Krange      <- 1:kmax
	Krange.plot <- 2:kmax
	# empty list to hold admixture plots
	admixturePlot <- list(); length(admixturePlot)   <- length(Krange.plot)
	for(K in 2:kmax){
		i=(K-1)
		if(K <= 15){
			myCols          <- goodcolors2(n=15)[1:K]
		}
		if(K>15){
			myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
		}
		q.matrix  <- slist2[[K]]
		#q.matrix  <- slist[[K]]
		# test <- ape::ladderize(phangorn::NJ(dist(q.matrix)))
		### Attempt to reorder individuals in the barplot by their population assignment proportions.
		indv.pop           <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))[1]})
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
		#if(debug) message(cat("\r",paste0("K=",K," step 7.3")))
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Admixture Proportion",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
		admixturePlot[[i]]  <- posterior.gg
	}
	if(save){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
			lapply(X=admixturePlot, FUN=print)
		dev.off()
	}
	admixturePlot
}
#' @examples
#' outdir.temp    <- "PATH/TO/STRUCTURE/OUTPUT"
#' names.df.path  <- list.files(outdir.temp, full.names=T, pattern="_sampleIDs.txt$")
#' samplenames.df <- read.table(names.df.path, header=T)
#' samplenames    <- samplenames.df$IndvNames
#' # Character vector with paths to output structure files
#' qfiles         <- list.files(outdir.temp, full.names=T, pattern="log_f$")
#' admixturePlots(x=qfiles,labels=samplenames,userun=c(1:5))
#'
#' admixturePlots(xdir="PATH/TO/STRUCTURE/OUTPUT",userun=c(1:5))
#' 




#' @title Generate a pdf with assignment plots for each K
#' 
#' Takes as input the output files from strucure runs and generates a pdf with plots of proportion of each population represented in each individual.
#' 
#' @param x Character vector with paths to input files. Ignored if 'xdir' is non-NULL. These are either '*.log_f' files generated by STRUCTURE, or a sinle '*.Qlog' file generated by 'run_DAPC', 'run_fastStructure', or 'run_SNMF' functions.
#' @param xdir Optional path to the directory containing all STRUCTURE output files plus the '.*_sampleIDs.txt' file generated by run_structure. Default NULL. If supplied, this overides save.as, labels, and x arguments.
#' @param labels Character vector with names of individuals. Default NULL. Ignored if 'xdir' non-NULL.
#' @param save.as Character string with path/name to use for the output PDF file. Default is to save the output in the current directory with the name "assignmentPlots.pdf".
#' @param userun Number or numerical vector indicating which runs should be used for assignment plots. Default is 1 (the first run). When multiple runs are used, the mean is used across runs after aligning clusters.
#' @return NULL; generates pdf with assignment plots for each K
#' @export admixturePlots
assignmentPlots <- function(x, xdir=NULL, labels=NULL, save.as=file.path(getwd(),"assignmentPlots.pdf"), userun=1){
	if(!is.null(xdir)){
		# Character vector with paths to input files
		qfiles         <- c(list.files(xdir, full.names=T, pattern="log_f$"), list.files(xdir, full.names=T, pattern="Qlog$"))
		save.as        <- file.path(xdir,"assignmentPlots.pdf")
		labels <- NULL
	} else {
		qfiles <- x
	}
	if(length(grep("log_f$", qfiles))>0){
		# Read qfiles into R
		qlist      <- pophelper::readQ(files=qfiles)
		names.df.path  <- list.files(xdir, full.names=T, pattern="_sampleIDs.txt$")
		samplenames.df <- read.table(names.df.path, header=T)
		samplenames    <- samplenames.df$IndvNames
	} else {
		if(length(grep("Qlog$", qfiles))==1){
			qtab <- read.table(qfiles, header=T,sep="\t")
			if(all(is.na(qtab[,"replicate"]))){
				qtab$replicate <- 1
			}
			qid         <- unique(qtab[,c("K","replicate")])
			### List of Q matrices (data frames) from the rectangular Qlog data frame
			qlist        <- lapply(X=1:nrow(qid), FUN=function(x){resB=do.call(cbind,lapply(X=1:qid[x,"K"], FUN=function(z) { A=qtab[which(qtab[,"replicate"]==qid[x,"replicate"] & qtab[,"K"]==qid[x,"K"] ),]; A[A[,"cluster"]==z, "assignment"]})); colnames(resB)=paste0("Cluster",1:qid[x,"K"]); rownames(resB)=1:nrow(resB); resC <- as.data.frame(resB); resC})
			names(qlist) <- sapply(X=1:nrow(qid),FUN=function(x) paste0("structure_K",qid[x,"K"],"r",qid[x,"replicate"],".log_f"))
			samplenames  <- qtab[1:length(unique(qtab[,"individual"])),"individual"]
			#### Need to set the following attributes for each data frame in qlist or else alignK function will fail.
			# ind, k, loci, reps, elpd, mvll, vll # each is an integer except elpd and mvll, which are numeric floats
			# This loop adds arbitrary values for required attributes, allowing us to use the alignK function for qmatrices programs other than STRUCTURE
			for(i in 1:length(qlist)){
				attr(qlist[[i]], 'ind')    <- nrow(qlist[[i]])
				attr(qlist[[i]], 'k')      <- ncol(qlist[[i]])
				attr(qlist[[i]], 'loci')   <- 5000
				attr(qlist[[i]], 'burnin') <- 9999
				attr(qlist[[i]], 'reps')   <- 99999
				attr(qlist[[i]], 'elpd')   <- -99999.9
				attr(qlist[[i]], 'mvll')   <- -99999.9
				attr(qlist[[i]], 'vll')    <- 2000
			}
		} else {
			stop("Q matrices must be supplied in either 'log_f' files produced by STRUCTURE (or 'run_structure' wrapper function) or 'Qlog' files like those produced by 'run_SNMF', 'run_fastStructure', or 'run_DAPC' functions")
		}
	}
	numind     <- nrow(qlist[[1]])
	if(!is.null(labels)){
		samplenames <- labels
	} else {
		if(is.null(xdir)){
			samplenames <- 1:numind
		}
	}
	label.size <- min((288/numind),7)
	# Set rownames of each matrix in qlist as the names of samples
	qlist2 <- lapply(X=1:length(qlist), FUN=function(x) {A=qlist[[x]]; rownames(A)=samplenames; A})
	# Aligned qlist
	aqlist <- pophelper::alignK(qlist)
	# List of aligned q matrices with rownames as samplenames
	aqlist2   <- lapply(X=1:length(aqlist), FUN=function(x) {A=aqlist[[x]]; rownames(A)=samplenames; A})
	aqlist2.K <- sapply(aqlist2, ncol)
	kmax      <- max(aqlist2.K)
	numruns   <- max(table(aqlist2.K))
	### Just use the first run if userun are al greater than the number of runs
	if(any(userun <= numruns)){
		userun <- userun[userun %in% 1:numruns]
	} else {
		userun <- 1
	}
	#for(i in 1:unique(aqlist2.K)){
	#	aqlist.i <- do.call(cbind,aqlist2[aqlist2.K == i])
	#	margin()
	#}
	# List of aligned q matrices for first run of each K. Will edit this to use instead the mean across runs.
	slist      <- aqlist2[match(1:kmax, sapply(aqlist2, ncol))]
	slist.list <- list(); length(slist.list) <- numruns
#	slist      <- qlist2[match(1:kmax,sapply(qlist2,ncol))+(userun-1)]
	#slist2     <- list(); length(slist2) <- length(slist)
	for(i in 1:numruns){
		slist.list[[i]] <- aqlist2[match(1:kmax,sapply(aqlist2, ncol)) + (i-1)]
	}
	slist.list2 <- slist.list[userun]
	if(length(slist.list2)>1){
		### Average across runs
		for(i in 1:kmax){
			X    <- lapply(X=slist.list2, FUN=function(x){ c(t(x[[i]])) })
			Y    <- do.call(rbind, X)
			Z    <- colMeans(Y)
			zmat <- matrix(Z,ncol=i,byrow=T)
			rownames(zmat) <- samplenames
			colnames(zmat) <- colnames(slist[[i]])
			zdf            <- as.data.frame(zmat)
			slist2[[i]]    <- zdf
			#Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
			#apply(X=Y, MARGIN=c(1, 2), FUN=mean, na.rm = TRUE)
		}
	} else {
		slist2    <- slist.list2[[1]]
	}
	Krange      <- 1:kmax
	Krange.plot <- 2:kmax
	# empty list to hold assignment plots
	assignmentPlot <- list(); length(assignmentPlot)   <- length(Krange.plot)
	for(K in 2:kmax){
		i=(K-1)
		if(K <= 15){
			myCols          <- goodcolors2(n=15)[1:K]
		}
		if(K>15){
			myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
		}
		q.matrix  <- slist2[[K]]
		#q.matrix  <- slist[[K]]
		# test <- ape::ladderize(phangorn::NJ(dist(q.matrix)))
		### Attempt to reorder individuals in the barplot by their population assignment proportions.
		indv.pop           <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))[1]})
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
		#if(debug) message(cat("\r",paste0("K=",K," step 7.3")))
		indv.maxPosterior    <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels               <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignmentPlot[[i]]  <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv, fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
	}
	pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=assignmentPlot, FUN=print)
	dev.off()
	assignmentPlot
}
#' @examples
#' outdir.temp    <- "PATH/TO/STRUCTURE/OUTPUT"
#' names.df.path  <- list.files(outdir.temp, full.names=T, pattern="_sampleIDs.txt$")
#' samplenames.df <- read.table(names.df.path, header=T)
#' samplenames    <- samplenames.df$IndvNames
#' # Character vector with paths to output structure files
#' qfiles         <- list.files(outdir.temp, full.names=T, pattern="log_f$")
#' assignmentPlots(x=qfiles,labels=samplenames,userun=c(1:5))
#'
#' assignmentPlots(xdir="PATH/TO/STRUCTURE/OUTPUT",userun=c(1:5))
#' 






#' @title Generate a popfile from STRUCTURE output
#' 
#' Takes as input the output files from strucutre runs and generates a pdf with plots of admixture
#' 
#' @param xdir Path to the directory containing all STRUCTURE output files plus the '.*_sampleIDs.txt' file generated by run_structure. Default NULL. If supplied, this overides save.as, method, labels, and x arguments.
#' @param K Number indicating the value of K to use.
#' @return NULL; generates pdf with a barplot of admixture for each K
#' @export create_popfile
create_popfile <- function(xdir=NULL, K=2){
	names.df.path  <- list.files(xdir, full.names=T, pattern="_sampleIDs.txt$")
	samplenames.df <- read.table(names.df.path, header=T)
	samplenames    <- samplenames.df$IndvNames
	# Character vector with paths to output structure files
	qfiles         <- list.files(xdir, full.names=T, pattern="log_f$")
	# Read qfiles into R
	qlist      <- pophelper::readQ(files=qfiles)
	numind     <- nrow(qlist[[1]])
	# Set rownames of each matrix in qlist as the names of samples
	qlist2 <- lapply(X=1:length(qlist), FUN=function(x) {A=qlist[[x]]; rownames(A)=samplenames; A})
	# Aligned qlist
	aqlist    <- pophelper::alignK(qlist)
	# List of aligned q matrices with rownames as samplenames
	aqlist2   <- lapply(X=1:length(aqlist), FUN=function(x) {A=aqlist[[x]]; rownames(A)=samplenames; A})
	aqlist2.K <- sapply(aqlist2, ncol)
	kmax      <- max(aqlist2.K)
	numruns   <- length(qlist)/kmax
	slist.list <- list(); length(slist.list) <- numruns
	for(i in 1:numruns){
		slist.list[[i]] <- aqlist2[match(1:kmax,sapply(aqlist2, ncol)) + (i-1)]
	}
	slist2 <- list(); length(slist2) <- kmax
	if(length(slist.list)>1){
		### Average across runs
		for(i in 1:kmax){
			X    <- lapply(X=slist.list, FUN=function(x){ c(t(x[[i]])) })
			Y    <- do.call(rbind, X)
			Z    <- colMeans(Y)
			zmat <- matrix(Z,ncol=i,byrow=T)
			rownames(zmat) <- samplenames
			colnames(zmat) <- paste0("Cluster",1:ncol(zmat))
			zdf            <- as.data.frame(zmat)
			slist2[[i]]    <- zdf
		}
	} else {
		slist2    <- slist.list[[1]]
	}
	result.list <- list(); length(result.list) <- length(K)
	for(i in K){
		save.as      <- file.path(xdir,paste0("popfile_K",i,"_from_structure.txt"))
		q.matrix     <- slist2[[i]]
		indv.pop     <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))[1]})
		if(length(unique(indv.pop)) < i){
			warning(paste0("For K=",i," some populations are represented only as a minority in admixed individuals. skipping this K."))
			next
		}
		result       <- data.frame(indv=names(indv.pop),assignment=paste0("pop",indv.pop))
		write.table(x=result,file=save.as,col.names=F,row.names=F,quote=F,sep="\t")
		result.list[[i]] <- result
	}
	result.list
}
#' @examples
#' create_popfile(xdir="PATH/TO/STRUCTURE/OUTPUT",K=2)






#' @title Generate a pdf with interpolated admixture coefficients for each K
#' 
#' Takes as input a coordinates file and one or more output files ('log_f' files) from structure, or the 'Qlog' file from run_DAPC, run_SNMF, or run_fastStructure, and generates a pdf with plots of admixture interpolated.
#' 
#' @param x Character vector with paths to input files. Ignored if 'xdir' is non-NULL. These are either '*.log_f' files generated by STRUCTURE, or a sinle '*.Qlog' file generated by 'run_DAPC', 'run_fastStructure', or 'run_SNMF' functions.
#' @param coords Character string with path to coordinates file
#' @param xdir Optional path to the directory containing all STRUCTURE output files plus the '.*_sampleIDs.txt' file generated by run_structure. Default NULL. If supplied, this overides save.as, labels, and x arguments.
#' @param labels Character vector with names of individuals. Default NULL. Ignored if 'xdir' non-NULL.
#' @param save.as Character string with path/name to use for the output PDF file. Default is to save the output in the current directory with the name "admixturePlots.pdf".
#' @param userun Number or numerical vector indicating which runs should be used for admixture plots. Default is 1 (the first run). When multiple runs are used, the mean is used across runs after aligning clusters.
#' @return NULL; generates pdf with a barplot of admixture for each K
#' @export admixtureMap
admixtureMap <- function(x=NULL, xdir=NULL, coords, labels=NULL, save.as=file.path(getwd(),"admixtureMaps.pdf"), userun=1){
	### Begin importing qlist
	if(!is.null(xdir)){
		# Character vector with paths to input files
		qfiles         <- c(list.files(xdir, full.names=T, pattern="log_f$"), list.files(xdir, full.names=T, pattern="Qlog$"))
		save.as        <- file.path(xdir,"admixtureMaps.pdf")
		labels <- NULL
	} else {
		qfiles <- x
	}
	if(length(grep("log_f$", qfiles))>0){
		# Read qfiles into R
		qlist      <- pophelper::readQ(files=qfiles)
		names.df.path  <- list.files(xdir, full.names=T, pattern="_sampleIDs.txt$")
		samplenames.df <- read.table(names.df.path, header=T)
		samplenames    <- samplenames.df$IndvNames
	} else {
		if(length(grep("Qlog$", qfiles))==1){
			qtab <- read.table(qfiles, header=T,sep="\t")
			if(all(is.na(qtab[,"replicate"]))){
				qtab$replicate <- 1
			}
			qid         <- unique(qtab[,c("K","replicate")])
			### List of Q matrices (data frames) from the rectangular Qlog data frame
			qlist        <- lapply(X=1:nrow(qid), FUN=function(x){resB=do.call(cbind,lapply(X=1:qid[x,"K"], FUN=function(z) { A=qtab[which(qtab[,"replicate"]==qid[x,"replicate"] & qtab[,"K"]==qid[x,"K"] ),]; A[A[,"cluster"]==z, "assignment"]})); colnames(resB)=paste0("Cluster",1:qid[x,"K"]); rownames(resB)=1:nrow(resB); resC <- as.data.frame(resB); resC})
			names(qlist) <- sapply(X=1:nrow(qid),FUN=function(x) paste0("structure_K",qid[x,"K"],"r",qid[x,"replicate"],".log_f"))
			samplenames  <- qtab[1:length(unique(qtab[,"individual"])),"individual"]
			#### Need to set the following attributes for each data frame in qlist or else alignK function will fail.
			# ind, k, loci, reps, elpd, mvll, vll # each is an integer except elpd and mvll, which are numeric floats
			# This loop adds arbitrary values for required attributes, allowing us to use the alignK function for qmatrices programs other than STRUCTURE
			for(i in 1:length(qlist)){
				attr(qlist[[i]], 'ind')    <- nrow(qlist[[i]])
				attr(qlist[[i]], 'k')      <- ncol(qlist[[i]])
				attr(qlist[[i]], 'loci')   <- 5000
				attr(qlist[[i]], 'burnin') <- 9999
				attr(qlist[[i]], 'reps')   <- 99999
				attr(qlist[[i]], 'elpd')   <- -99999.9
				attr(qlist[[i]], 'mvll')   <- -99999.9
				attr(qlist[[i]], 'vll')    <- 2000
			}
		} else {
			stop("Q matrices must be supplied in either 'log_f' files produced by STRUCTURE (or 'run_structure' wrapper function) or 'Qlog' files like those produced by 'run_SNMF', 'run_fastStructure', or 'run_DAPC' functions")
		}
	} ### Done importing qlist
	### Begin processing the qlist to obtain the qmatrix that will be used for spatial interpolation
	numind     <- nrow(qlist[[1]])
	if(!is.null(labels)){
		samplenames <- labels
	} else {
		if(is.null(xdir)){
			samplenames <- 1:numind
		}
	}
	label.size <- min((288/numind),7)
	# Set rownames of each matrix in qlist as the names of samples
	qlist2 <- lapply(X=1:length(qlist), FUN=function(x) {A=qlist[[x]]; rownames(A)=samplenames; A})
	# Aligned qlist
	aqlist <- pophelper::alignK(qlist)
	# List of aligned q matrices with rownames as samplenames
	aqlist2   <- lapply(X=1:length(aqlist), FUN=function(x) {A=aqlist[[x]]; rownames(A)=samplenames; A})
	aqlist2.K <- sapply(aqlist2, ncol)
	kmax      <- max(aqlist2.K)
	numruns   <- max(table(aqlist2.K))
	### Just use the first run if userun are all greater than the number of runs
	if(any(userun <= numruns)){
		userun <- userun[userun %in% 1:numruns]
	} else {
		userun <- 1
	}
	#for(i in 1:unique(aqlist2.K)){
	#	aqlist.i <- do.call(cbind,aqlist2[aqlist2.K == i])
	#	margin()
	#}
	# List of aligned q matrices for first run of each K. Will edit this to use instead the mean across runs.
	slist      <- aqlist2[match(1:kmax, sapply(aqlist2, ncol))]
	slist.list <- list(); length(slist.list) <- numruns
#	slist      <- qlist2[match(1:kmax,sapply(qlist2,ncol))+(userun-1)]
	slist2     <- list(); length(slist2) <- length(slist)
	for(i in 1:numruns){
		slist.list[[i]] <- aqlist2[match(1:kmax,sapply(aqlist2, ncol)) + (i-1)]
	}
	slist.list2 <- slist.list[userun]
	if(length(slist.list2)>1){
		### Average across runs
		for(i in 1:kmax){
			X    <- lapply(X=slist.list2, FUN=function(x){ c(t(x[[i]])) })
			Y    <- do.call(rbind, X)
			Z    <- colMeans(Y)
			zmat <- matrix(Z,ncol=i,byrow=T)
			rownames(zmat) <- samplenames
			colnames(zmat) <- colnames(slist[[i]])
			zdf            <- as.data.frame(zmat)
			slist2[[i]]    <- zdf
		}
	} else {
		slist2    <- slist.list2[[1]]
	}
	Krange      <- 1:kmax
	Krange.plot <- 2:kmax
	## Two columns, longitude then latitude. Optional rownames as individual names.
	coords           <- read.table(coords)
	colnames(coords) <- c("Lon","Lat")
	x.min    <- min((coords[,1]-0.5))
	x.max    <- max((coords[,1]+0.5))
	y.min    <- min((coords[,2]-0.5))
	y.max    <- max((coords[,2]+0.5))
	world_sf <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
	world_sp <- rnaturalearth::ne_countries(scale=10,returnclass="sp")
	# empty list to hold admixture maps
	mapplot  <- list(); length(mapplot)     <- length(Krange.plot)
	for(K in 2:kmax){
		i=(K-1)
		if(K <= 15){
			myCols          <- goodcolors2(n=15)[1:K]
		}
		if(K>15){
			myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
		}
		q.matrix     <- slist2[[K]]
		my.palette   <- tess3r::CreatePalette(myCols, 9)
		mapplot.i    <- suppressWarnings(tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10), resolution = c(300,300), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max), background=TRUE, map.polygon=world_sp))
		mapplot2.i   <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA) + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
		mapplot[[i]] <- mapplot2.i + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
	}
	pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=mapplot, FUN=print)
	dev.off()
	mapplot
}

