#' @title Run run_fastStructure from SNP data in a VCF file and and plot results.
#'
#' This fucntion is a wrapper that enables running fastStructure with SNP data in either a VCF file or vcfR object and coordinates in file or a matrix or data frame object.
#' 
#' Notes for running fastStructure:
#'    Requires python 2
#'    fastStructure 'cv' (cross-validation) feature is not yet implemented here.
#' 
#' @param x 'vcfR' object (see package::vcfR) or a character string with path to a SNPs dataset formatted according to the 'format' argument. Currently VCF or 'fastStructure' (a type of STRUCTURE format) can be used.
#' @param format Character string indicating the format of the data. Currently only "VCF" or "fastStructure" allowed. Other types may be added. Ignored if x is a vcfR object.
#' @param coords Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.
#' @param samplenames NULL or a character string vector with names of samples in the input data, and coords file if supplied. If NULL (the default), sample names are extracted from the SNPs datafile.
#' @param kmax Numerical vector with set of values to use for K. Default 40.
#' @param reps Number of repititions. Default 100.
##' @param save.as Where to save the output PDF. Default is NULL.
#' @param save.in Character string with path to directory where output files should be saved.
#' @param tolerance Tolerance for convergence, i.e., the change in marginal likelihood required to continue.
#' @param prior Type of prior to use. Default "simple".
#' @param full Whether or not to generate output files holding variation of Q, P, and marginal likelihood, in addition to the files holding means. Default FALSE.
#' @param seed Value to use as a seed for reproducing results. Default NULL.
#' @param python.path Character string with path to python 2 with fastStructure dependencies Numpy, Scipy, Cython, GNU Scientific Library
#' @param fastStructure.path Character string with path to folder containing the fastStructure python executable called 'structure.py'
#' @param cleanup Whether or not the original fastStructure output files (*.log, *.meanQ, *meanP file for each replicate of each K) should be deleted after the data from those files are compiled and saved in three tables. Default TRUE.
#' @param include.out Character vector indicating which type of files should be included as output. Default is c(".pdf",".Qlog",".margLlog"). An additional file ".Plog" can be included but can be very large.
#' @param debug Logical indicating whether or not to print messages indicating the internal step of the function.
#' @param overwrite Logical indicating whether or not to allow new output files to overwrite existing ones. Default FALSE.
#' @return List of plots
#' @export run_fastStructure
run_fastStructure <- function(x,format="VCF",coords=NULL,samplenames=NULL,kmax=10,save.in=NULL,reps=30,tolerance=10e-6,prior="simple",full=FALSE,seed=NULL,python.path=NULL,fastStructure.path=NULL,cleanup=TRUE,include.out=c(".pdf",".Qlog",".margLlog"), debug=FALSE,overwrite=FALSE){
	#if(is.null(save.as)){
	#	save.as <- file.path(getwd(),"result_fastStructure.pdf")
	#}
	if(is.null(save.in)){
		save.in <- getwd()
	}
	
	if(!is.null(include.out)){
		if(".pdf" %in% include.out){
			save.as.pdf <- file.path(save.in,"result_fastStructure.pdf")
		} else {
			save.as.pdf <- NULL
		}
		if(".Qlog" %in% include.out){
			save.as.Qlog <- file.path(save.in,"result_fastStructure.Qlog")
		} else {
			save.as.Qlog <- NULL
		}
		if(".margLlog" %in% include.out){
			save.as.margLlog <- file.path(save.in,"result_fastStructure.margLlog")
		} else {
			save.as.margLlog <- NULL
		}
		if(".Plog" %in% include.out){
			save.as.Plog <- file.path(save.in,"result_fastStructure.Plog")
		} else {
			save.as.Plog <- NULL
		}
	} else {
		save.as.pdf <- save.as.Qlog <- save.as.margLlog <- save.as.Plog <- NULL
	}
	if(!overwrite){
		files.to.check <- c(save.as.pdf,save.as.Qlog,save.as.margLlog,save.as.Plog)
		if(!is.null(files.to.check)){
			if(any(files.to.check %in% save.in)){
				stop("One or more output files already exist in directory indicated by 'save.in'. Choose a different output directory or change 'overwrite' to TRUE")
			}
		}
	}
	#if(file.exists(save.as)){
	#	stop("One or more output files already exists. Use a different name for 'save.as' argument.")
	#}
	Krange=1:kmax
	if(format=="VCF" | is(x,"vcfR")){
		if(is(x,"vcfR")){
			vcf.obj <- vcf <- x
		} else {
			vcf <- x
			vcf.obj     <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
		}
		gt.mat      <- gsub(":.+","",vcf.obj@gt[,-1])
		# Detect ploidy from genotype matrix of vcf
		test.sample <- unlist(gt.mat)[!is.na(unlist(gt.mat))][1]
		ploidy      <- length(unlist(strsplit(gt.mat[1],split="[/,|]",fixed=T)))
		if(is.null(samplenames)){
			samplenames <- colnames(vcf.obj@gt)[-1]
		}
		### Generate fastStrusture file from vcfR object
		str.path0   <- vcf2fastStructure(vcf=vcf.obj)
		str.path    <- tools::file_path_sans_ext(str.path0)
	} else {
		if(format=="fastStructure"){
			str.path0   <- x
			str.path    <- tools::file_path_sans_ext(str.path0)
			} else {
				stop("Data must be supplied in VCF or fastStructure format")
			}
	}
	numind      <- length(samplenames)
	label.size  <- min((288/numind),7)
	if(debug) message("step 0")
	### Checking that python2 exists and is executable
	if(is.null(python.path)){
		python.testpath <- find.exe.path("python")
		if(python.testpath!=1){
			python.path <- python.testpath
		} else {
			stop("No valid python identified. Set 'python.path' to location of python2.")
		}
	}
	### Checking that the fastStructure program 'structure.py' exists and is executable
	if(is.null(fastStructure.path)){
		fastStructure.testpath <- find.exe.path("structure.py")
		if(fastStructure.testpath!=1){
			fastStructure.path <- fastStructure.testpath
		} else {
			stop("No valid path to fastStructure identified. Set 'fastStructure.path' to location of 'structure.py' executable file")
		}
	}
	if(!file.exists(fastStructure.path)){
		stop(paste0(fastStructure.path," does not exist"))
	} else {
		if(file.info(fastStructure.path)$isdir){
			structure.py.path <- file.path(fastStructure.path,"structure.py")
		} else {
			structure.py.path <- fastStructure.path
		}
	}
	if(!file.exists(structure.py.path)){
		stop(paste0(structure.py.path," does not exist"))
	} else {
		if(check.if.executable(structure.py.path)!=0){
			stop("fastStructure 'structure.py' not executable")
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
	###### Defining other settings for fastStructure
	if(full){
		full <- " --full"
	} else {
		full <- NULL
	}
	if(!is.null(seed)){
		seed <- paste0(" --seed=",seed)
	} else {
		seed <- NULL
	}
	if(debug) message("step 1")
	outdir.temp  <- file.path(save.in,paste(sample(c(letters,LETTERS,rep(0:9,3)),10,replace=T),collapse=""))
	#outdir.temp <- file.path(dirname(save.as),paste(sample(c(letters,LETTERS,rep(0:9,3)),10,replace=T),collapse=""))
	dir.create(outdir.temp)
	for(i in 1:reps){
		outfile.temp.i <- file.path(outdir.temp,paste0("rep",i))
		for(K in Krange){
			command1     <- paste0(python.path," ",structure.py.path," -K ",K," --input=",str.path," --tol=",tolerance," --prior=",prior,full,seed," --format=str --output=", outfile.temp.i)
			run.command1 <- system(command1)
		}
		#if(FALSE){
		#	command2     <- paste0(python.path," ",chooseK.py.path," --input=",outfile.temp.i)
		#	run.command2 <- system(command2,intern=TRUE)
		#}
	}
	if(debug) message("step 2")
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
	colnames(margL.mat) <- paste0("rep",1:reps)
	mean.margL      <- apply(margL.mat,MARGIN=1,FUN=mean,na.rm=TRUE)
	range.margL.mat <- do.call(rbind,lapply(X=1:nrow(margL.mat),FUN=function(x){range(margL.mat[x,],na.rm=TRUE)}))
	if(debug) message("step 3")
	margL.df      <- data.frame(marginalLikelihood=unname(unlist(c(margL.mat))),Kval=rep(Krange,reps))
	margL.df$Kval <- factor(margL.df$Kval, levels=c(1:nrow(margL.df)))
	margLPlot     <- ggplot2::ggplot(margL.df, ggplot2::aes(x=Kval, y=marginalLikelihood)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("marginal likelihood (",reps," replicates) vs. number of ancestral populations (K)"), x="Number of ancestral populations", y = "marginalLikelihood") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
	if(debug) message("step 4")
	#### Load the Q matrices
	qpaths          <- list.files(outdir.temp,pattern="^.+\\.meanQ$",full.names=TRUE)
	qmats.list      <- lapply(qpaths,read.table)
	qfiles          <- basename(qpaths)
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
			myCols          <- goodcolors2(n=15)[1:K]
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
		pdf(height=6,width=10,file=save.as.pdf,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	if(debug) message("step 10")
	if(format=="VCF"){
		# Delete the temporary fastStructure file
		if(file.exists(str.path0)){
			remove.str <- file.remove(str.path0)
		}
	}
	if(debug) message("step 11")
	#### Generate compiled output files so that the other output files can be deleted.
	if(".margLlog" %in% include.out){
		### Organizing a single data frame to hold all data from all "*.log" files generated by fastStructure. This data frame will be written to a single file to replace the many '*.log' files.
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
		write.table(x=log.df, file=save.as.margLlog, row.names=F, col.names=T, quote=F, sep="\t")
	}
	if(debug) message("step 12")
	if(".Qlog" %in% include.out){
		### Organizing a single data frame to hold all q matrices for all replicates of every K. This data frame will be written to a single file so that the many '*.meanQ' files produced by fastStrucure can be deleted.
		q.df <- NULL
		for(i in 1:length(KQmats)){
			qmatrix.i    <- qmats.list[[i]]
			rownames(qmatrix.i) <- samplenames
			colnames(qmatrix.i) <- paste0("cluster",1:ncol(qmatrix.i))
			q.i.df <- data.frame(individual=rep(rownames(qmatrix.i),ncol(qmatrix.i)), cluster=as.numeric(gsub("^cluster","",rep(colnames(qmatrix.i),each=nrow(qmatrix.i)))), assignment=c(unlist(unname(qmatrix.i))), K=rep(KQmats[i],(KQmats[i]*numind)), replicate=rep(repsQmats[i], (KQmats[i]*numind)))
			# posterior.df <- data.frame(indv=rep(rownames(qmatrix.i),ncol(qmatrix.i)), pop=rep(colnames(qmatrix.i),each=nrow(qmatrix.i)), assignment=c(unlist(unname(qmatrix.i))),replicate=rep(repsQmats[i],(KQmats[i]*numind)))
			q.df <- rbind(q.df,q.i.df)
		}
		write.table(x=q.df,file=save.as.Qlog,row.names=F,col.names=T,quote=F,sep="\t")
	}
	if(debug) message("step 13")
	if(".Plog" %in% include.out){
		### Organizing a single data frame to hold all data from all "*.meanP" files generated by fastStructure. This data frame will be written to a single file to replace the many '*.meanP' files.
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
			write.table(x=p.df,file=save.as.Plog,row.names=F,col.names=T,quote=F,sep="\t")
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
#' # Run fastStructure 30 times each for K=1-10
#' run_fastStructure(x=vcf.path,kmax=10,reps=30,save.as="fs_simK4.pdf",include.out=c(".pdf"))
#' 
#' ## Example 2: Same SNP dataset as example 1, but here we also provide Lon/Lat coordinates of individuals to geographically interpolate admixture coefficients.
#' vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4.vcf.gz")
#' coords.path <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4_coords.txt")
#' run_fastStructure(x=vcf.path, coords=coords.path, kmax=10, reps=30, save.as="fs_simK4_withCoords.pdf", include.out=c(".pdf"))




#' @title Find path to program
#' 
#' This function searches the system path and the settings file of the misc.wrappers program for a usable path to a program with a specified name.
#' 
#' @param program Character string with name of program to find
#' @return Character string with path to the program's executable file, or 1 if path not found. 
#' @export find.exe.path
find.exe.path <- function(program){
	linked.paths <- misc.wrappers::config_miscwrappers()
	if(any(linked.paths$program == program)){
		trypath <- linked.paths$exe_path[which(linked.paths$program==program)[1]]
		if(misc.wrappers::check.if.executable(trypath)==0){
			return(trypath)
		}
	}
	check.paths <- system(paste("which",program),ignore.stdout=TRUE,ignore.stderr=TRUE)
	if(check.paths==0){
		check.paths <- system(paste("which",program),intern=TRUE)[1]
		if(misc.wrappers::check.if.executable(check.paths)==0){
			return(check.paths)
		} else {
			return(1)
		}
	} else {
		return(1)
	}
}


#' @title Convert VCF file/object to fastStructure input file (for SNP data)
#' 
#' Converts a VCF file or vcfR object with SNP data into a fastStructure file.
#' 
#' @param vcf Character string with path to input VCF, or an object of class vcfR.
#' @param IndvNames Logical indicating if the first column of the output file should include the names of individuals. Default TRUE.
#' @param OtherData NULL (the default) or a matrix with up to five (if IndvNames TRUE) or six (if IndvNames FALSE) columns of metadata. If non-NULL, the number of rows must equal number of individuals in the VCF.
#' @param out Path where output file should be written.
#' @return Character string with value of 'out'.
#' @export vcf2fastStructure
vcf2fastStructure <- function(vcf, IndvNames=TRUE, out=NULL, OtherData=NULL){
	if(is.null(out)){
		out  <- file.path(getwd(),paste0(paste(sample(c(letters,LETTERS,rep(0:9,3)),30,replace=T),collapse=""),".str"))
	}
	# Check class of vcf argument value.
	if(is(vcf,"vcfR")){
		vcf.obj <- vcf
	} else {
		if(is(vcf,"character")){
			# read VCF file into R as a vcfR object
			vcf.obj     <- vcfR::read.vcfR(vcf)
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
	# Number of individuals
	numind      <- length(samplenames)
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
	#### Metadata columns
	if(!is.null(OtherData)){
		if(is(OtherData,c("matrix","array","data.frame"))){
			if(nrow(OtherData)!=numind){
				stop("When 'OtherData' is non-NULL, it must be a matrix nrow = number of individuals")
			} else {
				if(!is.null(IndvNames)){
					ExtraCols <- ncol(OtherData) + 1
				} else {
					ExtraCols <- ncol(OtherData)
				}
				if(ExtraCols>6){
					stop("Only six columns of metadata allowed before genotype columns")
				} else {
					OtherData <- do.call(rbind,lapply(X=1:nrow(OtherData),FUN=function(x){paste(OtherData[x,],collapse=" ")}))
				}
			}
		} else {
			stop("When 'OtherData' is non-NULL, it must be a matrix nrow = number of individuals")
		}
	} else {
		if(!is.null(IndvNames)){
			OtherData <- matrix(data=NA,nrow=(numind*ploidy),ncol=5)
		} else {
			OtherData <- matrix(data=NA,nrow=(numind*ploidy),ncol=6)
		}
	}
	### split genotype strings by VCF haplotype separators ("/" or "|") and then reorganize into matrix with n columns per individual for n-ploidy individuals.
	mat.temp0 <- t(unname(do.call(rbind,lapply(1:nrow(gt.mat), FUN=function(x){unlist(strsplit(gt.mat[x,],split="[/,|]"))}))))
	### replace VCF missing data value (".") with fastStructure missing data value ("-9")
	mat.temp1 <- gsub(".","-9",mat.temp0,fixed=TRUE)
	### Merge metadata and genotype matrices
	result.mat <- cbind(IndvNames,OtherData,mat.temp1)
	write.table(x=result.mat,file=out,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
	### Function returns the path where the fastStructure file was saved
	out
}
#' @examples
#' vcf2fastStructure(vcf="Ahaetulla-prasina_AllPops_BestSNP.vcf",out="example.str")

