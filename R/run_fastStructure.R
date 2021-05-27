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
#' @param save.as Where to save the output PDF. Default is NULL.
#' @param tolerance Tolerance for convergence, i.e., the change in marginal likelihood required to continue.
#' @param prior Type of prior to use. Default "simple"
#' @param full Whether or not to generate output files holding variation of Q, P, and marginal likelihood, in addition to the files holding means. Default FALSE.
#' @param seed Value to use as a seed for reproducing results. Default NULL.
#' @param python.path Character string with path to python 2 with fastStructure dependencies Numpy, Scipy, Cython, GNU Scientific Library
#' @param fastStructure.path Character string with path to folder containing the fastStructure python executable called 'structure.py'
#' @param cleanup Whether or not the original fastStructure output files (*.log, *.meanQ, *meanP file for each replicate of each K) should be deleted after the data from those files are compiled and saved in three tables. Default TRUE.
#' @param include.out Character vector indicating which type of files should be included as output. Default is c(".pdf",".Qlog",".margLlog"). An additional file ".Plog" can be included but can be very large.
#' @return List of plots
#' @export run_fastStructure
run_fastStructure <- function(x,format="VCF",coords=NULL,samplenames=NULL,kmax=40,reps=100,save.as=NULL,tolerance=10e-6,prior="simple",full=FALSE,seed=NULL,python.path=NULL,fastStructure.path=NULL,cleanup=TRUE,include.out=c(".pdf",".Qlog",".margLlog")){
	if(is.null(save.as)){
		save.as <- file.path(getwd(),"result_fastStructure.pdf")
	}
	if(file.exists(save.as)){
		stop("Output file already exists. Use a different name for 'save.as' argument.")
	}
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
	outdir.temp <- file.path(dirname(save.as),paste(sample(c(letters,LETTERS,rep(0:9,3)),10,replace=T),collapse=""))
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
	if(any(diff(mean.margL)<0)){
		bestK <- unname(which(diff(mean.margL)<0)[1])
	} else {
		# bestK <- unname(Krange[1])
		### might change this back later
		bestK <- unname(which(mean.margL==max(mean.margL))) 
	}
	margL.df      <- data.frame(marginalLikelihood=unname(unlist(c(margL.mat))),Kval=rep(Krange,reps))
	margL.df$Kval <- factor(margL.df$Kval, levels=c(1:nrow(margL.df)))
	margLPlot     <- ggplot2::ggplot(margL.df, ggplot2::aes(x=Kval, y=marginalLikelihood)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("marginal likelihood (",reps," replicates) vs. number of ancestral populations (K)"), x="Number of ancestral populations", y = "marginalLikelihood") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
	
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
	slist          <- qmats.list.best
	Krange.plot    <- setdiff(Krange,1)
	admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
	assignmentPlot <- list(); length(assignmentPlot)  <- length(Krange.plot)
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
	for(K in Krange.plot){
		i=(K-1)
		q.matrix  <- slist[[K]]
		rownames(q.matrix) <- samplenames
		colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
		indv.pop            <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))})
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
		if(K <= 15){
			myCols          <- goodcolors2(n=K)
		}
		if(K>15){
			myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
		}
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Admixture Proportion",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
		admixturePlot[[i]]  <- posterior.gg
		indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels             <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignment.K         <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
		assignmentPlot[[i]]  <- assignment.K
		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols, 9)
			mapplot.i       <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE, map.polygon=world_sp)
			mapplot2.i      <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA) + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
			mapplot[[i]]    <- mapplot2.i + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			#mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
		} else {
			mapplot[[i]] <- NULL
		}
	}
	result <- c(list(margLPlot),admixturePlot,assignmentPlot,mapplot)
	#### Save the PDF
	if(".pdf" %in% include.out){
	#if(!is.null(save.as)){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	if(format=="VCF"){
		# Delete the temporary fastStructure file
		remove.str <- file.remove(str.path0)
	}
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
		write.table(x=log.df,file=paste0(tools::file_path_sans_ext(save.as),".margLlog"),row.names=F,col.names=T,quote=F,sep="\t")
	}
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
			write.table(x=p.df,file=paste0(tools::file_path_sans_ext(save.as),".Plog"),row.names=F,col.names=T,quote=F,sep="\t")
		}
	}
	# Delete the folder with .log, .meanQ, and .meanP files; the info from these are compiled in a single file.
	if(cleanup){
		system(paste0("rm -R ",outdir.temp))
	}
	result
}
#' @examples
#' library(misc.wrappers)
#' # Define path to input VCF file containing similated data for 500 SNPs from 50 individuals in three populations.
#' example_vcf_path <- file.path(system.file("extdata", package = "misc.wrappers"),"simulated_K4.vcf.gz")
#' # Perform sNMF analyses on the simulated dataset for for K=2–10 and 30 replicates.
#' test_fastStructure_K4 <- run_fastStructure(x=example_vcf_path, kmax=10, reps=30, save.as="fastStructure_example_K4.pdf",include.out=c(".pdf"))




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
#' @param out Path where output structure file should be written. The mainparams file is also written using the same name but with the extension '.params'.
#' @return A list with [[1]] Value of 'out', which contains the path to the output file, and [[2]] values that should be used for some of the parameters in the mainparams file.
#' @export vcf2structure
vcf2structure <- function(vcf, IndvNames=TRUE, OneRowPerIndv=TRUE, MarkerNames=TRUE, MissingData=c(-9), out, InterMarkerDists=FALSE, PopData=NULL, PopFlag=NULL, LocData=NULL, Phenotype=NULL, OtherData=NULL){
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
	if(OneRowPerIndv){
		### Transpose the genotype matrix
		t.gt.mat <- t(gt.mat)
		### Split genotype strings by VCF haplotype separators ("/" or "|") and then reorganize into matrix with n columns per individual for n-ploidy individuals.
		mat.temp0 <- unname(do.call(rbind,lapply(1:nrow(t.gt.mat), FUN=function(x){unlist(strsplit(t.gt.mat[x,],split="[/,|]"))})))
		### Replace VCF missing data value "." with value supplied by MissingData argument.
		mat.temp1 <- gsub(".",MissingData,mat.temp0,fixed=TRUE)
	} else {
		### split genotype strings by VCF haplotype separators ("/" or "|") and then reorganize into matrix with n columns per individual for n-ploidy individuals.
		mat.temp0 <- t(unname(do.call(rbind,lapply(1:nrow(gt.mat), FUN=function(x){unlist(strsplit(gt.mat[x,],split="[/,|]"))}))))
		### replace VCF missing data value "." with value supplied by MissingData argument
		mat.temp1 <- gsub(".",MissingData,mat.temp0,fixed=TRUE)
		#### at this point the structure genotype matrix (mat.temp1) is in the correct format but columns and rows are not named. Now the other optional information rows and columns need to be considered for inclusion on output.
		# Repeat each value for the non-genotype columns twice 
		IndvNames   <- rep(IndvNames,ploidy,each=T)
		PopData     <- rep(PopData,ploidy,each=T)
		PopFlag     <- rep(PopFlag,ploidy,each=T)
		LocData     <- rep(LocData,ploidy,each=T)
		Phenotype   <- rep(Phenotype,ploidy,each=T)
		OtherData   <- rep(OtherData,ploidy,each=T)
	}
	### Individual/Genotype matrix
	indv.gt.mat <- cbind(IndvNames,PopData,PopFlag,LocData,Phenotype,OtherData,mat.temp1)
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
	# For now just return a list holding the values of the parameters in mainparams
	mainparams <- matrix(data=NA,nrow=23,ncol=1)
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
#' vcf2structure(vcf="Ahaetulla-prasina_AllPops_BestSNP.vcf",out="example.str")


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

