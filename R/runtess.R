#' @title Run tess3r and plot results
#'
#' Pipeline with wrappers for several functions from the tess3r package.
#' Data can be supplied as a VCF file and results are plotted in multiple useful ways.
#' 
#' @param x 'vcfR' object (see package::vcfR) or character string with path to SNPs dataset fornatted to match the 'format' argument. Currently only VCF can be used, but others may be added soon.
#' @param format Character string indicating the format of the data. Currently only "VCF" or "fastStructure" allowed. Other types may be added. Ignored if x is a vcfR object.
#' @param coords Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.
#' @param samplenames NULL or a character string vector with names of samples (in same order) as data in x and coords. If NULL (the default), sample names are extracted from x.
#' @param kmax Numerical vector with set of values to use for K. Default 40.
#' @param reps Number of repititions. Default 100.
#' @param save.in Where to save output files. Default is NULL.
#' @param include.out Character vector indicating which type of files should be included as output. Default is c(".pdf",".Qlog", and ".entropyLog").
#' @param overwrite Logical indicating whether or not to allow new output files to overwrite existing ones. Default FALSE.
#' @param mask Proportion of input data to mask during each replicate when tess3 function is called. Default 0.05.
#' @param max.iteration Max iterations. Default 500.
#' @return List of plots
#' @export runtess
runtess <- function(x, format="VCF", coords, samplenames=NULL, kmax=10, reps=30, save.in=NULL, mask=0.05, max.iteration=500, include.out=c(".pdf",".Qlog",".entropyLog"), overwrite=FALSE){
	#if(is.null(save.as)){
	#	save.as <- file.path(getwd(),"result_tess.pdf")
	#}
	#if(file.exists(save.as)){
	#	stop("Output file already exists. Use a different name for 'save.as' argument.")
	#}
	

	if(is.null(save.in)){
		save.in <- getwd()
	}
	
	if(!is.null(include.out)){
		if(".pdf" %in% include.out){
			save.as.pdf <- file.path(save.in,"tess.pdf")
		} else {
			save.as.pdf <- NULL
		}
		if(".Qlog" %in% include.out){
			save.as.Qlog <- file.path(save.in,"tess.Qlog")
		} else {
			save.as.Qlog <- NULL
		}
		if(".entropyLog" %in% include.out){
			save.as.entropyLog <- file.path(save.in,"tess.entropyLog")
		} else {
			save.as.entropyLog <- NULL
		}
	} else {
		save.as.pdf <- save.as.Qlog <- save.as.entropyLog  <- NULL
	}
	if(!overwrite){
		files.to.check <- c(save.as.pdf,save.as.Qlog,save.as.entropyLog)
		if(!is.null(files.to.check)){
			if(any(files.to.check %in% save.in)){
				stop("One or more output files already exist in directory indicated by 'save.in'. Choose a different output directory or change 'overwrite' to TRUE")
			}
		}
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
		ploidy      <- length(unlist(strsplit(gt.mat[1],split="[/,|]")))
		if(is.null(samplenames)){
			samplenames <- colnames(vcf.obj@gt)[-1]
		}
		lfmm.obj    <- vcfR2lfmm(vcf=vcf)
	} else {
		stop("Currently, 'format' must be 'VCF'")
	}
	numind      <- length(samplenames)
	label.size  <- min((288/numind),7)
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
	if(max(Krange) > maxK){
		Krange <- 1:maxK
	}
	kmax <- max(Krange)
	tess.obj                   <- tess3r::tess3(X = lfmm.obj, coord = as.matrix(coords), K=Krange, ploidy = ploidy, verbose=FALSE ,mask=mask, rep=reps, max.iteration=max.iteration,keep="all")
	crossentropy.mat           <- do.call(rbind,lapply(X=1:length(tess.obj),FUN=function(x){matrix(unlist(tess.obj[[x]]["crossentropy"]),nrow=1)}))
	rownames(crossentropy.mat) <- Krange
	colnames(crossentropy.mat) <- paste0("rep",1:reps)
	mean.entropy      <- apply(crossentropy.mat, MARGIN=1, FUN=mean, na.rm=TRUE)
	range.entropy.mat <- do.call(rbind,lapply(X=1:nrow(crossentropy.mat),FUN=function(x){range(crossentropy.mat[x,],na.rm=TRUE)}))
	if(any(diff(mean.entropy)>0)){
		bestK <- unname(which(diff(mean.entropy)>0)[1])
	} else {
		bestK <- unname(Krange[1])
	}
	crossentropy.df <- data.frame(crossentropy=unname(unlist(c(crossentropy.mat))),K=rep(Krange,reps),replicate=rep(1:reps, each=kmax))
	if(!is.null(include.out)){
		if(".entropyLog" %in% include.out){
			write.table(x=crossentropy.df, file=save.as.entropyLog, row.names=T, col.names=T, quote=F, sep="\t")
		}
	}
	## List holding population assignment probabilities for each K
	slist <- lapply(X=Krange, FUN=function(x){as.data.frame(tess3r::qmatrix(tess3=tess.obj, K = x))})
	if(".Qlog" %in% include.out){
		q.df <- NULL
		for(K in 2:kmax){
			i=(K-1)
			q.matrix           <- slist[[K]]
			rownames(q.matrix) <- samplenames
			colnames(q.matrix) <- paste0("cluster", 1:ncol(q.matrix))
			q.i.df <- data.frame(individual=rep(rownames(q.matrix),ncol(q.matrix)), cluster=as.numeric(gsub("^cluster","",rep(colnames(q.matrix), each=nrow(q.matrix)))), assignment=c(unlist(unname(q.matrix))), K=rep(K,(K*numind)))
			q.df   <- rbind(q.df, q.i.df)
		}
		q.df[,"replicate"] <- NA
		write.table(x=q.df, file=save.as.Qlog, col.names=T, row.names=FALSE, quote=FALSE, sep="\t")
	}
	if(".pdf" %in% include.out){
		crossentropy.df$K <- factor(crossentropy.df$K, levels=c(1:nrow(crossentropy.df)))
		entropyPlot    <- ggplot2::ggplot(crossentropy.df, ggplot2::aes(x=K, y=crossentropy)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("Cross-entropy (",reps," replicates) vs. number of ancestral populations (K)"), x="Number of ancestral populations", y = "Cross-entropy") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) #+ ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
		coords.path    <- file.path(save.in,"coords.txt")
		write.table(x=coords,file=coords.path,sep="\t",quote=F)
		if(".Qlog" %in% include.out){
			admixturePlot  <- admixturePlots(xdir=save.in,userun=1:reps)
			assignmentPlot <- assignmentPlots(xdir=save.in,userun=1:reps)
			mapplot        <- admixtureMap(xdir=save.in,coords=coords.path,userun=1:reps)
			result         <- c(list(entropyPlot), admixturePlot, assignmentPlot, mapplot)
			return(result)
		} else {
			### Do things as before...
			Krange.plot    <- setdiff(Krange,1)
			admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
			assignmentPlot <- list(); length(assignmentPlot)  <- length(Krange.plot)
			mapplot        <- list(); length(mapplot)         <- length(Krange.plot)
			x.min <- min((coords[,1]-0.5))
			x.max <- max((coords[,1]+0.5))
			y.min <- min((coords[,2]-0.5))
			y.max <- max((coords[,2]+0.5))
			world_sf      <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
			world_sp      <- rnaturalearth::ne_countries(scale=10,returnclass="sp")
#			current_sf    <- sf::st_crop(world_sf,xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max)
#			current.gg.sf <- ggplot2::geom_sf(data=current_sf,colour = "black", fill = NA)
			for(K in Krange.plot){
				i=(K-1)
				q.matrix           <- slist[[K]]
				rownames(q.matrix) <- samplenames
				colnames(q.matrix) <- paste0("cluster", 1:ncol(q.matrix))
				indv.pop           <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))})
				posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
				posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
				if(K <= 15){
					myCols         <- goodcolors2(n=15)[1:K]
				}
				if(K>15){
					myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
				}
				# admixture barplots
				posterior.gg         <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Admixture Proportion",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
				admixturePlot[[i]]   <- posterior.gg
				# Best assignment for each individual
				indv.maxPosterior    <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
				labels               <- rep("",nrow(posterior.df))
				labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
				assignmentPlot[[i]]  <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv, fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
				### spatially interpolated admixture
				my.palette    <- tess3r::CreatePalette(myCols, 9)
				mapplot.i     <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE,map.polygon=world_sp)
				mapplot2.i    <- mapplot.i + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA) + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
				mapplot[[i]]  <- mapplot2.i + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			}
			if(".pdf" %in% include.out){
				pdf(height=6,width=10,file=save.as.pdf,onefile=TRUE)
					lapply(X=result,FUN=print)
				dev.off()
			}
			result <- c(list(entropyPlot), admixturePlot, assignmentPlot, mapplot)
			return(result)
		}
	} else {
		return(NULL)
	}
}
#' @examples
#' library(misc.wrappers)
#' 
#' # Path to VCF with SNPs
#' vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4.vcf.gz")
#' # Path to file with longitude and latitude of sampling locality of each individual
#' coords.path <- file.path(system.file("extdata", package = "misc.wrappers"), "simK4_coords.txt")
#' # Run tess3r 30 times each for K=1-10
#' runtess(x=vcf.path,coords=coords.path,kmax=10,reps=30,save.as="tess3r_simK4.pdf")


#' @title Filter (query) genotype by variation
#' 
#' Functon to test that a column of the genotypic matrix is a snp. Using default values for min.freqs this function tests if a site is a SNP.
#' This function is used internally in several other functions.
#' 
#' @param alignment.site Character vector containing alleles for each individual at a site.
#' @param min.freqs Numerical vector with minimum frequency of the most common allele, minimum frequency of the second most common allele, ect.. The length of this vector determines how many individuals must have non-missing data.
#' @return TRUE if site meets criteria specified by min.freqs, otherwise false.
#' @export filter.var
filter.var <- function(alignment.site,min.freqs=c(2,1,0,0)){
	## Number of each type of character in a given column
	allele.counts      <- table(tolower(alignment.site))
	#ACGT.counts        <- allele.counts[(names(allele.counts) %in% c("a","c","g","t"))]
	allele.counts.sorted <- sort(allele.counts,decreasing=T)
	### If the number of alleles (i.e., length of allele.counts) is more than the length of min.freqs,
	### go to the else argument and return FALSE because the site is too variable.
	if(length(allele.counts) <= length(min.freqs)){
		### If min.freqs only contains zeros, keep the site because even invariant sites are allowed.
		if(all(min.freqs==0)){
			result <- TRUE
		} else {
			### If the number of alleles is equal or greater than the number required by min.freqs (i.e., the number of nonzeros entries), continue, otherwise return FALSE.
			if(length(allele.counts) >= length(min.freqs[min.freqs>0])){
				### If min.freqs subtracted from allele.counts.sorted is a vector containing only zeros and/or positive numbers, then return TRUE, otherwise return FALSE.
				if(all((allele.counts.sorted[min.freqs>0]-min.freqs[min.freqs>0])>=0)){
					result <- TRUE
				} else {
					result <- FALSE
				}
			} else {
				result <- FALSE
			}
		} 
	} else {
		result <- FALSE
	}
	result
}

#' @title Reads a VCF object (including those not usually compatible with LEA) and returns a genotypic matrix equivalent to the lfmm-format used by LEA
#' 
#' The object returned by vcfR2lfmm can be used as input in the LEA function tess3.
#' Optionally supply a character string with path where lfmm object will be saved
#' This function is used in the function runtess
#' 
#' @param vcf Character string with path to input VCF file.
#' @param out Character string with path where output lfmm file should be saved. Default is NULL.
#' @return Matrix with genotypes in lfmm format
#' @export vcfR2lfmm
vcfR2lfmm <- function(vcf,out=NULL){
	vcf.obj   <- vcfR::read.vcfR(vcf,verbose=FALSE)
	gt.mat    <- gsub(":.+","",vcf.obj@gt[,-1])
	mat.temp0 <- gsub("|","/",gt.mat,fixed=TRUE)
	mat.temp1 <- gsub("^0/0$","0",mat.temp0)
	mat.temp2 <- gsub("^1/1$","2",mat.temp1)
	mat.temp3 <- gsub("^0/1$","1",mat.temp2)
	mat.temp4 <- gsub("./.",NA,mat.temp3,fixed=T)
	mat.temp5 <- unname(t(mat.temp4))
	mode(mat.temp5)     <- "numeric"
	colnames(mat.temp5) <- paste0("V",1:ncol(mat.temp5))
	if(length(unique(c(mat.temp5)))>4){
		stop("Some sites with >2 alleles")
	}
	if(!is.null(out)){
		write.table(mat.temp5,file=out,col.names=FALSE,row.names=FALSE,quote=FALSE,sep=" ")
	}
	mat.temp5
}







