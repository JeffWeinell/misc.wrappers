#' @title Wrapper for LEA snmf
#' 
#' Runs LEA function snmf function and generate a pdf with various useful plots; optionally interpolate cluster assignments onto a map.
#' 
#' @param x 'vcfR' object (see package::vcfR) or a character string with path to a SNPs dataset formatted according to the 'format' argument. Currently VCF or 'fastStructure' (a type of STRUCTURE format) can be used.
#' @param format Character string indicating the format of the data. Currently only "VCF" or "fastStructure" allowed. Other types may be added. Ignored if x is a vcfR object.
#' @param coords Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.
#' @param samplenames NULL or a character string vector with names of samples (in same order) as data in x and coords. If NULL (the default), sample names are extracted from x.
#' @param kmax Number indicating the maximum number of clusters to evaluate. Default is 10, which is converted using kmax = min(kmax, number of individuals-1)
#' @param reps Number of repititions. Default 30.
#' @param entropy Default TRUE
#' @param project Defualt 'new'
#' @param iter Default 500
##' @param save.as Character string with where to save the output PDF with plots of results. Default is NULL.
#' @param save.in Character string with path to directory where output files should be saved.
#' @param include.out Character strings indicating what output files should be generated. If "entropy.mat" is provided, the function writes the crossentropy matrix and is done.
#' @param overwrite Logical indicating whether or not to allow new output files to overwrite existing ones. Default FALSE.
#' @param ... Additional arguments passed to 'snmf' function of LEA package
#' @return List of plots
#' @export run_sNMF
run_sNMF <- function(x,format="VCF",coords=NULL, samplenames=NULL, kmax=10, reps=30, entropy=TRUE, project="new", iter=500, save.in=NULL, include.out=c(".pdf",".Qlog",".entropyLog"), overwrite=FALSE, ...){
	if(is.null(save.in)){
		save.in <- getwd()
	}
	if(!is.null(include.out)){
		if(".pdf" %in% include.out){
			save.as.pdf <- file.path(save.in,"result_snmf.pdf")
		} else {
			save.as.pdf <- NULL
		}
		if(".Qlog" %in% include.out){
			save.as.Qlog <- file.path(save.in,"result_snmf.Qlog")
		} else {
			save.as.Qlog <- NULL
		}
		if(".entropyLog" %in% include.out){
			save.as.entropyLog <- file.path(save.in,"result_snmf.entropyLog")
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
			#vcf.obj     <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
			vcf.obj     <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F,convertNA=FALSE)
		}
		if(is.null(samplenames)){
			samplenames <- colnames(vcf.obj@gt)[-1]
		}
	} else {
		stop("Currently, 'format' must be 'VCF'")
	}
	numind      <- length(samplenames)
	label.size  <- min((288/numind),7)
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
		x.min <- min((coords[,1]-0.5))
		x.max <- max((coords[,1]+0.5))
		y.min <- min((coords[,2]-0.5))
		y.max <- max((coords[,2]+0.5))
		maxK  <- min(nrow(unique(coords)),(numind-1))
		world_sf      <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
		world_sp      <- rnaturalearth::ne_countries(scale=10,returnclass="sp")
	} else {
		maxK <- (numind-1)
	}
	if(max(Krange) > maxK){
		Krange <- 1:maxK
	}
	kmax <- max(Krange)
	geno.temp.path   <- paste0(tempfile(),".geno")
	geno.obj         <- vcfR2geno(vcf=vcf, out=geno.temp.path)
	num.CPUs         <- parallel::detectCores()
	if(!is.numeric(num.CPUs)){
		num.CPUs <- 2
	}
	snmf.obj         <- LEA::snmf(geno.temp.path, K=Krange, repetitions=reps, entropy=entropy, project="new", iterations=iter, CPU=num.CPUs,...)
	crossentropy.mat <- t(do.call(cbind, lapply(X=Krange, FUN=function(x){LEA::cross.entropy(snmf.obj, K = x)})))
	rownames(crossentropy.mat) <- Krange
	colnames(crossentropy.mat) <- paste0("rep", 1:reps)
	crossentropy.df <- data.frame(crossentropy=unname(unlist(c(crossentropy.mat))), K=rep(Krange,reps), replicate=rep(1:reps, each=kmax))
	if(!is.null(include.out)){
		if(".entropyLog" %in% include.out){
			write.table(x=crossentropy.df, file=save.as.entropyLog, row.names=F, col.names=T, quote=F, sep="\t")
		}
	}
	mean.entropy   <- apply(crossentropy.mat, MARGIN=1, FUN=mean, na.rm=TRUE)
	if(any(diff(mean.entropy)>0)){
		bestK <- unname(which(diff(mean.entropy)>0)[1])
	} else {
		bestK <- unname(Krange[1])
	}
	Kbest.criteria1 <- bestK
	if(".Qlog" %in% include.out){
		#qlist <- list(); length(qlist) <- length(2:Krange)
		q.df <- NULL
		for(K in 2:kmax){
			i=(K-1)
			ce                 <- LEA::cross.entropy(snmf.obj, K = K)
			best               <- which.min(ce)
			q.matrix           <- LEA::Q(snmf.obj, K=K, run=best)
			rownames(q.matrix) <- samplenames
			colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
			q.i.df <- data.frame(individual=rep(rownames(q.matrix),ncol(q.matrix)), cluster=as.numeric(gsub("^cluster","",rep(colnames(q.matrix), each=nrow(q.matrix)))), assignment=c(unlist(unname(q.matrix))), K=rep(K,(K*numind)))
			q.df   <- rbind(q.df,q.i.df)
		}
		q.df[,"replicate"] <- NA
		write.table(x=q.df, file=save.as.Qlog, col.names=T, row.names=FALSE, quote=FALSE, sep="\t")
	}
	if(".pdf" %in% include.out){
		crossentropy.df$K <- factor(crossentropy.df$K, levels=c(1:nrow(crossentropy.df)))
		entropyPlot <- ggplot2::ggplot(crossentropy.df, ggplot2::aes(x=K, y=crossentropy)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("Cross-entropy (",reps," replicates) vs. number of ancestral populations (K)"), x="Number of ancestral populations", y = "Cross-entropy") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
		Krange.plot    <- setdiff(Krange,1)
		admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
		assignmentPlot <- list(); length(assignmentPlot)  <- length(Krange.plot)
		mapplot        <- list(); length(mapplot)         <- length(Krange.plot)
		for(K in Krange.plot){
			i=(K-1)
			ce           <- LEA::cross.entropy(snmf.obj, K = K)
			best         <- which.min(ce)
			if(K <= 15){
				myCols          <- goodcolors2(n=15)[1:K]
			}
			if(K>15){
				myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=K-15))
			}
			q.matrix           <- LEA::Q(snmf.obj,K=K,run=best)
			rownames(q.matrix) <- samplenames
			colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
			indv.pop           <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))})
			posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
			posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
			posterior.gg       <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Admixture Proportion",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
			admixturePlot[[i]]   <- posterior.gg
			indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
			labels             <- rep("",nrow(posterior.df))
			labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
			assignment.K       <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
			assignmentPlot[[i]]  <- assignment.K
			if(!is.null(coords)){
				my.palette      <- tess3r::CreatePalette(myCols, 9)
				mapplot.i     <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE,map.polygon=world_sp)
				mapplot2.i    <- mapplot.i + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA) + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
				mapplot[[i]]  <- mapplot2.i + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			}
		}
		if(!is.null(coords)){
			result <- c(list(entropyPlot),admixturePlot,assignmentPlot,mapplot)
		} else {
			result <- c(list(entropyPlot),admixturePlot,assignmentPlot)
		}
		if(!is.null(save.as.pdf)){
			pdf(height=6,width=10,file=save.as.pdf,onefile=TRUE)
				lapply(X=result,FUN=print)
			dev.off()
		}
		result
		return(result)
	} else {
		return(NULL)
	}
}
#' @examples
#' library(misc.wrappers)
#' 
#' ## Example 1:
#' # Path to VCF with SNPs
#' vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
#' run_sNMF(x=vcf.path,format="VCF",kmax=10,reps=30,save.as="sNMF_simK4.pdf")
#' 
#' ## Example 2: Same SNP dataset as example 1, but here we also provide Lon/Lat coordinates of individuals to geographically interpolate admixture coefficients.
#' vcf.path    <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
#' coords.path <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4_coords.txt")
#' run_sNMF(x=vcf.path,format="VCF",coords=coords.path,samplenames=NULL,kmax=10,reps=30,entropy=TRUE,project="new",iter=500,save.as="sNMF_simK4_withCoords.pdf")


#' @title Convert VCF object/file to genotypic matrix
#' 
#' The object returned by vcfR2geno can be used as input in the LEA function snmf
#' Optionally supply a character string with path where geno object will be saved
#' This function is used in the function run_sNMFR
#' 
#' @param vcf Character string with path to input VCF file.
#' @param out Character string with path where output geno file should be saved. Default is NULL.
#' @return Matrix with genotypes in geno format
#' @export vcfR2geno
vcfR2geno <- function(vcf,out=NULL){
	#vcf.obj   <- vcfR::read.vcfR(vcf,verbose=FALSE,convertNA=FALSE)
	vcf.obj   <- vcfR::read.vcfR(vcf,verbose=FALSE,convertNA=FALSE)
	gt.mat    <- gsub(":.+","",vcf.obj@gt[,-1])
	mat.temp0 <- gsub("|","/",gt.mat,fixed=TRUE)
	mat.temp1 <- gsub("^0/0$","0",mat.temp0)
	mat.temp2 <- gsub("^1/1$","2",mat.temp1)
	mat.temp3 <- gsub("^0/1$","1",mat.temp2)
	mat.temp4 <- gsub("./.","9",mat.temp3,fixed=T)
	if(length(unique(c(mat.temp4)))>4){
		stop("Some sites with >2 alleles")
	}
	if(!is.null(out)){
		lines.temp <- apply(mat.temp4,MARGIN=1,FUN=paste,collapse="")
		writeLines(text=lines.temp,con=out)
	}
	mat.temp4
}



