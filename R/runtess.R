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
	mat.temp1 <- gsub("^0/0$","0",gt.mat)
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

#' @title Run tess3r and plot results
#'
#' This is a wrapper for several functions from the tess3r package. The benefit of using this function is that data can be supplied as a VCF file, and results are plotted in multiple useful ways.
#' 
#' @param vcf Path to input VCF file with SNP data
#' @param coords Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.
#' @param Krange Numerical vector with set of values to use for K.
#' @param ploidy Number of copies of each chromosome. Default 2.
#' @param mask Proportion of input data to mask during each replicate when tess3 function is called
#' @param reps Number of repititions. Default 100.
#' @param max.iteration Max iterations. Default 500.
#' @param save.as Where to save the output PDF. Default is NULL. **This argument is ignored in some environments. Instead, use dev.new(file="Where/To/Save/Output.pdf",height=6,width=10,noRStudioGD=TRUE) before calling 'runtess'. Then dev.off().
#' @return List of plots
#' @export runtess
runtess <- function(vcf,coords=NULL,Krange=1:40,ploidy=2,mask=0.05,reps=100,max.iteration=500,save.as=NULL){
	vcf.obj     <- vcfR::read.vcfR(vcf)
	samplenames <- colnames(vcf.obj@gt)[-1]
	numind      <- length(samplenames)
	label.size  <- min((288/numind),7)
	lfmm.obj    <- vcfR2lfmm(vcf=vcf)
	if(is(coords,"array") | is(coords,"data.frame")){
		coords <-  coords
	} else {
		if(file.exists(coords)){
			coords   <- read.table(coords)
		}
	}
	maxK <- min(nrow(unique(coords)),(numind-1))
	if(max(Krange) > maxK){
		Krange <- 1:maxK
	}
	tess.obj <- tess3r::tess3(X = lfmm.obj, coord = as.matrix(coords), K=Krange, ploidy = ploidy, verbose=FALSE ,mask=mask, rep=reps, max.iteration=max.iteration,keep="all")
	crossentropy.mat <- do.call(rbind,lapply(X=1:length(tess.obj),FUN=function(x){matrix(unlist(tess.obj[[x]]["crossentropy"]),nrow=1)}))
	rownames(crossentropy.mat) <- Krange
	colnames(crossentropy.mat) <- paste0("rep",1:reps)
	par(mfrow=c(1,1))
	mean.entropy <- apply(crossentropy.mat,MARGIN=1,FUN=mean,na.rm=TRUE)
	range.entropy.mat <- do.call(rbind,lapply(X=1:nrow(crossentropy.mat),FUN=function(x){range(crossentropy.mat[x,],na.rm=TRUE)}))
	boxplot(t(crossentropy.mat))
	if(any(diff(mean.entropy)>0)){
		bestK <- unname(which(diff(mean.entropy)>0)[1])
	} else {
		bestK <- unname(Krange[1])
	}
	### Criteria 3: Which K value (for K>=2) yields the least variable entropy scores.
	Entropy.variation <- apply(X=crossentropy.mat,MARGIN=1,FUN=var,na.rm=TRUE)
	Kbest.criteria3   <- which(Entropy.variation==min(Entropy.variation[-1]))
	### Criteria 4: t-tests for entropy of each pairwise adjacent K
	for(i in 2:nrow(crossentropy.mat)){
		if(Entropy.variation[Kbest.criteria3]==0){
			Kbest.criteria4 <- NULL
			break
		}
		t.test.i <- t.test(crossentropy.mat[i-1,],crossentropy.mat[i,])
		pval.i   <- t.test.i$p.value
		stat.i   <- t.test.i$statistic
		if(pval.i < 0.05 & stat.i > 0){
			next
		} else {
			if(i==nrow(crossentropy.mat)){
				Kbest.criteria4 <- NULL
			} else{
				Kbest.criteria4 <- (i-1)
				break
			}
		}
	}
	###
	if(bestK>1){
		segments(x0=bestK,y0=par("usr")[3],y1=par("usr")[4],lty=2)
	}
	mtext(side=1,"Number of ancestral populations",line=2.2)
	mtext(side=2,"Cross-validation score",line=2.2)
	mtext(side=3,paste0("Cross-validation score (",reps," replicates) vs. number of ancestral populations (K)"),line=1)
	axis(1,at=Krange)
	entropyPlot <- recordPlot()
	## List holding population assignment probabilities for each K
	slist <- lapply(X=Krange,FUN=function(x){as.data.frame(tess3r::qmatrix(tess3=tess.obj, K = x))})
	par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
	Krange.plot    <- setdiff(Krange,1)
	admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
	mapplot        <- list(); length(mapplot)         <- length(Krange.plot)
	x.min <- min((coords[,1]-0.5))
	x.max <- max((coords[,1]+0.5))
	y.min <- min((coords[,2]-0.5))
	y.max <- max((coords[,2]+0.5))
	for(K in Krange.plot){
		i=(K-1)
		q.matrix  <- slist[[K]]
		rownames(q.matrix) <- samplenames
		colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		if(K < 5){
			myCols          <- goodcolors(K,thresh=100)
		}
		if(K >= 5 & K < 7){
			myCols          <- goodcolors(K,thresh=100,cbspace="deut")
		}
		if(K >= 7 & K < 15){
			myCols          <- goodcolors(K,thresh=100,cbspace="")
		}
		if(K>=15){
			myCols          <- c(goodcolors(14,thresh=100,cbspace=""), sample(adegenet::funky(100), size=K-14))
		}
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Membership Probability",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
		plot(posterior.gg)
		admixturePlot[[i]]   <- recordPlot()
		my.palette      <- tess3r::CreatePalette(myCols, 9)
		xdist           <- geosphere::distm(x=c(x.min,0),y=c(x.max,0))
		ydist           <- geosphere::distm(x=c(0,y.min),y=c(0,y.max))
		mapplot.initial <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), main = "", xlab = "", ylab = "",resolution = c(2,2), col.palette = lapply(X=1:K,FUN=function(x){rep("#FFFFFF",9)}), cex=0,window=c(x.min,x.max,y.min,y.max),asp=xdist/ydist,add=FALSE)
		mapplot.i       <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(500,500), cex = 0.4, col.palette = my.palette, window=par("usr"),asp=xdist/ydist,add=FALSE)
		maps::map(add=TRUE)
#		mapplot.initial <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(100,100), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE)
#		plot(mapplot.initial + ggplot2::theme_classic())
#		map(xlim=c(x.min,x.max),ylim=c(y.min,y.max),mar=c(5.1,4.1,4.1,2.1))
#		mapplot.i       <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=par("usr"),background=TRUE)
#		plot(mapplot.i + ggplot2::theme_classic())
		mapplot[[i]]    <- recordPlot()
	}
	result <- c(list(entropyPlot),admixturePlot,mapplot)
	if(!is.null(save.as)){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	result
}
#' @examples
####
#' library(maps)
#' library(tess3r)
#' libary(ade)
#' library(adegenet)
#' library(ggplot2)
#' library(rworldmap)
#' library(LEA)
#' library(vcfR)
#' library(geosphere)
#' #source("~/DAPC_adegenet.R")
#' #source("~/runtess.R")
#' dev.new(width=10,height=6)
#' Oxyrhabdium_AllSpecies <- runtess(vcf="Oxyrhabdium_AllSpecies_BestSNP.vcf",
#'                      coords="Oxyrhabdium_AllSpecies_coords.txt",
#'                      Krange=1:15,
#'                      max.iteration=500,
#'                      save.as="Oxyrhabdium_AllSpecies_BestSNP_tess3r_run2.pdf")
#'
#' leporinum <- runtess(vcf="Oxyrhabdium-leporinum_BestSNP.vcf",
#'                      coords="Oxyrhabdium-leporinum_coords.txt",
#'                      Krange=1:15,
#'                      max.iteration=500,
#'                      save.as="Oxyrhabdium-leporinum_BestSNP_tess3r_run3.pdf")
#'
#' modestum        <- runtess(vcf="Oxyrhabdium-modestum_BestSNP.vcf",
#'                            coords="Oxyrhabdium-modestum_coords.txt",
#'                            Krange=1:15,
#'                            max.iteration=500,
#'                            save.as="Oxyrhabdium-modestum_BestSNP_tess3r_run3.pdf")
#'
#' cfmodestum        <- runtess(vcf="Oxyrhabdium-cf.modestum_BestSNP.vcf",
#'                            coords="Oxyrhabdium-cfmodestum_coords.txt",
#'                            Krange=1:15,
#'                            max.iteration=500,
#'                            save.as="Oxyrhabdium-cf.modestum_BestSNP_tess3r_run2.pdf")
#'
#' leporinum_Luzon <- runtess(vcf="Oxyrhabdium-leporinum_Luzon_BestSNP.vcf",
#'                            coords="Oxyrhabdium-leporinum_Luzon_coords.txt",
#'                            Krange=1:15,
#'                            max.iteration=500,
#'                            save.as="Oxyrhabdium-leporinum_Luzon_BestSNP_tess3r_run2.pdf")
#'
#' cfmodestum_Luzon        <- runtess(vcf="Oxyrhabdium-cf.modestum_Luzon_BestSNP.vcf",
#'                            coords="Oxyrhabdium-cfmodestum_Luzon_coords.txt",
#'                            Krange=1:15,
#'                            max.iteration=500,
#'                            save.as="Oxyrhabdium-cf.modestum_Luzon_BestSNP_tess3r_run3.pdf")








