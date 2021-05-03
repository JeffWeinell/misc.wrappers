#' Wrapper for LEA snmf
#' 
#' Runs LEA function snmf function from VCF input file and optionally interpolates cluster assignments onto a map (if 'coords' argument is non-NULL).
#' 
#' @param vcf Path to vcf input file
#' @param coords Default NULL
#' @param Krange Default 1:40
#' @param reps Default 100
#' @param entropy Default TRUE
#' @param project Defualt 'new'
#' @param iter Default 500
#' @param CPU Defualt 2
#' @param save.as Path where output PDF should be saved. Default NULL.
#' @return List of plots
#' @export 
run_SNMF <- function(vcf,coords=NULL,Krange=1:40,reps=100,entropy=TRUE,project="new",iter=500,CPU=2,save.as=NULL){
	vcf.obj     <- vcfR::read.vcfR(vcf)
	samplenames <- colnames(vcf.obj@gt)[-1]
	numind      <- length(samplenames)
	label.size  <- min((288/numind),7)
	
	if(!is.null(coords)){
		if(is(coords,"array") | is(coords,"data.frame")){
			coords <-  coords
		} else {
			if(file.exists(coords)){
				coords   <- read.table(coords)
			}
		}
		### Check that all individuals with coords are in the vcf file, and vice versa.
		if(!all(samplenames %in% rownames(coords) & rownames(coords) %in% samplenames)){
			stop("All individuals in coords file must be in vcf")
		}
		x.min <- min((coords[,1]-0.5))
		x.max <- max((coords[,1]+0.5))
		y.min <- min((coords[,2]-0.5))
		y.max <- max((coords[,2]+0.5))
		maxK <- min(nrow(unique(coords)),(numind-1))
	} else {
		maxK <- (numind-1)
	}
	if(max(Krange) > maxK){
		Krange <- 1:maxK
	}
	geno.temp.path <- paste0(tempfile(),".geno")
	geno.obj       <- vcfR2geno(vcf=vcf,out=geno.temp.path)
	snmf.obj       <- LEA::snmf(geno.temp.path,K=Krange, repetitions=reps,entropy=entropy,project="new",iterations=iter,CPU=CPU)
	crossentropy.mat  <- t(do.call(cbind,lapply(X=Krange,FUN=function(x){LEA::cross.entropy(snmf.obj,K = x)})))
	rownames(crossentropy.mat) <- Krange
	colnames(crossentropy.mat) <- paste0("rep",1:reps)
	boxplot(t(crossentropy.mat))
	mean.entropy   <- apply(crossentropy.mat,MARGIN=1,FUN=mean,na.rm=TRUE)
#	plot(Krange,mean.entropy,pch=21,col="blue",xlab="",ylab="",xlim=range(Krange), ylim=range(range.entropy.mat))
#	arrows(x0=Krange,y0=range.entropy.mat[,1],y1=range.entropy.mat[,2],length=0.07,col="black",angle=90,code=3)
#	lines(Krange,mean.entropy,col="blue")
	if(any(diff(mean.entropy)>0)){
		bestK <- unname(which(diff(mean.entropy)>0)[1])
	} else {
		bestK <- unname(Krange[1])
	}
	Kbest.criteria1   <- bestK
	### Criteria 2: Checking if the max BIC of the next K value is less than the min BIC of the previous K value.
	range.entropy.mat <- do.call(rbind,lapply(X=1:nrow(crossentropy.mat),FUN=function(x){range(crossentropy.mat[x,],na.rm=TRUE)}))
	entropy.is.nonoverlap <- NULL
	entropy.is.reduced    <- NULL
	for(i in 2:nrow(range.entropy.mat)){
		entropy.is.nonoverlap  <- c(entropy.is.nonoverlap, range.entropy.mat[i,1] > range.entropy.mat[(i-1),2] | range.entropy.mat[i,2] < range.entropy.mat[(i-1),1])
		entropy.is.reduced     <- c(entropy.is.reduced, range.entropy.mat[i,2] < range.entropy.mat[(i-1),1])
	}
	# If max BIC for K=2 is better than BIC K=1 and if some K values are not better (all BIC lower) than K-1, then find the first K value in which K+1 is not better.
	if(any(!entropy.is.reduced) & entropy.is.reduced[1]){
		Kbest.criteria2 <- which(!entropy.is.reduced)[1]
	} else {
		Kbest.criteria2 <- 1
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
	if(bestK>1){
		segments(x0=bestK,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="black")
	}
	# segments(x0=bestK,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="green")
	# segments(x0=Kbest.criteria2,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="blue")
	# segments(x0=Kbest.criteria3,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="red")
	mtext(side=1,"Number of ancestral populations",line=2.2)
	mtext(side=2,"Cross-entropy",line=2.2)
	mtext(side=3,paste0("Cross-entropy (",reps," replicates) vs. number of ancestral populations (K)"),line=1)
	axis(1,at=Krange)
	entropyPlot    <- recordPlot()
	Krange.plot    <- setdiff(Krange,1)
	admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
	mapplot        <- list(); length(mapplot)         <- length(Krange.plot)
	for(K in Krange.plot){
		i=(K-1)
		ce           <- LEA::cross.entropy(snmf.obj, K = K)
		best         <- which.min(ce)
	#	q.matrix.best <- suppressWarnings(tess3r::as.qmatrix(LEA::Q(snmf.obj,K=bestK,run=best)))
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
		
		q.matrix           <- LEA::Q(snmf.obj,K=K,run=best)
		rownames(q.matrix) <- samplenames
		colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		posterior.gg       <- ggplot2::ggplot(posterior.df, aes(fill= pop, x= assignment, y=indv)) + geom_bar(position="stack", stat="identity") + theme_classic() + theme(axis.text.y = element_text(size = label.size), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + ggplot2::labs(x = "Membership Probability",y="",fill="Cluster",title=paste0("K = ",K)) + scale_fill_manual(values=myCols[1:K])
		plot(posterior.gg)
		admixturePlot[[i]]   <- recordPlot()
		
		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols, 9)
			xdist           <- geosphere::distm(x=c(x.min,0),y=c(x.max,0))
			ydist           <- geosphere::distm(x=c(0,y.min),y=c(0,y.max))
			mapplot.initial <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), method = "map.max", interpol = FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(100,100), cex = 0.4,col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),asp=xdist/ydist)
			mapplot.i       <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), method = "map.max", interpol = FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(500,500), cex = 0.4,col.palette = my.palette, window=c(par("usr")[1],par("usr")[2],par("usr")[3],par("usr")[4]),asp=xdist/ydist)
			mapplot[[i]]    <- recordPlot()
		}
	}
	if(!is.null(coords)){
		result <- c(list(entropyPlot),admixturePlot,mapplot)
	} else {
		result <- c(list(entropyPlot),admixturePlot)
	}
	if(!is.null(save.as)){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	result
}


####### Runs
#if(FALSE){
#	library(ade4)
#	library(adegenet)
#	library(vcfR)
#	library(ggplot2)
#	library(maps)
#	library(geosphere)
#	library(tess3r)
#	library(LEA)
#	library(rworldmap)
#	#library(dartR)
#	source("/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/DAPC_adegenet.R")
#	source("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/runtess.R")
#	source("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/LEA/run_SNMF.R")
#	####
#	dev.new(width=10,height=6)
#	leporinum <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium-leporinum_BestSNP.vcf",
#	                      coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/tess3r/Oxyrhabdium-leporinum_coords.txt",
#	                      save.as="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/LEA/Oxyrhabdium-leporinum_BestSNP_sNMF_run3.pdf")
#
#	modestum  <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium-modestum_BestSNP.vcf",
#	                     coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/tess3r/Oxyrhabdium-modestum_coords.txt",
#	                     save.as="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/LEA/Oxyrhabdium-modestum_BestSNP_sNMF_run1.pdf")
#
#	cfmodestum  <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium-cf.modestum_BestSNP.vcf",
#	                     coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/tess3r/Oxyrhabdium-cfmodestum_coords.txt",
#	                     save.as="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/LEA/Oxyrhabdium-cfmodestum_BestSNP_sNMF_run1.pdf")
#
#	cfmodestum_Luzon        <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium-cf.modestum_Luzon_BestSNP.vcf",
#	                           coords="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/tess3r/Oxyrhabdium-cfmodestum_Luzon_coords.txt",
#	                           save.as="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/LEA/Oxyrhabdium-cf.modestum_Luzon_BestSNP_sNMF_run1.pdf")
#
#	leporinum_Luzon <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium-leporinum_Luzon_BestSNP.vcf",
#	                           coords="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/tess3r/Oxyrhabdium-leporinum_Luzon_coords.txt",
#	                           save.as="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/LEA/Oxyrhabdium-leporinum_Luzon_BestSNP_sNMF_run1.pdf")
#
#	bothmodestum <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium_both-modestum_BestSNP.vcf",
#	                           coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/tess3r/Oxyrhabdium_bothmodestum_coords.txt",
#	                           save.as="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/LEA/Oxyrhabdium_both-modestum_BestSNP_sNMF_run1.pdf")
#
#	Oxyrhabdium <- run_SNMF(vcf="/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/Oxyrhabdium_AllSpecies_BestSNP.vcf",
#	                            coords="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/tess3r/Oxyrhabdium_AllSpecies_coords.txt",
#	                            save.as="/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/LEA/Oxyrhabdium_AllSpecies_BestSNP_sNMF_v1.pdf")
#
#	setwd("/panfs/pfs.local/home/j926w878/scratch/scratch_v2/LEA/")
#	Oxyrhabdium <- run_SNMF(vcf="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Oxyrhabdium_AllSpecies_BestSNP.vcf",
#	                         coords="/panfs/pfs.local/home/j926w878/work/ddRAD/DAPC/Oxyrhabdium_AllSpecies_coords.txt",
#	                         save.as="/panfs/pfs.local/home/j926w878/scratch/scratch_v2/LEA/Oxyrhabdium_AllSpecies_BestSNP_sNMF_run1.pdf")
#
#
##' #library(sNMF)
#' library(tess3r)
#' library(ggplot2)
#' library(rworldmap)
#' library(LEA)
#' library(geosphere)
#' source("/Users/alyssaleinweber/Documents/Chapter4_Bicol-vs-Luzon/DAPC/DAPC_adegenet.R")
#' source("/Users/alyssaleinweber/Documents/Chapter3_Oxyrhabdium/runtess.R")
#}
#




