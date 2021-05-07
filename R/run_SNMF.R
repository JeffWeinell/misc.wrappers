#' @title Wrapper for LEA snmf
#' 
#' Runs LEA function snmf function from VCF input file and optionally interpolates cluster assignments onto a map (if 'coords' argument is non-NULL).
#' 
#' @param vcf Path to input VCF file with SNP data
#' @param coords Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.
#' @param kmax Number indicating the maximum number of clusters to evaluate. Default is 40, which is converted using kmax = min(kmax, number of individuals-1)
#' @param reps Number of repititions. Default 100.
#' @param entropy Default TRUE
#' @param project Defualt 'new'
#' @param iter Default 500
#' @param CPU Defualt 2
#' @param save.as Character string with where to save the output PDF with plots of results. Default is NULL.
#' @return List of plots
#' @export run_sNMF
run_sNMF <- function(vcf,coords=NULL,kmax=40,reps=100,entropy=TRUE,project="new",iter=500,save.as=NULL){
	if(is.null(save.as)){
		save.as <- file.path(getwd(),"result_LEA-sNMF.pdf")
	}
	if(file.exists(save.as)){
		stop("Output file already exists. Use a different name for 'save.as' argument.")
	}
	Krange=1:kmax
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
		world_sf      <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
		world_sp      <- rnaturalearth::ne_countries(scale=10,returnclass="sp")
		current_sf    <- sf::st_crop(world_sf,xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max)
		current.gg.sf <- ggplot2::geom_sf(data=current_sf,colour = "black", fill = NA)
		
	} else {
		maxK <- (numind-1)
	}
	if(max(Krange) > maxK){
		Krange <- 1:maxK
	}
	geno.temp.path   <- paste0(tempfile(),".geno")
	geno.obj         <- vcfR2geno(vcf=vcf,out=geno.temp.path)
	num.CPUs         <- parallel::detectCores()
	if(!is.numeric(num.CPUs)){
		num.CPUs <- 2
	}
	snmf.obj         <- LEA::snmf(geno.temp.path,K=Krange, repetitions=reps,entropy=entropy,project="new",iterations=iter,CPU=num.CPUs)
	crossentropy.mat <- t(do.call(cbind,lapply(X=Krange,FUN=function(x){LEA::cross.entropy(snmf.obj,K = x)})))
	rownames(crossentropy.mat) <- Krange
	colnames(crossentropy.mat) <- paste0("rep",1:reps)
	mean.entropy   <- apply(crossentropy.mat,MARGIN=1,FUN=mean,na.rm=TRUE)
	if(any(diff(mean.entropy)>0)){
		bestK <- unname(which(diff(mean.entropy)>0)[1])
	} else {
		bestK <- unname(Krange[1])
	}
	Kbest.criteria1 <- bestK
	crossentropy.df <- data.frame(crossentropy=unname(unlist(c(crossentropy.mat))),Kval=rep(Krange,reps))
#	mode(crossentropy.df$Kval) <- "character"
	crossentropy.df$Kval <- factor(crossentropy.df$Kval, levels=c(1:nrow(crossentropy.df)))
	entropyPlot <- ggplot2::ggplot(crossentropy.df, ggplot2::aes(x=Kval, y=crossentropy)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("Cross-entropy (",reps," replicates) vs. number of ancestral populations (K)"), x="Number of ancestral populations", y = "Cross-entropy") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
#	boxplot(t(crossentropy.mat))
#	plot(Krange,mean.entropy,pch=21,col="blue",xlab="",ylab="",xlim=range(Krange), ylim=range(range.entropy.mat))
#	arrows(x0=Krange,y0=range.entropy.mat[,1],y1=range.entropy.mat[,2],length=0.07,col="black",angle=90,code=3)
#	lines(Krange,mean.entropy,col="blue")
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
#	if(bestK>1){
#		segments(x0=bestK,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="black")
#	}
	# segments(x0=bestK,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="green")
	# segments(x0=Kbest.criteria2,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="blue")
	# segments(x0=Kbest.criteria3,y0=par("usr")[3],y1=par("usr")[4],lty=2,col="red")
#	mtext(side=1,"Number of ancestral populations",line=2.2)
#	mtext(side=2,"Cross-entropy",line=2.2)
#	mtext(side=3,paste0("Cross-entropy (",reps," replicates) vs. number of ancestral populations (K)"),line=1)
#	axis(1,at=Krange)
	
#	entropyPlot    <- recordPlot()
	Krange.plot    <- setdiff(Krange,1)
	admixturePlot  <- list(); length(admixturePlot)   <- length(Krange.plot)
	assignmentPlot <- list(); length(assignmentPlot)  <- length(Krange.plot)
	mapplot        <- list(); length(mapplot)         <- length(Krange.plot)
	# par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
	for(K in Krange.plot){
		i=(K-1)
		ce           <- LEA::cross.entropy(snmf.obj, K = K)
		best         <- which.min(ce)
	#	q.matrix.best <- suppressWarnings(tess3r::as.qmatrix(LEA::Q(snmf.obj,K=bestK,run=best)))
		if(FALSE){
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
		}
		if(K <= 15){
			myCols          <- goodcolors2(n=K)
		}
		if(K>15){
			myCols          <- c(goodcolors2(n=K), sample(adegenet::funky(100), size=K-15))
		}
		

		q.matrix           <- LEA::Q(snmf.obj,K=K,run=best)
		rownames(q.matrix) <- samplenames
		colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
		posterior.df       <- data.frame(indv=rep(rownames(q.matrix),ncol(q.matrix)), pop=rep(colnames(q.matrix),each=nrow(q.matrix)), assignment=c(unlist(unname(q.matrix))))
		posterior.gg       <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Admixture Proportion",y="",fill="Cluster",title=paste0("K = ",K)) + ggplot2::scale_fill_manual(values=myCols[1:K])
#		plot(posterior.gg)
#		admixturePlot[[i]]   <- recordPlot()
		admixturePlot[[i]]   <- posterior.gg
		

		indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels             <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignment.K       <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
#		assignment.K        <- adegenet::assignplot(dapc.pcabest.K,cex.lab=(label.size/10))
#		mtext(text=paste0("K = ",K,"; PCs retained = ",best.npca[i]))
		assignmentPlot[[i]]  <- assignment.K
	#	

		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols, 9)
			xdist           <- geosphere::distm(x=c(x.min,0),y=c(x.max,0))
			ydist           <- geosphere::distm(x=c(0,y.min),y=c(0,y.max))
			#mapplot.initial <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(100,100), cex = 0.4,col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),asp=xdist/ydist)
			#mapplot.i       <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(500,500), cex = 0.4,col.palette = my.palette, window=c(par("usr")[1],par("usr")[2],par("usr")[3],par("usr")[4]),asp=xdist/ydist)
#			mapplot.initial <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), main = "", xlab = "", ylab = "",resolution = c(2,2), col.palette = lapply(X=1:K,FUN=function(x){rep("#FFFFFF",9)}), cex=0,window=c(x.min,x.max,y.min,y.max),asp=xdist/ydist,add=FALSE)
#			mapplot.i       <- plot(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(500,500), cex = 0.4, col.palette = my.palette, window=par("usr"),asp=xdist/ydist,add=FALSE)
#			maps::map(add=TRUE)
			# mapplot.i       <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=par("usr"),background=TRUE)
			#mapplot.i       <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE)
			mapplot.i       <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE,map.polygon=world_sp)
			# mapplot.i2      <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + borders(database="Philippines",xlim=c(x.min,x.max),ylim=c(y.min,y.max),colour="black") + geom_point(data = coords, aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			# new.x <- ggplot_build(mapplot.i2)$layout$panel_scales_x[[1]]$range$range
			# new.y <- ggplot_build(mapplot.i2)$layout$panel_scales_y[[1]]$range$range
			# mapplot.i3      <- tess3r::ggtess3Q(suppressWarnings(tess3r::as.qmatrix(q.matrix)), as.matrix(coords), interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(new.x[1],new.x[2],new.y[1],new.y[2]),background=TRUE)
			# mapplot.i4      <- mapplot.i3 + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + borders(database="world", xlim=c(x.min,x.max),ylim=c(y.min,y.max),colour="black")
			# mapplot.i4 + geom_point(data = coords, aes(x = Lon, y = Lat), size = 4, shape = 23, fill = "darkred")
			# mapplot.i2 + geom_point(data = coords, aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			# mapplot[[i]]    <- mapplot.i + borders(xlim=c(x.min,x.max),ylim=c(y.min,y.max),colour="black")
#			mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::borders(database="world", xlim=c(x.min,x.max), ylim=c(y.min,y.max), colour="black") + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			#mapplot[[i]]    <- recordPlot()
			# ggplot2::ggplot(data = world) + ggplot2::geom_sf() + ggplot2::coord_sf(xlim = c(x.min, x.max), ylim = c(y.min, y.max), expand = FALSE)
		}
	}
	if(!is.null(coords)){
		result <- c(list(entropyPlot),admixturePlot,assignmentPlot,mapplot)
	} else {
		result <- c(list(entropyPlot),admixturePlot,assignmentPlot)
	}
	if(!is.null(save.as)){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	# result
	# result.grob <- lapply(X=result,FUN=grob)
	# pl <- lapply(X=list(entropyPlot, mapplot),FUN=grob)
	# ml <- marrangeGrob(pl, nrow=1, ncol=1)
	# ggsave(save.as, ml)


}
#' @examples
#'	library(ade4)
#'	library(adegenet)
#'	library(vcfR)
#'	library(ggplot2)
#'	library(maps)
#'	library(geosphere)
#'	library(tess3r)
#'	library(LEA)
#'	library(rworldmap)
#'	library(JeffWeinell/misc.wrappers)
#'	#source("DAPC_adegenet.R")
#'	#source("runtess.R")
#'	#source("run_sNMF.R")
#'
#'	dev.new(width=10,height=6)
#'	leporinum <- run_sNMF(vcf="Oxyrhabdium-leporinum_BestSNP.vcf",
#'	                      coords="Oxyrhabdium-leporinum_coords.txt",
#'	                      save.as="Oxyrhabdium-leporinum_BestSNP_sNMF_run3.pdf")
#'
#'	modestum  <- run_sNMF(vcf="Oxyrhabdium-modestum_BestSNP.vcf",
#'	                     coords="Oxyrhabdium-modestum_coords.txt",
#'	                     save.as="Oxyrhabdium-modestum_BestSNP_sNMF_run1.pdf")
#'
#'	cfmodestum  <- run_sNMF(vcf="Oxyrhabdium-cf.modestum_BestSNP.vcf",
#'	                     coords="Oxyrhabdium-cfmodestum_coords.txt",
#'	                     save.as="Oxyrhabdium-cfmodestum_BestSNP_sNMF_run1.pdf")
#'
#'	cfmodestum_Luzon     <- run_sNMF(vcf="Oxyrhabdium-cf.modestum_Luzon_BestSNP.vcf",
#'	                        coords="Oxyrhabdium-cfmodestum_Luzon_coords.txt",
#'	                        save.as="Oxyrhabdium-cf.modestum_Luzon_BestSNP_sNMF_run1.pdf")
#'
#'	leporinum_Luzon <- run_sNMF(vcf="Oxyrhabdium-leporinum_Luzon_BestSNP.vcf",
#'	                           coords="Oxyrhabdium-leporinum_Luzon_coords.txt",
#'	                           save.as="Oxyrhabdium-leporinum_Luzon_BestSNP_sNMF_run1.pdf")
#'
#'	bothmodestum <- run_sNMF(vcf="Oxyrhabdium_both-modestum_BestSNP.vcf",
#'	                           coords="Oxyrhabdium_bothmodestum_coords.txt",
#'	                           save.as="Oxyrhabdium_both-modestum_BestSNP_sNMF_run1.pdf")
#'
#'	Oxyrhabdium <- run_sNMF(vcf="Oxyrhabdium_AllSpecies_BestSNP.vcf",
#'	                            coords="Oxyrhabdium_AllSpecies_coords.txt",
#'	                            save.as="Oxyrhabdium_AllSpecies_BestSNP_sNMF_v1.pdf")
#'
#'	Oxyrhabdium_AllSpecies <- run_sNMF(vcf="Oxyrhabdium_AllSpecies_BestSNP.vcf",
#'	                         coords="Oxyrhabdium_AllSpecies_coords.txt",
#'	                         save.as="Oxyrhabdium_AllSpecies_BestSNP_sNMF_run1.pdf")


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
	vcf.obj   <- vcfR::read.vcfR(vcf,verbose=FALSE)
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



