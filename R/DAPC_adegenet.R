#' @title Run DAPC and Plot
#' 
#' Run adagenet DAPC analyses and generate plots of BIC vs number of clusters, baplot of alpha optimum number of principle components at each K, admixture for each K, assignments at each K.
#' If sample coordinates are supplied, this function interpolates cluster membership probabilities on a map, using functions from tess3r package.
#' 
#' @param vcf Character string with path to vcf file containing snp data
#' @param kmax Number indicating the maximum number of clusters to evaluate. Default is NULL, in which case kmax is set to one less than the number of individuals.
#' @param coords Optional character string with path to a table with longitude and latitude of individuals in the vcf file, or a matrix or data frame with longitude and latitude columns. Default is NULL, in which case membership probabilities are not interpolated onto a map.
#' @param reps Number indicating the number of replicates of 'find.clusters'. Default 100.
#' @param probs.out NULL or a character string with location where to write a table containing the membership probabilities for the best K and alpha-optimized number of PCAs.
#' @param save.as Character string with where to save the output PDF with plots of results. Default is NULL. **Important! This argument is ignored in some environments. Instead, use dev.new(file="Where/To/Save/Output.pdf",height=6,width=10,noRStudioGD=TRUE) before using run_DAPC. Then dev.off().
#' @return A list of plots.
#' @export run_DAPC
run_DAPC <- function(vcf, kmax=50, coords=NULL, reps=100,probs.out=NULL,save.as=NULL){
	#dev.new(width=10,height=6)
	vcf.obj     <- vcfR::read.vcfR(vcf,verbose=F)
	samplenames <- colnames(vcf.obj@gt)[-1]
	genind      <- suppressWarnings(vcfR::vcfR2genind(vcf.obj))
	numind      <- (dim(attributes(vcf.obj)[["gt"]])[2])-1
	label.size  <- min((288/numind),7)
	if(!is.null(coords)){
		if(is(coords,"array") | is(coords,"data.frame")){
			coords <-  coords
		} else {
			if(file.exists(coords)){
				coords   <- read.table(coords)
			}
		}
		x.min <- min((coords[,1]-0.5))
		x.max <- max((coords[,1]+0.5))
		y.min <- min((coords[,2]-0.5))
		y.max <- max((coords[,2]+0.5))
		maxK  <- min(nrow(unique(coords)),(numind-1))
		#world_df <- ggplot2::map_data(rnaturalearth::ne_countries(scale=10))
		world_sf      <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
		world_sp      <- rnaturalearth::ne_countries(scale=10,returnclass="sp")
		current_sf    <- sf::st_crop(world_sf,xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max)
		current.gg.sf <- ggplot2::geom_sf(data=current_sf,colour = "black", fill = NA)
		#world.gg.sf   <- ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA)
	} else {
		maxK <- (numind-1)
	}
	if(is.null(kmax)){
		kmax <- maxK
	} else {
		if(kmax > maxK){
			kmax <- maxK
		}
	}
	max.clusters <- kmax
	Krange       <- 1:max.clusters
	grp          <- adegenet::find.clusters(genind, max.n.clust=max.clusters,n.pca=max.clusters,choose.n.clust=F)
	grp.list <- list(); length(grp.list) <- reps
	par(mar=c(3.5,4,3,2.1))
	for(i in 1:reps){
		grp.list[[i]] <- adegenet::find.clusters(genind, max.n.clust=max.clusters,n.pca=max.clusters,choose.n.clust=F)
	}
	BIC.mat           <- do.call(cbind,lapply(X=1:reps,FUN=function(x){grp.list[[x]]$Kstat}))
	rownames(BIC.mat) <- 1:max.clusters
	colnames(BIC.mat) <- paste0("rep",1:reps)
	mean.BIC        <- apply(BIC.mat,MARGIN=1,FUN=mean,na.rm=TRUE)
	### Lowest K with a lower mean BIC than K+1 mean BIC.
	if(any(diff(mean.BIC)>0)){
		bestK <- unname(which(diff(mean.BIC)>0)[1])
	} else {
		bestK <- 1
	}
	Kbest.criteria1   <- bestK
#	boxplot(t(BIC.mat))
	BIC.df      <- data.frame(BIC=unname(unlist(c(BIC.mat))),Kval=rep(Krange,reps))
	#BIC.df <- BIC.df[order(BIC.df.temp[,"Kval"]),]
	#mode(BIC.df$Kval) <- "character"
	BIC.df$Kval <- factor(BIC.df$Kval, levels=c(1:nrow(BIC.df)))
	BICPlot     <- ggplot2::ggplot(data=BIC.df,ggplot2::aes(x=Kval, y=BIC)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("BIC (",reps," replicates of find.clusters) vs. number of clusters (K)"), x="Number of ancestral populations", y = "BIC") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
	### Checking if the max BIC of the next K value is less than the min BIC of the previous K value.
	range.BIC.mat     <- do.call(rbind,lapply(X=1:nrow(BIC.mat),FUN=function(x){range(BIC.mat[x,],na.rm=TRUE)}))
	BIC.is.nonoverlap <- NULL
	BIC.is.reduced    <- NULL
	for(i in 2:nrow(range.BIC.mat)){
		BIC.is.nonoverlap     <- c(BIC.is.nonoverlap,range.BIC.mat[i,1] > range.BIC.mat[(i-1),2] | range.BIC.mat[i,2] < range.BIC.mat[(i-1),1])
		BIC.is.reduced     <- c(BIC.is.reduced,range.BIC.mat[i,2] < range.BIC.mat[(i-1),1])
	}
	# If max BIC for K=2 is better than BIC K=1 and if some K values are not better (all BIC lower) than K-1, then find the first K value in which K+1 is not better.
	if(any(!BIC.is.reduced) & BIC.is.reduced[1]){
		Kbest.criteria2 <- which(!BIC.is.reduced)[1]
	} else {
		Kbest.criteria2 <- 1
	}
	### Which K value (for K>=2) yields the least variable BIC scores.
	BICvK.variation <- apply(X=BIC.mat,MARGIN=1,FUN=var)
	KminVarBIC      <- which(BICvK.variation==min(BICvK.variation[-1]))
	Kbest.criteria3 <- KminVarBIC[1]
	### Criteria 4: t-tests for BIC of each pairwise adjacent K
	for(i in 2:nrow(BIC.mat)){
		if(BICvK.variation[Kbest.criteria3]==0){
			Kbest.criteria4 <- NULL
			break
		}
		t.test.i <- t.test(BIC.mat[i-1,],BIC.mat[i,])
		pval.i   <- t.test.i$p.value
		stat.i   <- t.test.i$statistic
		if(pval.i < 0.05 & stat.i > 0){
			next
		} else {
			if(i==nrow(BIC.mat)){
				Kbest.criteria4 <- NULL
			} else{
				Kbest.criteria4 <- (i-1)
				break
			}
		}
	}
#	if(bestK>1){
#		segments(x0=Kbest.criteria1, y0=par("usr")[3], y1=par("usr")[4],lty=2,col="black")
#	}
	# segments(x0=Kbest.criteria1, y0=par("usr")[3], y1=par("usr")[4],lty=2,col="green")
	# segments(x0=Kbest.criteria2, y0=par("usr")[3], y1=par("usr")[4],lty=2,col="blue")
	# segments(x0=Kbest.criteria3, y0=par("usr")[3], y1=par("usr")[4],lty=2,col="orange")
#	mtext(side=1,"Number of ancestral populations",line=2.2)
#	mtext(side=2,"BIC",line=2.2)
#	mtext(side=3,paste0("BIC (",reps," replicates of find.clusters) vs. number of clusters (K)"),line=1)
#	axis(1,at=1:max.clusters)
#	BICPlot    <- recordPlot()
	best.npca  <- NULL
	grp.mat    <- matrix(data=0,nrow=length(grp$grp),ncol=(max.clusters-1))
	for(K in 2:max.clusters){
		i=(K-1)
		grp.K         <- adegenet::find.clusters(genind, max.n.clust=max.clusters,n.pca=max.clusters,n.clust=K)
		grp.mat[,i]   <- grp.K$grp
		dapc.pcamax.K <- suppressWarnings(adegenet::dapc(genind, grp.K$grp,n.pca=max.clusters,n.da=5))
		alpha_optim.K <- suppressWarnings(adegenet::optim.a.score(dapc.pcamax.K,plot=FALSE))
		best.npca     <- c(best.npca,alpha_optim.K$best)
	}
	##### Plot 2: BIC vs. K when number of PCs retained = alpha optimized 
	par(mar=c(4,4,3,2.1))
	# names(best.npca) <- 2:max.clusters
	best.npca.df <- data.frame(best.npca=best.npca,Kval=2:max.clusters)
	best.npca.df$Kval <- factor(best.npca.df$Kval)
	grp.plot2    <- ggplot2::ggplot(data=best.npca.df, ggplot2::aes(x=Kval,y=best.npca)) + ggplot2::geom_bar(stat="identity",fill="lightgray") + ggplot2::labs(title= "alpha optimized # of PCs vs. number of clusters", x="Number of clusters", y = "Alpha optimized number of principle components to retain") + theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
#	barplot(best.npca)
#	mtext(text="Alpha optimized number of principle components to retain",side=2,line=2)
#	mtext(text="Number of clusters",side=1,line=2)
#	mtext(text="alpha optimized # of PCs vs. number of clusters",side=3,line=1)
#	grp.plot2      <- recordPlot()
	admixturePlot  <- list(); length(admixturePlot)   <- max.clusters-1
	assignmentPlot <- list(); length(assignmentPlot)  <- max.clusters-1
	posterior.list <- list(); length(posterior.list)  <- max.clusters-1
	mapplot        <- list(); length(mapplot)  <- max.clusters-1
	
	for(K in 2:max.clusters){
		i=(K-1)
		par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
		dapc.pcabest.K  <- adegenet::dapc(genind, grp.mat[,i],n.pca=best.npca[i],n.da=5)
		posterior       <- dapc.pcabest.K$posterior
		q.matrix        <- posterior
		posterior.list[[i]] <- posterior
		posterior.df    <- data.frame(indv=rep(rownames(posterior),ncol(posterior)), pop=rep(colnames(posterior),each=nrow(posterior)), assignment=c(posterior))
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
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Membership Probability",y="",fill="Cluster",title=paste0("K = ",K,"; PCs retained = ",best.npca[i])) + ggplot2::scale_fill_manual(values=myCols[1:K])
		admixturePlot[[i]]  <- posterior.gg
		par(mar=c(5,20,2,2.1))
		# test
		#assignment.K  <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv, fill= assignment)) + ggplot2::geom_tile() + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(title = paste0("K = ",K,"; PCs retained = ",best.npca[i]), x="cluster", y="") + ggplot2::scale_colour_gradientn(colours = c("yellow", "orange", "red")) # + scale_fill_brewer(palette = "YlOrRd",trans="probability")
		#assignment.K   <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv), fill= assignment) + ggplot2::geom_tile() + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(title = paste0("K = ",K,"; PCs retained = ",best.npca[i]), x="cluster", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5)  # ggplot2::scale_colour_gradient2(colours = c("yellow", "red")) # + scale_fill_brewer(palette = "YlOrRd",trans="probability")
		indv.KmaxPosterior <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))})
		indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels             <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignment.K       <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K,"; PCs retained = ", best.npca[i]), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)

#		assignment.K        <- adegenet::assignplot(dapc.pcabest.K,cex.lab=(label.size/10))
#		mtext(text=paste0("K = ",K,"; PCs retained = ",best.npca[i]))
		assignmentPlot[[i]]  <- assignment.K
		par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols, 9)
		#	xdist           <- geosphere::distm(x=c(x.min,0),y=c(x.max,0))
		#	ydist           <- geosphere::distm(x=c(0,y.min),y=c(0,y.max))
		#	xdist2          <- (ydist*(10/6))
		#	xbuff.deg       <- ((xdist2-xdist)/2)/(111111*cos(((y.max-y.min)*pi)/180))
		#	x.min2          <- x.min-xbuff.deg
		#	x.max2          <- x.max+xbuff.deg
			tess3r.qmat     <- suppressWarnings(tess3r::as.qmatrix(q.matrix))
			coords.mat      <- as.matrix(coords)
		#	borders.aes     <- borders(xlim = c(x.min,x.max), ylim=c(y.min,y.max),colour="black")
		#	gg.extent       <- ggplot() + borders(xlim = c(x.min,x.max), ylim=c(y.min,y.max),colour="black")
		#	new.x           <- ggplot_build(new.extent)$layout$panel_scales_x[[1]]$range$range
		#	new.y           <- ggplot_build(new.extent)$layout$panel_scales_y[[1]]$range$range
			#mapplot.initial <- plot(tess3r.qmat, coords.mat, method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(100,100), cex = 0.4,col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),asp=xdist/ydist)
			#mapplot.i       <- plot(tess3r.qmat, coords.mat, method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(500,500), cex = 0.4,col.palette = my.palette, window=c(par("usr")[1],par("usr")[2],par("usr")[3],par("usr")[4]),asp=xdist/ydist)
		#	mapplot.initial <- plot(tess3r.qmat, coords.mat, main = "", xlab = "", ylab = "",resolution = c(2,2), col.palette = lapply(X=1:K,FUN=function(x){rep("#FFFFFF",9)}), cex=0,window=c(x.min,x.max,y.min,y.max),asp=xdist/ydist,add=FALSE)
		#	mapplot.i       <- plot(tess3r.qmat, coords.mat, method = "map.max", interpol = tess3r::FieldsKrigModel(10), main = paste0("Ancestry coefficients; K=",K), xlab = "", ylab = "",resolution = c(500,500), cex = 0.4, col.palette = my.palette, window=par("usr"),asp=xdist/ydist,add=FALSE)
		#	maps::map(add=TRUE)
		#	mapplot[[i]]    <- recordPlot()
		#	mapplot.i       <- tess3r::ggtess3Q(tess3r.qmat,coords.mat, interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(new.x,new.y),background=TRUE,map.polygon=world_sp)
			mapplot.i       <- tess3r::ggtess3Q(tess3r.qmat,coords.mat, interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE,map.polygon=world_sp)
		#	mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + borders.aes + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			# test <- map_data("world")
			#world_sf   <- sf::st_as_sf(rnaturalearth::ne_countries(scale=10))[1]
			#world_sf      <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
			#gg.current <- mapplot.i + current.gg.sf
			#gg.world   <- world.gg.sf
		}
	}

	if(bestK>1){
		posterior.bestK <- posterior.list[[bestK-1]]
		colnames(posterior.bestK) <- paste0("K",colnames(posterior.bestK))
		if(!is.null(probs.out)){
			write.table(posterior.bestK,file=probs.out,quote=F,col.names=T,row.names=T)
		}
	} else {
		posterior.bestK <- matrix(data=rep(1,numind),ncol=1)
		rownames(posterior.bestK) <- samplenames
	}
	#dev.off()
	if(!is.null(coords)){
		result <- c(list(BICPlot,grp.plot2),admixturePlot,assignmentPlot,mapplot)
	} else {
		result <- c(list(BICPlot,grp.plot2),admixturePlot,assignmentPlot)
	}
	if(!is.null(save.as)){
		pdf(height=6,width=10,file=save.as,onefile=TRUE)
		lapply(X=result,FUN=print)
		dev.off()
	}
	result
}
#######
# module load R/4.0
# module load gdal
# source("rundelimitR_functions.R")
# source("DAPC_adegenet.R")
# library(ade4)
# library(adegenet)
# library(vcfR)
# library(ggplot2)
#' 
#' dev.new(height=10,width=6)
#' Ahaetulla_Luzon <- run_DAPC(vcf="Ahaetulla-prasina_Luzon_1snpPerLocus.vcf",
#'                             probs.out="Ahaetulla-prasina_Luzon_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                             save.out="Ahaetulla-prasina_Luzon_DAPC_v2.pdf")

#' Calamaria_allpops <- run_DAPC(vcf="Calamaria-gervaisii_snps.vcf",
#'                             probs.out="Calamaria-gervaisii_allpops_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                             save.as="Calamaria-gervaisii_allpops_DAPC.pdf")

#' dev.new(height=6,width=10)
#' cfmodestum_Luzon <- run_DAPC(vcf="Oxyrhabdium-cf.modestum_Luzon_BestSNP.vcf",
#'                             probs.out="Oxyrhabdium-cf.modestum_Luzon_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                             coords="tess3r/Oxyrhabdium-cfmodestum_Luzon_coords.txt",
#'                             save.as="Oxyrhabdium-cfmodestum_Luzon_BestSNP_DAPC_run2.pdf")

#' leporinum_Luzon <- run_DAPC(vcf="Oxyrhabdium-leporinum_Luzon_BestSNP.vcf".
#'                             probs.out="Oxyrhabdium-leporinum_Luzon_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                             coords="Oxyrhabdium-leporinum_Luzon_coords.txt",
#'                             save.as="Oxyrhabdium-leporinum_Luzon_BestSNP_DAPC.pdf")

#' leporinum  <- run_DAPC(vcf="Oxyrhabdium-leporinum_BestSNP.vcf",
#'                             probs.out="Oxyrhabdium-leporinum_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                             coords="Oxyrhabdium-leporinum_coords.txt",
#'                             save.as="Oxyrhabdium-leporinum_BestSNP_DAPC_run2.pdf")

#' modestum   <- run_DAPC(vcf="Oxyrhabdium-modestum_BestSNP.vcf",
#'                             probs.out="Oxyrhabdium-modestum_BestSNP_DAPC_BestK_Membership.txt_run2",
#'                             coords="Oxyrhabdium-modestum_coords.txt",
#'                             save.as="Oxyrhabdium-modestum_BestSNP_DAPC_run2.pdf")

#' cfmodestum <- run_DAPC(vcf="Oxyrhabdium-cf.modestum_BestSNP.vcf",
#'                             probs.out="Oxyrhabdium-cf.modestum_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                             coords="Oxyrhabdium-cfmodestum_coords.txt",
#'                             save.as="Oxyrhabdium-cf.modestum_BestSNP_DAPC_run2.pdf")

#'Oxyrhabdium <- run_DAPC(vcf="Oxyrhabdium_AllSpecies_BestSNP.vcf",
#'                            probs.out="Oxyrhabdium_AllSpecies_BestSNP_DAPC_BestK_Membership_run2.txt",
#'                            coords="Oxyrhabdium_AllSpecies_coords.txt",
#'                            save.as="Oxyrhabdium_AllSpecies_BestSNP_DAPC_v3.pdf")

#' @title Get Best SNP for each locus from VCF
#' 
#' From a VCF with multiple sites/locus, create a VCF with only the best site/locus.
#' The best site is the one with the least missing data. To break ties, take the first site among the best sites.
#' 
#' @param vcf Character string with path to input vcf.
#' @param vcftools.path Character string with path to the vcftools executable.
#' @param out Character string where to write output vcf.
#' @param indv.keep Character string with names of individuals to keep. Default is NULL (all individuals kept).
#' @param min.n Minimum number of individuals required to keep a site. Default = 4.
#' @param min.n0 Minimum number of individuals required to have the major allele to keep a site. Default = 2.
#' @param min.n1 Minimum number of individuals required to have the minor allele to keep a site. Default = 1.
#' @param which.site Character string indicating the method for choosing a site to keep for each locus (or chromosome). Default = "best", which is considered the one with the least missing data, or the first among sites tied for least missing data. Other options are "all.passing", which retains all sites (positions) that pass variation filters (min.n, min.0.n.0, min.1.n), "first" (first site kept at each locus), or "random".
#' @return Character vector with values of input arguments
#' @export vcf_getSNP
vcf_getSNP     <- function(vcftools.path,vcf,out,indv.keep=NULL,min.n=4,min.n0=2,min.n1=1,which.site="best"){
	vcf.obj    <- vcfR::read.vcfR(vcf)
	# matrix with "fixed" columns, which are the columns with site-specific stats across all samples
	# fix.mat    <- attributes(vcf.obj)[["fix"]]
	fix.mat    <- vcf.obj@fix
	# matrix with genotype columns, which are the columns with site-specific stats across all samples
	# gt.mat     <- attributes(vcf.obj)[["gt"]]
	# gt.mat     <- vcf.obj@gt[,-1]
	# Remove first column of gt
	#gt.mat     <- gt.mat[,-1]
	# Remove everything after ":" in strings
	#gt.mat     <- gsub(":.+","",gt.mat)
	gt.mat     <- gsub(":.+","",vcf.obj@gt[,-1])
	if(!is.null(indv.keep)){
		if(is(indv.keep,"character") & length(indv.keep)==1){
			if(file.exists(indv.keep)){
				#indv.keep <- unname(unlist(read.table(indv.keep)))
				indv.keep <- suppressWarnings(readLines(indv.keep))
			}
		}
		### Check that all names in indv.keep actually exist in the vcf
		if(all(indv.keep %in% colnames(gt.mat))){
			gt.mat <- gt.mat[,indv.keep,drop=F]
		} else {
			stop(paste(paste(setdiff(indv.keep,colnames(gt.mat)),collapse=","), "not in VCF"))
		}
	}
	# For each site, the number of non-missing alleles
	site.NS      <- vapply(X=1:nrow(gt.mat),FUN=function(x){length(grep(".", unlist(strsplit(gt.mat[x,],split="/",fixed=T)),fixed=T,invert=T))},FUN.VALUE=1)
	# For each site, the number of individuals with at least one copy of the major allele
	site.0.NS    <- vapply(X=1:nrow(gt.mat),FUN=function(x){length(grep("0", gt.mat[x,],fixed=T,invert=F))},FUN.VALUE=1)
	# For each site, the number of individuals with at least one copy of the minor allele
	site.1.NS    <- vapply(X=1:nrow(gt.mat),FUN=function(x){length(grep("1", gt.mat[x,],fixed=T,invert=F))},FUN.VALUE=1)
	# Matrix with "CHROM" and "POS" from fix.mat, plus columns containing site.NS, site.0.NS, and site.1.NS
	chrom.pos.mat <- cbind(fix.mat[,c("CHROM","POS")],site.NS,site.0.NS,site.1.NS)
	
	chrom.pos.df <- data.frame(CHROM=fix.mat[,"CHROM"],POS=fix.mat[,"POS"],site.NS=site.NS,site.0.NS=site.0.NS,site.1.NS=site.1.NS)
	mode(chrom.pos.df[,"CHROM"])     <- "character"
	mode(chrom.pos.df[,"POS"])       <- "character"
	mode(chrom.pos.df[,"site.NS"])   <- "numeric"
	mode(chrom.pos.df[,"site.0.NS"]) <- "numeric"
	mode(chrom.pos.df[,"site.1.NS"]) <- "numeric"
	if(any(chrom.pos.df[,"site.NS"] >= min.n & chrom.pos.df[,"site.0.NS"] >= min.n0 & chrom.pos.df[,"site.1.NS"] >= min.n1)){
		# filter1 <- which(chrom.pos.df[,"site.1.NS"] >=2)
		filter1 <- which(chrom.pos.df[,"site.NS"] >= min.n & chrom.pos.df[,"site.0.NS"] >= min.n0 & chrom.pos.df[,"site.1.NS"] >= min.n1)
	} else {
		stop("no sites in which the minor allele occurs in more than 1 individual")
	}
	chrom.pos.df.filtered1 <- chrom.pos.df[filter1,]
	# vector with loci names for each site retained
	loci        <- chrom.pos.df.filtered1[,"CHROM"]
	# vector with unique loci names
	loci.unique <- unique(loci)
	# For each unique locus ("CHROM") in chrom.pos.df.filtered1, return the row number for the site with the greatest site.NS
	chrom.pos.list.filtered2 <- list(); length(chrom.pos.list.filtered2) <- length(loci.unique)
	for(i in 1:length(loci.unique)){
		locus.i <- loci.unique[i]
		chrom.pos.df.i <- chrom.pos.df.filtered1[which(chrom.pos.df.filtered1[,"CHROM"]==locus.i),,drop=F]
		site.NS.i <- chrom.pos.df.i[,"site.NS"]
		max.i     <- max(site.NS.i)
		if(which.site=="best"){
			chrom.pos.list.filtered2[[i]]    <- chrom.pos.df.i[(which(site.NS.i==max.i)[1]),]
		}
		if(which.site=="first"){
			chrom.pos.list.filtered2[[i]]    <- chrom.pos.df.i[1,]
		}
		if(which.site=="random"){
			chrom.pos.list.filtered2[[i]]    <- chrom.pos.df.i[sample(1:length(site.NS.i),1),]
		}
		if(which.site=="all.passing"){
			chrom.pos.list.filtered2[[i]]    <- chrom.pos.df.i
		}
	}
	chrom.pos.df.filtered2 <- do.call(rbind, chrom.pos.list.filtered2)
	# Remove rows of chrom.pos.filtered.mat if the best site for the locus has fewer than four individuals with data
	chrom.pos.filtered.df <- chrom.pos.df.filtered2[,c("CHROM","POS")]
	# Write chrom.pos.filtered.mat to a file as a tab-separated table. This file will be used by vcftools to extract only these sites from the input vcf.
	temp.file <- file.path(getwd(),"temp_filtertable.txt")
	write.table(chrom.pos.filtered.df, temp.file, quote=F, sep="\t", col.names=F, row.names=F)
	# Write list of individuals that should be included in the final set
	temp.file2 <- file.path(getwd(),"temp_filtertable2.txt")
	writeLines(text=colnames(gt.mat),con=temp.file2)
	vcf.command <- paste(vcftools.path,"--vcf",vcf,"--positions",temp.file,"--keep",temp.file2,"--recode -c > ",out)
	system(vcf.command)
	system(paste("rm",temp.file))
	system(paste("rm",temp.file2))
	log.path <- file.path(getwd(),"out.log")
	system(paste("rm",log.path))
	return(c(vcftools.path=vcftools.path, vcf=vcf, out=out, nloci.out=length(unique(chrom.pos.filtered.df[,"CHROM"])), nsites.out=nrow(chrom.pos.filtered.df),numind.out=ncol(gt.mat)))
}
#' @examples
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium_AllSpecies_AllSNPs.vcf",out="Oxyrhabdium_AllSpecies_BestSNP.vcf")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium-modestum_AllSNPs.vcf",out="Oxyrhabdium-modestum_BestSNP.vcf")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium-cf.modestum_AllSNPs.vcf",out="Oxyrhabdium-cf.modestum_BestSNP.vcf")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium-cf.modestum_Luzon_AllSNPs.vcf",out="Oxyrhabdium-cf.modestum_Luzon_BestSNP.vcf")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium-leporinum_AllSNPs.vcf",out="Oxyrhabdium-leporinum_BestSNP.vcf")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium-leporinum_Luzon_AllSNPs.vcf",out="Oxyrhabdium-leporinum_Luzon_BestSNP.vcf")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium_both-modestum_AllSNPs.vcf",out="Oxyrhabdium_both-modestum_AllSNPs.vcf",which.site="all.passing")
#'
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium_both-modestum_AllSNPs.vcf",out="Oxyrhabdium_both-modestum_BestSNP.vcf",which.site="best")
#' 
#' vcf_getSNP(vcftools.path="/vcftools",vcf="Oxyrhabdium_AllSpecies_AllSNPs.vcf", indv.keep="indv_keep_Oxyrhabdium_both-modestum.txt", out="Oxyrhabdium_both-modestum_BestSNP.vcf", which.site="best")
#####

#' @title Hex to xyz colors
#' 
#' Hex to xyz color coordinates
#' 
#' @param hex Character string with hex color code
#' @return xyz color coordinates
#' @export hex2dec
hex2dec <- function(hex){
	vhex  <- unlist(strsplit(hex,split=""))
	if(length(vhex)==7){
		vhex  <- vhex[2:7]
	}
	part1 <- paste0("0x",paste0(vhex[1:2],collapse=""))
	part2 <- paste0("0x",paste0(vhex[3:4],collapse=""))
	part3 <- paste0("0x",paste0(vhex[5:6],collapse=""))
	return(strtoi(c(part1,part2,part3)))
}

#' @title Cartesian color distance
#'
#' Cartesian distance between two colors supplied in hexidecimal
#' 
#' @param col1 Hex color code for color #1
#' @param col2 Hex color code for color #2
#' @return Number with cartesional color distance
#' @export diffcol
diffcol <- function(col1,col2){
	if(is.character(col1) & is.character(col2)){
		xyz1 <- hex2dec(col1)
		xyz2 <- hex2dec(col2)
	} else {
		if(is.numeric(col1) & length(col1) == 3 & is.numeric(col2) & length(col2) == 3){
			xyz1 <- col1
			xyz2 <- col2
		} else {
			stop("col1 and col2 must be hex or rgb")
		}
	}
	diff <- sqrt(((xyz2[1]-xyz1[1])^2) + ((xyz2[2]-xyz1[2])^2) + ((xyz2[3]-xyz1[3])^2))
	return(diff)
}

#' @title Powerset
#' 
#' Function to return the powerset for a numeric or character vector of set elements.
#' 
#' @param x Numerical or character vector
#' @return List with powerset of x
#' @export pset
pset <- function (x) {
	m   <- length(x)
	out <- list(x[c()])
	for (i in seq_along(x)) {
		out <- c(out, lapply(X=out[lengths(out) < m], FUN=function(y){c(y,x[i])}))
	}
	out
}

#' @title Transform sRGB colorspace into colorblind colorspaces
#'
#' Uses the model described here: https://ixora.io/projects/colorblindness/color-blindness-simulation-research/
#'
#' @param srgb Numerical vector with R,G,B color coordinates
#' @return Numeric matrix with transformed RGB color coordinates. First row contains input coordinates. Second through fourth rows contain coordinates for Protanopia, Deuteranopia, and Tritanopia colors.
#' @export rgb.cb
rgb.cb <- function(srgb){
	lrgb.list <- list(); length(lrgb.list) <- 3
	if(any(srgb <= (0.04045*255))){
		lrgb.list[which(srgb <= (0.04045*255))] <- ((srgb[which(srgb <= (0.04045*255))]/255)/12.92)
	}
	if(any(srgb > (0.04045*255))){
		lrgb.list[which(srgb > (0.04045*255))] <- ((((srgb[which(srgb > (0.04045*255))]/255)+0.055)/1.055)^2.4)
	}
	lrgb    <- unlist(lrgb.list)
	Tlms    <- matrix(data=c(0.31399022, 0.15537241, 0.01775239, 0.63951294, 0.75789446, 0.10944209, 0.04649755, 0.08670142, 0.87256922),ncol=3,nrow=3,byrow=F)
	## inverse of Tlms
	Tinv    <- solve(Tlms)
	### The input colors in the long, medium, short wavelength (lms) colorspace
	lms  <- Tlms%*%lrgb
	S.prot <- matrix(data=c(0,0,0,1.05118294,1,0,-0.05116099,0,1),ncol=3,nrow=3,byrow=F)
	S.deut <- matrix(data=c(1,0.9513092,0,0,0,0,0,0.04866992,1),ncol=3,nrow=3,byrow=F)
	S.trit <- matrix(data=c(1,0,-0.86744736,0,1,1.86727089,0,0,0),ncol=3,nrow=3,byrow=F)
	### LMS transformed color space for Protanopia, Deuteranopia, and Tritanopia.
	lms.prot  <- S.prot%*%lms
	lms.deut  <- S.deut%*%lms
	lms.trit  <- S.trit%*%lms
	### Going back to lRGB world
	lrgb.prot <- Tinv%*%lms.prot
	lrgb.deut <- Tinv%*%lms.deut
	lrgb.trit <- Tinv%*%lms.trit
	### Back to RGB world
	inv.v <- function(vp){
		res <- list(); length(res) <- 3
		if(any(vp <= 0.0031308 )){
			res[which(vp <= 0.0031308)] <- ((vp[which(vp <= 0.0031308)]*255)*12.92)
		}
		if(any(vp > 0.0031308)){
			res[which(vp > 0.0031308)] <- ((((vp[which(vp > 0.0031308)]*1.055)^0.41666)-0.055)*255)
		}
		res <- unlist(res)
		res
	}
	rgb.prot <- inv.v(lrgb.prot)
	rgb.deut <- inv.v(lrgb.deut)
	rgb.trit <- inv.v(lrgb.trit)
	### 
	result <- rbind(srgb,rgb.prot,rgb.deut,rgb.trit)
	colnames(result) <- c("r","g","b")
	rownames(result) <- c("norm","prot","deut","trit")
	result
}

#' @title Return a good set of colors to use for clusters in admixture plots.
#' 
#' Generate an optimal set of colors for an admixture plot and particular K (number of clusters), while considering colorblindness.
#' 
#' @param n Number of colors to include in the output set.
#' @param thresh Minimum color distance between any pair of colors in the output set. This will be ignored if not met after 'iter' attempts. Default is 65. Values lower than this tend to result in some colors being too similar to easily distinguish. It will be difficult to find sets with more than 10 colors sets larger than 10 colors (with cbspace=NULL) if thresh is >80.
#' @param iter Number of times to try to find a set of colors that meet 'thresh'. 
#' @param cbspace Character vector indicating which colorblind-transformed sets of colors should also meet thresh. Default is "prot","deut","trit", but it becomes difficult to find a passing color set when n > 4 using the defaults.
#' @return Character vector of hex color codes.
#' @export goodcolors
goodcolors <- function(n,thresh=65,iter=50,cbspace=c("prot","deut","trit")){
	# "#33E8BB" "#C69A01" "#799525" "#030A34" "#091E57" "#4ED671" "#A4282F" "#2F08BA" "#5BD80A" "#CD0B50" "#B9C1F5"
	# funky100          <- adegenet::funky(100)
	pset.index  <- pset(1:n)[which(lengths(pset(1:n))==2)]
#	cols.list   <- lapply(X=1:n,FUN=function(x){sample(0:255,3)})
#	colors.temp <- do.call(rbind,cols.list)
#	pairs.colorsdiffs <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(cols.list[[A[1]]], cols.list[[A[2]]])}))
	norm.cols <- lapply(X=1:n,FUN=function(x){sample(0:255,3)})
	colors.temp <- do.call(rbind,norm.cols)
	all.cols  <- lapply(X=norm.cols,FUN=rgb.cb)
	prot.cols <- lapply(1:length(all.cols),FUN=function(x){all.cols[[x]]["prot",]})
	deut.cols <- lapply(1:length(all.cols),FUN=function(x){all.cols[[x]]["deut",]})
	trit.cols <- lapply(1:length(all.cols),FUN=function(x){all.cols[[x]]["trit",]})
	pairsdiffs.norm <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(norm.cols[[A[1]]], norm.cols[[A[2]]])}))
	pairsdiffs.prot <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(prot.cols[[A[1]]], prot.cols[[A[2]]])}))
	pairsdiffs.deut <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(deut.cols[[A[1]]], deut.cols[[A[2]]])}))
	pairsdiffs.trit <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(trit.cols[[A[1]]], trit.cols[[A[2]]])}))
	difftest <- c(pairsdiffs.norm)
	if(!is.null(cbspace)){
		if("prot" %in% cbspace){
			difftest <- c(difftest,pairsdiffs.prot)
		}
		if("deut" %in% cbspace){
			difftest <- c(difftest,pairsdiffs.deut)
		}
		if("trit" %in% cbspace){
			difftest <- c(difftest,pairsdiffs.trit)
		}
	}
	# difftest <- c(pairsdiffs.norm,pairsdiffs.prot,pairsdiffs.deut,pairsdiffs.trit)
	# difftest  <- c(pairsdiffs.norm,pairsdiffs.prot,pairsdiffs.deut)
	# difftest  <- c(pairsdiffs.norm,pairsdiffs.deut)
	# 
	# colors.temp       <- sample(funky100,n)
	# colors.temp.pset  <- pset(colors.temp)
	# pairs.colors.temp <- colors.temp.pset[which(lengths(colors.temp.pset)==2)]
	# pairs.colorsdiffs <- unlist(lapply(X=pairs.colors.temp,FUN=function(x){diffcol(x[1],x[2])}))
	max.tries     <- iter
	min.tries     <- 10
	try.i         <- 1
	colors.i      <- list(); length(colors.i) <- max.tries
	colors.i[[1]] <- colors.temp
	mindiffs      <- min(difftest)
	while((mindiffs[try.i] < thresh & try.i <= max.tries) | try.i < min.tries){
		try.i <- try.i + 1
		if(try.i > max.tries + 1){
			stop("")
		}
	#	colors.temp       <- sample(funky100,n)
	#	colors.temp.pset  <- pset(colors.temp)
	#	colors.temp <- lapply(X=1:n,FUN=function(x){sample(0:255,3)})
	#	#pairs.colors.temp <- colors.temp.pset[which(lengths(colors.temp.pset)==2)]
	#	#pairs.colorsdiffs <- unlist(lapply(X=pairs.colors.temp,FUN=function(x){diffcol(x[1],x[2])}))
	#	pairs.colorsdiffs <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(cols.list[[A[1]]], cols.list[[A[2]]])}))
	#	colors.i[[try.i]] <- do.call(rbind,colors.temp)
	#	min.i     <- min(pairs.colorsdiffs)
	#	mindiffs  <- c(mindiffs,min.i)
		norm.cols   <- lapply(X=1:n,FUN=function(x){sample(0:255,3)})
		colors.temp <- do.call(rbind,norm.cols)
		all.cols  <- lapply(X=norm.cols,FUN=rgb.cb)
		prot.cols <- lapply(1:length(all.cols),FUN=function(x){all.cols[[x]]["prot",]})
		deut.cols <- lapply(1:length(all.cols),FUN=function(x){all.cols[[x]]["deut",]})
		trit.cols <- lapply(1:length(all.cols),FUN=function(x){all.cols[[x]]["trit",]})
		pairsdiffs.norm <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(norm.cols[[A[1]]], norm.cols[[A[2]]])}))
		pairsdiffs.prot <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(prot.cols[[A[1]]], prot.cols[[A[2]]])}))
		pairsdiffs.deut <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(deut.cols[[A[1]]], deut.cols[[A[2]]])}))
		pairsdiffs.trit <- unlist(lapply(X=1:length(pset.index),FUN=function(x){A=pset.index[[x]]; diffcol(trit.cols[[A[1]]], trit.cols[[A[2]]])}))
		# difftest  <- c(pairsdiffs.norm,pairsdiffs.prot,pairsdiffs.deut,pairsdiffs.trit)
		# difftest  <- c(pairsdiffs.norm,pairsdiffs.prot,pairsdiffs.deut)
		# difftest  <- c(pairsdiffs.norm,pairsdiffs.deut)
		difftest <- c(pairsdiffs.norm)
		if(!is.null(cbspace)){
			if("prot" %in% cbspace){
				difftest <- c(difftest,pairsdiffs.prot)
			}
			if("deut" %in% cbspace){
				difftest <- c(difftest,pairsdiffs.deut)
			}
			if("trit" %in% cbspace){
				difftest <- c(difftest,pairsdiffs.trit)
			}
		}
		min.i     <- min(difftest)
		mindiffs  <- c(mindiffs,min.i)
		colors.i[[try.i]] <- colors.temp
	}
	if(try.i>1){
		result.rgb <- colors.i[[which(mindiffs == max(mindiffs) )[1]]]
		colnames(result.rgb) <- c("r","g","b")
		rownames(result.rgb) <- paste0("col",1:n)
	} else {
		result.rgb <- colors.temp
		colnames(result.rgb) <- c("r","g","b")
		rownames(result.rgb) <- paste0("col",1:n)
	}
	# print(try.i); print(max(mindiffs))
	result.hex <- vapply(X=1:nrow(result.rgb),FUN=function(x){rgb(result.rgb[x,1], result.rgb[x,2], result.rgb[x,3],maxColorValue=255)},FUN.VALUE="char")
	result.hex
}




