#' @title Run DAPC and Plot
#' 
#' Run adagenet DAPC analyses and generate plots of BIC vs number of clusters, baplot of alpha optimum number of principle components at each K, admixture for each K, assignments at each K.
#' If sample coordinates are supplied, this function interpolates cluster membership probabilities on a map, using functions from tess3r package.
#' 
#' @param vcf Character string with path to vcf file containing snp data
#' @param kmax Number indicating the maximum number of clusters to evaluate. Default is 40, which is converted using kmax = min(kmax, number of individuals-1)
#' @param coords Optional character string with path to a table with longitude and latitude of individuals in the vcf file, or a matrix or data frame with longitude and latitude columns. Default is NULL, in which case membership probabilities are not interpolated onto a map.
#' @param reps Number indicating the number of replicates of 'find.clusters'. Default 100.
#' @param include.out Character vector indicating which type of files should be included as output. Default is c(".pdf",".Qlog",".BIClog").
#' @param probs.out NULL or a character string with location where to write a table containing the membership probabilities for the best K and alpha-optimized number of PCAs.
#' @param save.as Character string with where to save the output PDF with plots of results. Default is NULL.
#' @return A list of plots.
#' @export run_DAPC
run_DAPC <- function(vcf, kmax=40, coords=NULL, reps=100,probs.out=NULL,save.as=NULL,include.out=c(".pdf",".Qlog",".BIClog")){
	if(is.null(save.as)){
		save.as <- file.path(getwd(),"result_DAPC.pdf")
	}
	if(file.exists(save.as)){
		stop("Output file already exists. Use a different name for 'save.as' argument.")
	}
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
	
	### Defining colors to use
	if(max.clusters <= 15){
			myCols          <- goodcolors2(n=max.clusters)
		}
		if(max.clusters>15){
			myCols          <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=max.clusters-15))
	}
	
	#### find.clusters
	grp          <- adegenet::find.clusters(genind, max.n.clust=max.clusters,n.pca=max.clusters,choose.n.clust=F)
	grp.list     <- list(); length(grp.list) <- reps
	# par(mar=c(3.5,4,3,2.1))
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
	### Construct data frame holding BIC scores for each replicate of each K
	BIC.df      <- data.frame(BIC=unname(unlist(c(BIC.mat))),K=rep(Krange,reps),replicate=rep(1:reps,each=length(Krange)))
	### save a copy of BIC scores
	if(".BIClog" %in% include.out){
		write.table(x=BIC.df,file=paste0(tools::file_path_sans_ext(save.as),".BIClog"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	}
	#BIC.df <- BIC.df[order(BIC.df.temp[,"Kval"]),]
	#mode(BIC.df$Kval) <- "character"
	BIC.df$K <- factor(BIC.df$K, levels=c(1:nrow(BIC.df)))
	BICPlot  <- ggplot2::ggplot(data=BIC.df,ggplot2::aes(x=K, y=BIC)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("BIC (",reps," replicates of find.clusters) vs. number of clusters (K)"), x="Number of ancestral populations", y = "BIC") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
	### Checking if the max BIC of the next K value is less than the min BIC of the previous K value.
	#range.BIC.mat     <- do.call(rbind,lapply(X=1:nrow(BIC.mat),FUN=function(x){range(BIC.mat[x,],na.rm=TRUE)}))
	#BIC.is.nonoverlap <- NULL
	#BIC.is.reduced    <- NULL
	#for(i in 2:nrow(range.BIC.mat)){
	#	BIC.is.nonoverlap     <- c(BIC.is.nonoverlap,range.BIC.mat[i,1] > range.BIC.mat[(i-1),2] | range.BIC.mat[i,2] < range.BIC.mat[(i-1),1])
	#	BIC.is.reduced     <- c(BIC.is.reduced,range.BIC.mat[i,2] < range.BIC.mat[(i-1),1])
	#}
	# If max BIC for K=2 is better than BIC K=1 and if some K values are not better (all BIC lower) than K-1, then find the first K value in which K+1 is not better.
	#if(any(!BIC.is.reduced) & BIC.is.reduced[1]){
	#	Kbest.criteria2 <- which(!BIC.is.reduced)[1]
	#} else {
	#	Kbest.criteria2 <- 1
	#}
	### Which K value (for K>=2) yields the least variable BIC scores.
	#BICvK.variation <- apply(X=BIC.mat,MARGIN=1,FUN=var)
	#KminVarBIC      <- which(BICvK.variation==min(BICvK.variation[-1]))
	#Kbest.criteria3 <- KminVarBIC[1]
	### Criteria 4: t-tests for BIC of each pairwise adjacent K
	# for(i in 2:nrow(BIC.mat)){
	# 	if(BICvK.variation[Kbest.criteria3]==0){
	# 		Kbest.criteria4 <- NULL
	# 		break
	# 	}
	# 	t.test.i <- t.test(BIC.mat[i-1,],BIC.mat[i,])
	# 	pval.i   <- t.test.i$p.value
	# 	stat.i   <- t.test.i$statistic
	# 	if(pval.i < 0.05 & stat.i > 0){
	# 		next
	# 	} else {
	# 		if(i==nrow(BIC.mat)){
	# 			Kbest.criteria4 <- NULL
	# 		} else{
	# 			Kbest.criteria4 <- (i-1)
	# 			break
	# 		}
	# 	}
	# }
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
#	par(mar=c(4,4,3,2.1))
	# names(best.npca) <- 2:max.clusters
	best.npca.df      <- data.frame(K=2:max.clusters,best.npca=best.npca)
	best.npca.df$K    <- factor(best.npca.df$K)
	grp.plot2         <- ggplot2::ggplot(data=best.npca.df, ggplot2::aes(x=K,y=best.npca)) + ggplot2::geom_bar(stat="identity",fill="lightgray") + ggplot2::labs(title= "alpha optimized # of PCs vs. number of clusters", x="Number of clusters", y = "Alpha optimized number of principle components to retain") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
	#### Creating lists to be filled furing the loop
	dapc.pcabest.list <- list(); length(dapc.pcabest.list) <- max.clusters-1
	admixturePlot  <- list(); length(admixturePlot)        <- max.clusters-1
	scatterPlot    <- list(); length(scatterPlot)          <- max.clusters-1
	da.densityPlot    <- list(); length(da.densityPlot)    <- max.clusters-1
	da.biPlot         <- list(); length(da.biPlot)         <- max.clusters-1
	pca.densityPlot    <- list(); length(pca.densityPlot)  <- max.clusters-1
	pca.biPlot         <- list(); length(pca.biPlot)       <- max.clusters-1
	assignmentPlot <- list(); length(assignmentPlot)       <- max.clusters-1
	posterior.list <- list(); length(posterior.list)       <- max.clusters-1
	mapplot        <- list(); length(mapplot)              <- max.clusters-1
	#q.df           <- NULL
	q.df           <- list(); length(q.df)         <- max.clusters-1
	dapc.df        <- list(); length(dapc.df)      <- max.clusters-1
	for(K in 2:max.clusters){
		i=(K-1)
		dapc.pcabest.K      <- adegenet::dapc(genind, grp.mat[,i],n.pca=best.npca[i],n.da=5)
		dapc.pcabest.list[[i]] <- dapc.pcabest.K
		posterior           <- dapc.pcabest.K$posterior
		q.matrix            <- posterior
		posterior.list[[i]] <- posterior
		posterior.df        <- data.frame(indv=rep(rownames(posterior),ncol(posterior)), pop=rep(colnames(posterior),each=nrow(posterior)), assignment=c(posterior),K=K)
		q.df[[i]]           <- posterior.df
		##### ggplot scatterplots of discriminant functions
		###
		scatterPlot.i       <- ggscatter.dapc(dapc.pcabest.K,col=myCols,legend=F,cstar=1,cpoint=4,label=T)
		scatterPlot[[i]]    <- scatterPlot.i
		### density plots of discriminant functions
		density.da.list.i <- list(); length(density.da.list.i) <- dapc.pcabest.K$n.da
		for(z in 1:dapc.pcabest.K$n.da){
			density.da.list.i[[z]] <- ggscatter.dapc(dapc.pcabest.K,xax=z,yax=z,col=myCols,legend=F,show.title=F,hideperimeter=T)
		}
		### biplots of discriminant functions
		if(dapc.pcabest.K$n.da>1){
			da.pairs.i     <- pset(x=1:dapc.pcabest.K$n.da,min.length=2,max.length=2)
			biplots.da.list.i <- list(); length(biplots.da.list.i) <- length(da.pairs.i)
			for(z in length(da.pairs.i)){
				da.pairs.i.z <- da.pairs.i[[z]]
				biplots.da.list.i[[z]] <- ggscatter.dapc(dapc.pcabest.K,xax=da.pairs.i.z[1],yax=da.pairs.i.z[2],col=myCols,legend=F,cstar=1,cpoint=4,label=F,show.title=F,hideperimeter=T)
			}
		} else {
			biplots.da.list.i <- NULL
		}
		##### ggplot density plots of principle components
		density.pca.list.i <- list(); length(density.pca.list.i) <- dapc.pcabest.K$n.pca
		for(z in 1:dapc.pcabest.K$n.pca){
			density.pca.list.i[[z]] <- ggscatter.dapc(dapc.pcabest.K,vartype="pc",xax=z,yax=z,col=myCols,legend=F,show.title=F,hideperimeter=T)
		}
		### ggplot biplots of principle components
		if(dapc.pcabest.K$n.pca>1){
			pca.pairs.i     <- pset(x=1:dapc.pcabest.K$n.pca,min.length=2,max.length=2)
			biplots.pca.list.i <- list(); length(biplots.pca.list.i) <- length(pca.pairs.i)
			for(z in 1:length(pca.pairs.i)){
				pca.pairs.i.z           <- pca.pairs.i[[z]]
				biplots.pca.list.i[[z]] <- ggscatter.dapc(dapc.pcabest.K,vartype="pc",xax=pca.pairs.i.z[1],yax=pca.pairs.i.z[2],col=myCols,legend=F,cstar=0,cpoint=4,label=F,hideperimeter=T,show.title=F)
			}
		} else {
			biplots.pca.list.i <- NULL
		}
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Membership Probability",y="",fill="Cluster",title=paste0("K = ",K,"; PCs retained = ",best.npca[i])) + ggplot2::scale_fill_manual(values=myCols[1:K])
		admixturePlot[[i]]  <- posterior.gg
		#indv.KmaxPosterior <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))})
		indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels             <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignment.K       <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K,"; PCs retained = ", best.npca[i]), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
		assignmentPlot[[i]]  <- assignment.K
		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols[1:K], 9)
			tess3r.qmat     <- suppressWarnings(tess3r::as.qmatrix(q.matrix))
			coords.mat      <- as.matrix(coords)
			mapplot.i       <- tess3r::ggtess3Q(tess3r.qmat,coords.mat, interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE,map.polygon=world_sp)
			mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
		}
		da.densityPlot[[i]]  <- density.da.list.i
#		da.biPlot[[i]]       <- biplots.da.list.i
		pca.densityPlot[[i]] <- density.pca.list.i
#		pca.biPlot[[i]]      <- biplots.pca.list.i
	}
	
	da.density.arranged      <- dapc.plot.arrange(da.densityPlot)
	pca.densityPlot.arranged <- dapc.plot.arrange(pca.densityPlot,variable="PC")
	### This part wont work. Need to modify function or make new function for biplots such that a page is generated for the combinations DF or PC for a particular K.
#	da.biplot.arranged       <- dapc.plot.arrange(da.biPlot)
#	pca.biPlot.arranged      <- dapc.plot.arrange(pca.biPlot)
	scatterPlot.arranged     <- lapply(1:length(scatterPlot),FUN=function(x){gridExtra::arrangeGrob(scatterPlot[[x]],top=paste0("K=",(x+1)))})
	#scatterPlot.arranged2    <- gridExtra::arrangeGrob(scatterPlot.arranged,nrow=4)
	 #gridExtra::arrangeGrob(scatterPlot,top=col1.names[z])
	########
	## Density plots of discriminant functions
	### For each K, the number of plots (i.e., number of discriminant functions
	#numplots.da.density <- lengths(da.densityPlot)
	#da.max <- max(numplots.da.density)
	#layout.mat          <- matrix(data=NA,nrow=length(numplots.da.density), ncol=da.max)
	## lapply(X=1:length(numplots.da.density),FUN=function(x){ c(1:numplots.da.density[x],rep(0,max(numplots.da.density)-numplots.da.density[x]))})
	#for(i in 1:length(numplots.da.density)){
	#	if(numplots.da.density[i]==0){
	#		next
	#	} else {
	#		layout.mat[i,1:numplots.da.density[i]] <- rep(1,numplots.da.density[i])
	#	}
	#}
	#ent.update <- which(c(layout.mat)!=0)
	#layout.mat[ent.update] <- 1:length(ent.update)
	#### flatten lists of lists
	#da.densityPlot <- do.call(c, da.densityPlot)
	#### convert list of ggplots into list of grobs
	#da.density.grobs       <- lapply(da.densityPlot,FUN=ggplot2::ggplotGrob)
	#### arrange the list of grobs into a gtable
	#da.density.arranged    <- gridExtra::arrangeGrob(grobs=da.density.grobs,layout_matrix=layout.mat)

	q.df    <- do.call(rbind,q.df)
	#dapc.df <- do.call(rbind,dapc.df)
	if(".Qlog" %in% include.out){
		write.table(x=q.df,file=paste0(tools::file_path_sans_ext(save.as),".Qlog"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	}
	#if(bestK>1){
	#	posterior.bestK <- posterior.list[[bestK-1]]
	#	colnames(posterior.bestK) <- paste0("K",colnames(posterior.bestK))
	#	if(!is.null(probs.out)){
	#		write.table(posterior.bestK,file=probs.out,quote=F,col.names=T,row.names=T)
	#	}
	#} else {
	#	posterior.bestK <- matrix(data=rep(1,numind),ncol=1)
	#	rownames(posterior.bestK) <- samplenames
	#}
	#dev.off()
	if(!is.null(coords)){
		result <- c(list(BICPlot,grp.plot2),scatterPlot,admixturePlot,assignmentPlot,mapplot)
	} else {
		result <- c(list(BICPlot,grp.plot2),scatterPlot,admixturePlot,assignmentPlot)
	}
	#if(!is.null(save.as)){
	if(".pdf" %in% include.out){
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
#' Note: I may update this to not require vcftools.
#' 
#' @param vcf Character string with path to input vcf.
#' @param vcftools.path Character string with path to the vcftools executable.
#' @param out Character string where to write output vcf.
#' @param indv.keep Character string with names of individuals to keep. Default is NULL (all individuals kept).
#' @param min.n Minimum number of non-missing alleles required to keep a site. Default = 4. If set to "all", then no sites with any missing data are removed (after first filtering individuals if indv.keep is non-NULL).
#' @param min.n0 Minimum number of individuals required to have at least one copy of the major allele to keep a site. Default = 2.
#' @param min.n1 Minimum number of individuals required to have at least one copy of the minor allele to keep a site. Default = 1.
#' @param which.site Character string indicating the method for choosing a site to keep for each locus (or chromosome). Default = "best", which is considered the one with the least missing data, or the first among sites tied for least missing data. Other options are "all.passing", which retains all sites (positions) that pass variation filters (min.n, min.0.n.0, min.1.n), "first" (first site kept at each locus), or "random".
#' @return List with [[1]] path to vcftools, [[2]] dataframe with input and output values for VCF filepaths, number of loci (chromosomes), sites (positions), and individuals (samples).
#' @export vcf_getSNP
vcf_getSNP      <- function(vcftools.path,vcf,out,indv.keep=NULL,which.site="best",min.n=4,min.n0=2,min.n1=1){
	vcf.obj     <- vcfR::read.vcfR(vcf)
	samplenames <- colnames(vcf.obj@gt)[-1]
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
	gt.mat <- gsub(":.+","",vcf.obj@gt[,-1])
	if(!is.null(indv.keep)){
		if(is(indv.keep,"character") & length(indv.keep)==1){
			if(file.exists(indv.keep)){
				#indv.keep <- unname(unlist(read.table(indv.keep)))
				indv.keep <- suppressWarnings(readLines(indv.keep))
			}
		}
		### Check that all names in indv.keep actually exist in the vcf, and then update gt.mat to only include individuals that pass
		if(all(indv.keep %in% colnames(gt.mat))){
			gt.mat <- gt.mat[,indv.keep,drop=F]
		} else {
			stop(paste(paste(setdiff(indv.keep,colnames(gt.mat)),collapse=","), "not in VCF"))
		}
	}
	# If min.n = "all", update min.all to equal the maximum possible number of alleles at a site.
	if(min.n=="all"){
		test.sample <- unlist(gt.mat)[!is.na(unlist(gt.mat))][1]
		ploidy <- length(unlist(strsplit(test.sample,split="/",fixed=T)))
		min.n  <- ploidy*ncol(gt.mat)
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
		stop("no sites pass filtering criteria")
		#stop("no sites in which the minor allele occurs in more than 1 individual")
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
	chrom.pos.filtered.df  <- chrom.pos.df.filtered2[,c("CHROM","POS")]
	### Need to test if any individuals have only missing data; if TRUE, stop and suggest that these individuals be removed or that a different filtering scheme be used.
	gt.mat2      <- gt.mat[which(fix.mat[,"CHROM"] %in% chrom.pos.filtered.df[,"CHROM"] & fix.mat[,"POS"] %in% chrom.pos.filtered.df[,"POS"]),]
	sitesPerIndv <- unlist(apply(X=gt.mat2,MARGIN=2,FUN=function(x){length(grep("/",x,fixed=TRUE))}))
	if(any(sitesPerIndv==0)){
		noDataIndvs <- names(which(sitesPerIndv==0))
		stop(paste("After filtering,",noDataIndvs,"individuals without data."))
	}
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
	InOut.df <- data.frame(input=c(vcf,length(unique(fix.mat[,"CHROM"])),nrow(vcf.obj@gt),length(samplenames)),output=c(out,length(unique(chrom.pos.filtered.df[,"CHROM"])),nrow(chrom.pos.filtered.df),ncol(gt.mat)))
	rownames(InOut.df) <- c("VCF_filepath","N_loci","N_sites","N_samples")
	colnames(InOut.df) <- c("Input","Output")
	return(list(vcftools.path=vcftools.path,InputOutput=InOut.df))
	#return(c(vcftools.path=vcftools.path, vcf=vcf, out=out, nloci.out=length(unique(chrom.pos.filtered.df[,"CHROM"])), nsites.out=nrow(chrom.pos.filtered.df),numind.out=ncol(gt.mat)))
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
#'
#' vcf_getSNP(vcftools.path="/panfs/pfs.local/home/j926w878/programs/vcftools_0.1.13/bin/vcftools",vcf="/panfs/pfs.local/home/j926w878/programs/easySFS/popfiles/Ahaetulla_snps.vcf",indv.keep="/panfs/pfs.local/home/j926w878/programs/easySFS/popfiles/indv_keep_Ahaetulla-prasina_AllLocalities_GoodData.txt",out="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Ahaetulla-prasina_AllPops_AllSNPs.vcf",which.site="all.passing")
#'
#' vcf_getSNP(vcftools.path="/panfs/pfs.local/home/j926w878/programs/vcftools_0.1.13/bin/vcftools",vcf="/panfs/pfs.local/home/j926w878/programs/easySFS/popfiles/Ahaetulla_snps.vcf",indv.keep="/panfs/pfs.local/home/j926w878/programs/easySFS/popfiles/indv_keep_Ahaetulla-prasina_AllLocalities_GoodData.txt",out="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Ahaetulla-prasina_AllPops_FirstSNP.vcf",which.site="first")
#'
#' vcf_getSNP(vcftools.path="/panfs/pfs.local/home/j926w878/programs/vcftools_0.1.13/bin/vcftools",vcf="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Ahaetulla-prasina_AllPops_AllSNPs.vcf",indv.keep="/panfs/pfs.local/home/j926w878/programs/easySFS/popfiles/indv_keep_Ahaetulla-prasina_Luzon.txt",out="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Ahaetulla-prasina_Luzon_BestSNP.vcf",which.site="best")
#' 
#' vcf_getSNP(vcftools.path="/panfs/pfs.local/home/j926w878/programs/vcftools_0.1.13/bin/vcftools",vcf="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Ahaetulla-prasina_AllPops_AllSNPs.vcf",indv.keep="/panfs/pfs.local/home/j926w878/programs/easySFS/popfiles/indv_keep_Ahaetulla-prasina_Luzon.txt",out="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Ahaetulla-prasina_Luzon_NoMissingData.vcf",which.site="all.passing",min.n = (19*2))
#' 
#' vcf_getSNP(vcftools.path="/panfs/pfs.local/home/j926w878/programs/vcftools_0.1.13/bin/vcftools",vcf="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Calamaria-gervaisii_AllPops_AllSNPs.vcf",indv.keep=NULL,out="/panfs/pfs.local/home/j926w878/work/ddRAD/snps_goodData/Calamaria-gervaisii_AllPops_BestSNP.vcf",which.site="best")

#####

#' @title Function to arrange plots of density for each DF or PC and each K
#' 
#' 
#' 
#' @param x 
#' @param variable Either "DF" or "PC"
#' @return A gtable object
#' @export dapc.plot.arrange
dapc.plot.arrange <- function(x,variable="DF"){
	numplots   <- lengths(x)
	stat.max   <- max(numplots)
	layout.mat0 <- matrix(data=NA,nrow=length(numplots), ncol=stat.max)
	for(i in 1:length(numplots)){
		if(numplots[i]==0){
			next
		} else {
			layout.mat0[i,1:numplots[i]] <- rep(1,numplots[i])
		}
	}
	vals       <- c(t(layout.mat0))
	ent.update <- which(vals!=0)
	vals[ent.update] <- 1:length(ent.update)
	gg.list    <- do.call(c, x)
	grobs.list <- lapply(gg.list, FUN=ggplot2::ggplotGrob)
		layout.mat <- matrix(data=1:length(vals),nrow=length(numplots), ncol=stat.max,byrow=TRUE)
		kmax=nrow(layout.mat)+1
		stat.max <- ncol(layout.mat)
		col1.vals  <- c(layout.mat[,1])[-1]
		row1.vals  <- c(layout.mat[1,])[-1]
		m1n1.names <- c("K=2",paste0(variable,"1"))
		if(stat.max>1){
			col1.names <- paste0(variable,2:stat.max)
		} else {
			col1.names <- NULL
		}
		if(kmax>2){
			row1.names <- paste0("K=",3:kmax)
		} else {
			row1.names <- NULL
		}
		grobsTable.list <- list(); length(grobsTable.list) <- length(vals)
		for(i in 1:length(vals)){
			if(is.na(vals[i])){
				grob.i          <- grid::rectGrob(gp=grid::gpar(col=NA))
			} else {
				grob.i          <- grobs.list[[vals[i]]]
			}
			if(i==1){
				grobTable.i <- gridExtra::arrangeGrob(grob.i,left=m1n1.names[1],top=m1n1.names[2])
			}
			if(i %in% col1.vals){
				z <- which(col1.vals %in% i)
				grobTable.i <- gridExtra::arrangeGrob(grob.i,left=row1.names[z])
			}
			if(i %in% row1.vals){
				z <- which(row1.vals %in% i)
				grobTable.i <- gridExtra::arrangeGrob(grob.i,top=col1.names[z])
			}
			if(!(i %in% c(1,col1.vals,row1.vals))){
				grobTable.i <- gridExtra::arrangeGrob(grob.i)
			}
			grobsTable.list[[i]] <- grobTable.i
		}
		grobs.arranged <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=layout.mat)
		grobs.arranged
}


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
#' Function to return the powerset for a numeric or character vector of set elements, with options for filtering sets by setlength.
#' Changing min.length or max.length arguments from default values produces a subset of the powerset.
#' 
#' @param x Numerical or character vector
#' @param min.length Minimum length of a set to return the set. Default = 0, therefore the NULL set is included by default.
#' @param max.length Maximum length of a set to return the set. Default = NULL, meaning any set at least as large as min.length will be returned.
#' @return List with powerset of x, or a subset of the powerset if min.length != 0 and/or max.length != NULL.
#' @export pset
pset <- function (x,min.length=0,max.length=NULL){
	m   <- length(x)
	out <- list(x[c()])
	for (i in seq_along(x)) {
		out <- c(out, lapply(X=out[lengths(out) < m], FUN=function(y){c(y,x[i])}))
	}
	if(min.length>0){
		out <- out[lengths(out)>=min.length]
	}
	if(!is.null(max.length)){
		out <- out[lengths(out)<=max.length]
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


#' @title Return my favorite set of colors to use for clusters in admixture plots, for a particular value of K.
#' 
#' For each K (2-15), a particular result from the 'goodcolors' function that I like.
#' 
#' @param n Number of colors to include in the output set. Equal to K if using this function for population structure analyses that produce admixture plots or population assignment probabilities.
#' @param plot.palette Whether or not to plot the color palette in the plotting window. Default FALSE.
#' @export goodcolors2
goodcolors2 <- function(n,plot.palette=FALSE){
	if(n==1 | n > 15){
		stop("'n' must be an integer >=2 and <=15")
	}
	cols2  <- c("#659DDE","#9A0C2B") # "#8BBC6F","#4B04E0"
	cols3  <- c("#3DF1C9","#492E50","#FADD56")  # "#F1529B" "#FEF51B" "#2C27E8"
	cols4  <- c("#009CC8","#F48E00","#4300E3","#9AE569") # c("#67F57A","#D7A003","#6E8C48","#C41D1B") # "#256D9E" "#EFFB0A" "#E41A3C" "#F98ACB"
	cols5  <- c("#789D95","#EDCBD2","#3F04C7","#1E9A00","#B95AEC") # "#6A96F1" "#CBE589" "#BA0803" "#7CDBE0" "#DEA14C"
	cols6  <- c("#BF3C2C","#6CD57D","#7D075C","#D7AB41","#E9D4F7","#481EE0") # c("#37175A","#DB8A83","#E91A00","#AA569D","#4385FD","#9AE8CF") # 
	cols7  <- c("#424B13","#47E5F4","#511D49","#3F9876","#182D84","#AB5CC4","#90A930") # c("#382E69","#F7C665","#1DFC96","#200311","#9DCD2D","#3202E0","#851B20")
	cols8  <- c("#EA7250","#9DDF5B","#B74CE0","#687235","#021025","#1BB1A1","#3F8187","#C7B7CE")
	cols9  <- c("#8838BD","#997275","#E24911","#ACEB4F","#1DB756","#611D97","#6C7BCF","#125034","#F4ABD7")
	cols10 <- c("#98ACD2","#6197ED","#7A478F","#986D65","#D9C452","#22553B","#B4993A","#0DF1FA","#E8971B","#36D8B9") # "#51BE1A" "#34A6A2" "#86D0E8" "#BDFD81" "#21594F" "#F1F9D6" "#678534" "#6B49C8" "#57191C" "#EC906B"
	cols11 <- c("#33E8BB","#C69A01","#799525","#030A34","#091E57","#4ED671","#A4282F","#2F08BA","#5BD80A","#CD0B50","#B9C1F5") # "#5BBDFD" "#200042" "#616B8F" "#578B7C" "#4D59C3" "#B0FE9F" "#0B335C" "#7F780B" "#C55BD7" "#9F2502" "#D9EF1C"
	cols12 <- c("#FE42DD","#68B206","#FFD4E9","#4FBF72","#567D96","#94245B","#226826","#B8CC94","#1752F7","#8B8C5B","#29A6CF","#2ADEC9") 
	cols13 <- c("#0ACCF8","#9F9D46","#2FCA60","#3D31C8","#50087B","#A17D78","#C0B27E","#EA6EA5","#B21D00","#31E701","#AC1479","#0B3D98","#5DBB15") # "#FCED80" "#336AB6" "#8F3A8C" "#FA9304" "#3545B1" "#83FDF4" "#B6575B" "#232E32" "#AF9746" "#4CA612" "#ED7B66" "#C3C600" "#47B9F5" # "#5B45C8" "#63F469" "#F3B4C8" "#D9EFE8" "#C74B77" "#F1A655" "#813402" "#3D27C4" "#CCD427" "#27AD47" "#13117F" "#97753D" "#09C1B2"
	cols14 <- c("#A9D12A","#1D5DB8","#F18249","#44EB9D","#3702AA","#0AF6FF","#48063C","#8DAFB0","#343E3D","#ED20FD","#D1B67D","#472B20","#F1A4CC","#CE1B4F") # "#FFA357" "#7AAD72" "#2D5AB4" "#E41B0B" "#51ECF5" "#007866" "#D43CB8" "#34E590" "#DC9F06" "#91827F" "#830638" "#B07527" "#F6BCC0" "#2B63EC" #### Change first color: "#62D236" "#FCED80" "#336AB6" "#8F3A8C" "#FA9304" "#3545B1" "#83FDF4" "#B6575B" "#232E32" "#AF9746" "#4CA612" "#ED7B66" "#C3C600" "#47B9F5"
	cols15 <- c("#810890","#EFD814","#2D793E","#796ED5","#699866","#343E3D","#F6233D","#877D0E","#C9BAD5","#8BB0B1","#00C3C0","#522AA9","#21FF61","#EA8316","#F5C387")
	cols16 <- NULL
	cols17 <- NULL
	cols18 <- NULL
	cols19 <- NULL
	cols20 <- NULL
	cols21 <- NULL
	cols22 <- NULL
	cols23 <- NULL
	cols24 <- NULL
	cols25 <- NULL
	colorList <- list(cols2 ,cols3 ,cols4 ,cols5 ,cols6 ,cols7 ,cols8 ,cols9 ,cols10,cols11,cols12,cols13,cols14,cols15,cols16,cols17,cols18,cols19,cols20,cols21,cols22,cols23,cols24,cols25)
	result    <- colorList[[n-1]]
	if(plot.palette){
		dummydata        <- rep(1,n)
		names(dummydata) <- 1:n
		bp <- barplot(dummydata,col=result,axes=FALSE,ylim=c(0,1.25))
		mtext(side=1,text="Color Palette for Clusters",line=2.5)
		text(x=c(bp), y=1.12, labels=result,srt=90)
	}
	result
}


#' @title Make ggplot with DAPC scatterplot
#' 
#' THIS FUNCTION IS LARGELY BASED ON THE FUNCTION 'scatter.dapc' FROM ADEGENET
#' This function produces plots very similar to those produced by the adegenet function 'scatter.dapc', except that the object returned is a ggplot. This is useful when making lists of plots.
#' SOME FEATURES STILL IN PROGRESS:
#' - ability to add screeplot of discriminant functions or PCAs.
#' 
#' @param x Object of class DAPC (adegenet package).
#' @param xax Which descriminant function to plot on the x axis.
#' @param yax Which descriminant function to plot on the y axis. Default 2. Ignored if only one discriminant function exists.
#' @param vartype Character string indicating which variables to plot. Either "df" for discriminant functions or "pc" for principle components. Default "df".
#' @param bg Color to use for the background of the plot. Default "white".
#' @param grp Object of class 'factor' with length equal to the number of individuals, which indicates individual assignments to clusters. By default posterior assignments are extracted from 'x'.
#' @param col Vector with color to use for each cluster.
#' @param pch Either one number indicating pch code of symbol to use for all points (see 'pch' in graphical parameters), or, a numeric vector with pch codes to each for clusters . Default 20.
#' @param cpoint Number with cex size to use for points.
#' @param solid Number between 0 and 1 indicating the level of transparency to use for colors; 0 = fully transparent; 1 = fully opaque; default = 0.7.
#' @param legend A logical indicating whether a legend for group colours should added to the plot. Default FALSE.
#' @param onedim.filled Logical indicating, when only one discriminant function is to be plotted, whether or not density plots should be filled or unfilled with the colors indicated by 'col'. Default TRUE.
#' @param hideperimeter Logical indicating whether or not to hide the x and y axes and labels (the perimeter area of the plottting area). Default FALSE.
#' @param show.title Logical indicating whether or not to include a plot title. Default TRUE. Overriden if hideperimeter=TRUE.
#' @param addaxes Logical indicating if reference axes should be drawn at x=a and y=b, with a and b supplied by the 'origin' argument. Default TRUE.
#' @param ltyaxes Line weight to use for edges of the minimum spanning tree linking the groups. Default 0.5.
#' @param lwdaxes Line type to use for edges of the minimum spanning tree linking the groups. Default 2 (dashed).
#' @param cellipse A positive coefficient for the inertia ellipse size. Default 1.5. Setting to zero removes ellipses.
#' @param cstar A number greater than 0 defining the length of the star size (i.e., the lines radiating from the center of clusters towards individuals). Default 1. Setting to zero removes the star lines; setting =1 connects points to cluster mean; setting > 1 extends lines through their points.
#' @param mstree A logical indicating whether a minimum spanning tree linking the groups and based on the squared distances between the groups inside the entire space should added to the plot (TRUE), or not (FALSE). Default FALSE.
#' @param lwd Line weight to use for edges of the minimum spanning tree linking the groups. Default 0.25.
#' @param lty Line type to use for edges of the minimum spanning tree linking the groups. Default 1 (solid).
#' @param segcol Color to use for edges of the minimum spanning tree linking the groups Default "black".
#' @param label Logical indicating whether or not clusters should be labeled in discriminant function biplots. Default TRUE.
#' @param clabel A number indicating the size of labels. Default 1.
#' @param axesell NOT YET IMPLEMENTED. A logical value indicating whether the ellipse axes should be drawn. Default FALSE.
#' @param txt.leg NOT YET IMPLEMENTED. Labels to use for clusters in the legend. 
#' @param scree.da NOT YET IMPLEMENTED. Logical indicated whether or not a screeplot of the discriminant functions should be included. Default TRUE.
#' @param scree.pca  NOT YET IMPLEMENTED. Logical indicated whether or not a screeplot of the PCs should be included. Default FALSE.
#' @param posi.da NOT YET IMPLEMENTED. The position of the screeplot of discriminant functions. Can match any combination of "top/bottom" and "left/right", or a set of x/y coordinates stored as a list (locator can be used). Default = "bottomright".
#' @param posi.pca NOT YET IMPLEMENTED. The position of the screeplot of discriminant functions. Can match any combination of "top/bottom" and "left/right", or a set of x/y coordinates stored as a list (locator can be used). Default = "bottomleft".
#' @param posi.leg NOT YET IMPLEMENTED. The position of the legend holding group colors. Can match any combination of "top/bottom" and "left/right", or a set of x/y coordinates stored as a list (locator can be used). Default = "topright".
#' @param bg.inset NOT YET IMPLEMENTED. Default "white".
#' @param ratio.da NOT YET IMPLEMENTED. Default 0.25
#' @param ratio.pca NOT YET IMPLEMENTED. Default 0.25
#' @param inset.da NOT YET IMPLEMENTED. Default 0.02
#' @param inset.pca NOT YET IMPLEMENTED. Default 0.02
#' @param inset.solid NOT YET IMPLEMENTED. Default 0.5
#' @param cleg NOT YET IMPLEMENTED. Size factor for the legend. Default 1.
#' @param xlim NOT YET IMPLEMENTED. Default NULL.
#' @param ylim NOT YET IMPLEMENTED. Default NULL.
#' @param grid NOT YET IMPLEMENTED. Whether or not to include a grid in the background. Default FALSE.
#' @param cgrid NOT YET IMPLEMENTED. A number used with par("cex")* cgrid to specify the mesh of the grid.
#' @param origin NOT YET IMPLEMENTED. Location of the origin.
#' @param include.origin NOT YET IMPLEMENTED. A logical value indicating whether the point "origin" should be belonged to the graph space. Default TRUE.
#' @param sub NOT YET IMPLEMENTED. A string of characters to be inserted as legend. Default "".
#' @param csub NOT YET IMPLEMENTED. Number specifying text size for 'sub'.
#' @param possub NOT YET IMPLEMENTED. The position of the subtitle ("topleft", "topright", "bottomleft", "bottomright"). Default "bottomleft".
#' @param pixmap NOT YET IMPLEMENTED. An object 'pixmap' displayed in the map background. Default NULL.
#' @param contour NOT YET IMPLEMENTED. A data frame with 4 columns to plot the contour of the map : each row gives a segment (x1,y1,x2,y2). Default NULL. 
#' @param area NOT YET IMPLEMENTED. A data frame of class 'area' to plot a set of surface units in contour. Default NULL.
#' @param label.inds NOT YET IMPLEMENTED. Default NULL. Named list of arguments passed to the orditorp function. This will label individual points witout overlapping. Arguments x and display are hardcoded and should not be specified by user.
#' @param new.pred NOT YET IMPLEMENTED. An optional list, as returned by the predict method for dapc objects; if provided, the individuals with unknown groups are added at the bottom of the plot. To visualize these individuals only, specify only.grp="unknown".
#' @param only.grp NOT YET IMPLEMENTED. Character vector indicating which groups should be displayed. Values should match values of x$grp. If NULL, all results are displayed.
#' @return A ggplot object
#' @export ggscatter.dapc
ggscatter.dapc <- function (x, xax = 1, yax = 2, vartype="df", grp = x$grp , cpoint=2, col = adegenet::seasun(length(levels(grp))), txt.leg = levels(grp), label = TRUE, pch = 20, solid = 0.7, hideperimeter=FALSE, show.title=TRUE,scree.da = TRUE, scree.pca = FALSE, posi.da = "bottomright", posi.pca = "bottomleft",bg="white", bg.inset = "white", ratio.da = 0.25, ratio.pca = 0.25, inset.da = 0.02, inset.pca = 0.02, inset.solid = 0.5, onedim.filled = TRUE, mstree = FALSE, lwd = 0.25, lty = 1, segcol = "black", legend = FALSE, posi.leg = "topright", cleg = 1, cstar = 1, cellipse = 1.5, axesell = FALSE, clabel = 1, xlim = NULL, ylim = NULL, grid = FALSE, addaxes = TRUE, ltyaxes=2, lwdaxes=0.5, origin = c(0,0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, label.inds = NULL, new.pred=NULL){
	if(vartype=="df"){
		ind.vals <- x$ind.coord
		varname  <- "Discriminant function"
	} else {
		if(vartype=="pc"){
			ind.vals <- x$tab
			varname  <- "Principle component"
			mstree   <- FALSE
		} else {
			stop("'vartype' must be either 'df' or 'pc'")
		}
	}
	### Logical indicating if only one dimension retained
	ONEDIM     <- xax == yax | ncol(ind.vals) == 1
	simple.col <- transp(col[1:length(levels(grp))],solid)
	col        <- transp(rep(col, length(levels(grp))), solid)
	pch        <- rep(pch, length(levels(grp)))
	bg.inset   <- transp(bg.inset, inset.solid)
	### Posterior assignments of individuals to groups
	if (is.null(grp)) {
		grp <- x$grp
	}
	### Further checking/updating if there is only one discriminant function
	if (is.null(xax) || is.null(yax)) {
		## a number followed by L indiciates that the class should be an integer
		xax    <- 1L
		yax    <- ifelse(ncol(ind.vals) == 1L, 1L, 2L)
		ONEDIM <- TRUE
	}
	### Things to do when more than one PC exists.
	if(!ONEDIM){
		coords.df  <- data.frame(x.coords=ind.vals[, xax],y.coords=ind.vals[, yax],Cluster=grp)
		xlim       <- c(-max(abs(coords.df[,"x.coords"])),max(abs(coords.df[,"x.coords"])))
		ylim       <- c(-max(abs(coords.df[,"y.coords"])),max(abs(coords.df[,"y.coords"])))
		### Creating columns to hold x and y mean of group that each individual belongs too
		unique.clusters <- levels(grp)
		for(z in unique.clusters){
			rows.temp <- which(coords.df$Cluster==z)
			coords.df[rows.temp,"grp.center.x"] <- mean(coords.df[rows.temp,"x.coords"])
			coords.df[rows.temp,"grp.center.y"] <- mean(coords.df[rows.temp,"y.coords"])
		}
		xy3.df           <- do.call(rbind,lapply(1:nrow(coords.df),FUN=function(z){newpoint(p0=coords.df[z,c("grp.center.x","grp.center.y")], p1=coords.df[z,c("x.coords","y.coords")], c=cstar)}))
		coords.df[,"x3"] <- xy3.df[,1]
		coords.df[,"y3"] <- xy3.df[,2]
		### blank plotting area
		ggscatter.tempA      <- ggplot2::ggplot(coords.df, ggplot2::aes(x=x.coords, y=y.coords,color=Cluster,shape=Cluster,fill=Cluster)) + ggplot2::scale_x_continuous(name=paste(varname,xax)) + ggplot2::scale_y_continuous(name=paste(varname,yax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = bg))  + ggplot2::stat_ellipse(color="white") + ggplot2::theme_classic() + ggplot2::geom_blank()
		# Includes box around plotting area
		ggscatter.tempB      <- ggscatter.tempA + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) 
		# Add reference lines (axes) at x=0 and y=0
		if(addaxes){
			ggscatter.tempB  <- ggscatter.tempB + ggplot2::geom_vline(ggplot2::aes(xintercept=0),linetype=ltyaxes,size=lwdaxes,color="lightgray") + ggplot2::geom_hline(ggplot2::aes(yintercept=0),linetype=ltyaxes,size=lwdaxes,color="lightgray")
		}
		# Hide axis ticks and labels
		if(hideperimeter){
			ggscatter.tempC  <- ggscatter.tempB + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
		} else {
			ggscatter.tempC  <- ggscatter.tempB
		}
		### Adding the points. The scale to use for colors of points was defined in ggscatter.tempA so no need to redefine colors here.
		ggscatter.temp0      <- ggscatter.tempC + ggplot2::geom_point(size=cpoint,show.legend=TRUE) + ggplot2::scale_color_manual(values=col) + ggplot2::scale_shape_manual(values=pch) #+ ggplot2::scale_fill_manual(values=col,fill=col)
		# Hide or show legend (guide)
		if(!legend){
			ggscatter.temp1  <- ggscatter.temp0 + ggplot2::theme(legend.position = "none")
		} else {
			ggscatter.temp1  <- ggscatter.temp0 + ggplot2::theme(legend.position = c(0.98,0.98), legend.justification=c("right","top"), legend.background = ggplot2::element_rect(fill="white", size=0.25, linetype="solid")) + ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(fill=simple.col,color="black",size=12,shape=22))) # + ggplot2::scale_fill_manual(values=simple.col) 
		}
		# Add ellipses around clusters
		if(cellipse > 0){
			ggscatter.temp2  <- ggscatter.temp1 + ggplot2::stat_ellipse(level=(cellipse*0.43),type="norm",show.legend=FALSE)
		} else {
			ggscatter.temp2  <- ggscatter.temp1
		}
		# Add 'star' lines from each cluster mean to the coordinates of individuals in the cluster.
		if(cstar > 0){
			ggscatter.temp3 <- ggscatter.temp2 + ggplot2::geom_segment(data = coords.df, ggplot2::aes(x = x3, y = y3, xend = grp.center.x, yend = grp.center.y, color = Cluster),show.legend=FALSE)
		} else {
			ggscatter.temp3 <- ggscatter.temp2
		}
		if(mstree){
			meanposi <- apply(x$tab, 2, tapply, grp, mean)
			axes     <- c(xax, yax)
			D        <- dist(meanposi)^2
			tre      <- ade4::mstree(D)
			x0       <- unname(x$grp.coord[tre[, 1], axes[1]])
			y0       <- unname(x$grp.coord[tre[, 1], axes[2]])
			x1       <- unname(x$grp.coord[tre[, 2], axes[1]])
			y1       <- unname(x$grp.coord[tre[, 2], axes[2]])
			tree.df  <- data.frame(xA=x0,yA=y0,xB=x1,yB=y1)
			tree.mat <- cbind(x0,y0,x1,y1)
			coords.df[,"tree.x0"] <- x0
			coords.df[,"tree.y0"] <- y0
			coords.df[,"tree.x1"] <- x1
			coords.df[,"tree.y1"] <- y1
			ggscatter.temp4 <- ggscatter.temp3 + ggplot2::geom_segment(data = coords.df, ggplot2::aes(x = tree.x0, y = tree.y0, xend = tree.x1, yend = tree.y1), color = segcol, size=lwd,linetype=lty,show.legend=FALSE)
		} else {
			ggscatter.temp4 <- ggscatter.temp3
		}
		#if(!is.null(label)){
		if(label){
			ggscatter.temp5 <- ggscatter.temp4 + ggplot2::geom_label(data=coords.df,ggplot2::aes(x=grp.center.x,y=grp.center.y,label=Cluster),fill="white",size=(clabel*3),show.legend=FALSE)
		} else {
			ggscatter.temp5 <- ggscatter.temp4
		}
		# return(ggscatter.temp5)
	} else {# If only one PC
		scree.da <- FALSE
		if(ncol(ind.vals) == 1) {
			pcLab <- 1
		} else {
			pcLab <- xax
		}
		### Data frame with coordinates of individuals and the posterior assignment of individuals to groups (clusters)
		coords.df  <- data.frame(coords=ind.vals[, pcLab],Cluster=grp)
		### Apply the density function to the individual coordinates for individuals in each group
		ldens <- tapply(X=ind.vals[, pcLab], INDEX=grp, FUN=density)
		### The actual x (coorinates) and y (density) values that are to be plotted
		allx  <- unlist(lapply(ldens, function(e) e$x))
		ally  <- unlist(lapply(ldens, function(e) e$y))
		## defining locations for x-axis ticks
		xat0 <- seq(from=round(min(allx)),to=round(max(allx)))
		xat  <- xat0[xat0/2 == round(xat0/2)]
		xpoints <- ind.vals[grp == levels(grp)[i], pcLab]
		ypoints <- rep(0, sum(grp == levels(grp)[i]))
		### ggplot of PC coordinates of groups.
		if(onedim.filled){
			gg.density.temp  <- ggplot2::ggplot(coords.df, ggplot2::aes(x=coords,color=Cluster,fill=Cluster)) + ggplot2::geom_density() + ggplot2::theme_classic() + ggplot2::scale_color_manual(values=col) + ggplot2::scale_fill_manual(values=col) + ggplot2::scale_x_continuous(breaks=xat,name=paste(varname,xax),limits=range(allx))
		} else {
			gg.density.temp  <- ggplot2::ggplot(coords.df, ggplot2::aes(x=coords,color=Cluster,fill=NA)) + ggplot2::geom_density() + ggplot2::theme_classic() + ggplot2::scale_color_manual(values=col) + ggplot2::scale_fill_manual(values=NA) + ggplot2::scale_x_continuous(breaks=xat,name=paste(varname,xax),limits=range(allx))
		}
		gg.density0       <- gg.density.temp + ggplot2::ylab("Density") + ggplot2::geom_text(ggplot2::aes(x=coords,y=rep(0,length(coords)),label=rep("|",length(coords))),show.legend=FALSE) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1), axis.line=ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size=12), axis.text.y = ggplot2::element_text(size=12))
		### Add a title
		if(show.title){
			gg.density1 <- gg.density0 + ggplot2::labs(title=paste0("K=",K))
		} else {
			gg.density1 <- gg.density0
		}
		### Remove legend if legend = FALSE
		if(!legend){
			gg.density2 <- gg.density1 + ggplot2::theme(legend.position = "none")
		} else {
			gg.density2 <- gg.density1 + ggplot2::theme(legend.position = c(0.98,0.98), legend.justification=c("right","top"))
		}
		# Hide axis ticks and labels
		if(hideperimeter){
			gg.density3  <- gg.density2 + ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())
		} else {
			gg.density3  <- gg.density2
		}

	}
	if(ONEDIM){
		return(gg.density3)
	} else {
		return(ggscatter.temp5)
	}
}


#' @title Coordinates of point on ray
#' 
#' Used in the function ggscatter.dapc to control segments radiating from cluster centers towards individual's pca coordinates.
#' Returns x,y coordinates of a point (p2) on a ray with origin at p0(x0,y0) and another point at p1(x1,y1), such that the distance from p0 to p2 = c*(distance from p0 to p1).
#' 
#' @param p0 Numberical vector with coordinates c(x0,y0) defining the origin of a ray
#' @param p1 Numberical vector with coordinates c(x1,y1) defining a point on the ray other than the origin.
#' @param c A number equal to (distance p0 to p2)/(distance p0 to p1)
#' @return Numberical vector with coordinates c(x2,y2) of p2.
#' @export newpoint
newpoint <-function(p0,p1,c){
	v  <- p1-p0
	d  <- sqrt(sum(v^2))
	u  <- v/d
	d2 <- d*c
	p2 <- p0 + (d2*u)
	p2
}






