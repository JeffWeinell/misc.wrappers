#' @title Run DAPC and Plot
#' 
#' Run adagenet DAPC analyses and generate plots of BIC vs number of clusters, baplot of alpha optimum number of principle components at each K, admixture for each K, assignments at each K.
#' If sample coordinates are supplied, this function interpolates cluster membership probabilities on a map, using functions from tess3r package.
#' 
#' @param x 'vcfR' object (see package::vcfR) or character string with path to file containing snp data, with data format specified by 'format' argument (currently only VCF files can be used).
#' @param format Character string indicating the format (filetype) of the data. Currently only "VCF" is allowed, but other types may be added.
#' @param kmax Number indicating the maximum number of clusters to evaluate. Default is 40, which is converted using kmax = min(kmax, number of individuals-1)
#' @param coords Optional character string with path to a table with longitude and latitude of individuals in the vcf file, or a matrix or data frame with longitude and latitude columns. Default is NULL, in which case membership probabilities are not interpolated onto a map.
#' @param reps Number indicating the number of replicates of 'find.clusters'. Default 100.
#' @param include.out Character vector indicating which type of files should be included as output. Default is c(".pdf",".Qlog",".BIClog").
#' @param plot.components FALSE. This is still in development.
#' @param save.as Character string with where to save the output PDF with plots of results. Default is NULL.
#' @param save.in Character string with path to directory where output files should be saved.
#' @param overwrite Logical indicating whether or not to allow new output files to overwrite existing ones. Default FALSE.
#' @return A list of plots.
#' @export run_DAPC
run_DAPC <- function(x, format="VCF", kmax=40, coords=NULL, samplenames=NULL, reps=100, save.as=NULL,save.in=NULL, plot.components=FALSE, include.out=c(".pdf",".Qlog",".BIClog"),overwrite=FALSE){
	debug <- FALSE
	if(is.null(save.as)){
		save.as <- file.path(getwd(),"result_DAPC.pdf")
	}
	if(file.exists(save.as)){
		stop("Output file already exists. Use a different name for 'save.as' argument.")
	}

	if(is.null(save.in)){
		save.in <- getwd()
	}
	if(!is.null(include.out)){
		if(".pdf" %in% include.out){
			save.as.pdf <- file.path(save.in,"result_dapc.pdf")
		} else {
			save.as.pdf <- NULL
		}
		if(".Qlog" %in% include.out){
			save.as.Qlog <- file.path(save.in,"result_dapc.Qlog")
		} else {
			save.as.Qlog <- NULL
		}
		if(".BIClog" %in% include.out){
			save.as.BICLog <- file.path(save.in,"result_dapc.BICLog")
		} else {
			save.as.BICLog <- NULL
		}
	} else {
		save.as.pdf <- save.as.Qlog <- save.as.BICLog  <- NULL
	}
	if(!overwrite){
		files.to.check <- c(save.as.pdf,save.as.Qlog,save.as.BICLog)
		if(!is.null(files.to.check)){
			if(any(files.to.check %in% save.in)){
				stop("One or more output files already exist in directory indicated by 'save.in'. Choose a different output directory or change 'overwrite' to TRUE")
			}
		}
	}
	





	if(format=="VCF" | is(x,"vcfR")){
		if(is(x,"vcfR")){
			vcf.obj <- vcf <- x
		} else {
			vcf <- x
			vcf.obj     <- vcfR::read.vcfR(vcf,verbose=F,checkFile=F)
		}
		if(is.null(samplenames)){
			samplenames <- colnames(vcf.obj@gt)[-1]
		}
		genind      <- suppressWarnings(vcfR::vcfR2genind(vcf.obj))
		#numind      <- (dim(attributes(vcf.obj)[["gt"]])[2])-1
	} else {
		stop("Currently, 'format' must be 'VCF'")
	}
	numind      <- length(samplenames)
	label.size  <- min((288/numind),7)
	if(debug) message("step 0")
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
		current_sf    <- suppressWarnings(suppressMessages(sf::st_crop(world_sf,xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max)))
		current.gg.sf <- ggplot2::geom_sf(data=current_sf,color = "black", fill = NA)
		#world.gg.sf   <- ggplot2::geom_sf(data=world_sf,colour = "black", fill = NA)
		#message("potential replacement for 'current_sf' and 'current.gg.sf' is on the next line")
		current.gg.sf2 <- ggplot2::geom_sf(data=world_sf) + ggplot2::coord_sf(xlim=c(x.min,x.max),ylim=c(y.min,ymax=y.max))
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
			myCols  <- goodcolors2(n=15)[1:max.clusters]
		}
		if(max.clusters>15){
			myCols  <- c(goodcolors2(n=15), sample(adegenet::funky(100), size=max.clusters-15))
	}
	if(debug) message("step 1")
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
	BIC.df      <- data.frame(BIC=unname(unlist(c(BIC.mat))), K=rep(Krange,reps), replicate=rep(1:reps,each=length(Krange)))
	### save a copy of BIC scores
	if(".BIClog" %in% include.out){
		write.table(x=BIC.df, file=paste0(tools::file_path_sans_ext(save.as),".BIClog"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	}
	#BIC.df <- BIC.df[order(BIC.df.temp[,"Kval"]),]
	#mode(BIC.df$Kval) <- "character"
	BIC.df$K <- factor(BIC.df$K, levels=c(1:nrow(BIC.df)))
	BICPlot  <- ggplot2::ggplot(data=BIC.df,ggplot2::aes(x=K, y=BIC)) + ggplot2::geom_boxplot(fill='lightgray', outlier.colour="black", outlier.shape=16,outlier.size=2, notch=FALSE) + ggplot2::theme_classic() + ggplot2::labs(title= paste0("BIC (",reps," replicates of find.clusters) vs. number of clusters (K)"), x="Number of ancestral populations", y = "BIC") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) # + ggplot2::geom_vline(xintercept=bestK, linetype=2, color="black", size=0.25)
	if(debug) message("step 1.5")
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
	if(debug) message("step 2")
	##### Plot 2: BIC vs. K when number of PCs retained = alpha optimized 
#	par(mar=c(4,4,3,2.1))
	# names(best.npca) <- 2:max.clusters
	best.npca.df      <- data.frame(K=2:max.clusters,best.npca=best.npca)
	best.npca.df$K    <- factor(best.npca.df$K)
	grp.plot2         <- ggplot2::ggplot(data=best.npca.df, ggplot2::aes(x=K,y=best.npca)) + ggplot2::geom_bar(stat="identity",fill="lightgray") + ggplot2::labs(title= "alpha optimized # of PCs vs. number of clusters", x="Number of clusters", y = "Alpha optimized number of principle components to retain") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

	#### Creating lists to be filled
	empty.set <- list(); length(empty.set) <- max.clusters-1
	dapc.pcabest.list <- admixturePlot <- da.densityPlot <- da.biPlot <- pca.densityPlot <- pca.biPlot <- da.psets <- da.layout.mat<- pca.psets <- pca.layout.mat <- assignmentPlot <- posterior.list <- mapplot <- q.df <- dapc.df <- empty.set
	
	dapc.list <- empty.set
	for(K in 2:max.clusters){
		i=K-1
		dapc.list[[(i)]] <- adegenet::dapc(genind, grp.mat[,i],n.pca=best.npca[i],n.da=5)
	}
	names(dapc.list) <- paste0("K",c(2:max.clusters))
	grp.assignments.mat <- do.call(cbind,lapply(dapc.list,FUN=function(x){c(x$assign)}))
	rownames(grp.assignments.mat) <- samplenames
	colnames(grp.assignments.mat) <- paste0("K",c(2:max.clusters))
	grp.minsizes <- apply(grp.assignments.mat,MARGIN=2,FUN=function(x){min(table(x))})
	### The values of K for which density plots can be made because every cluster has at least two individuals
	if(any((grp.minsizes > 1))){
		Ks.showdensity <- unname(which((grp.minsizes > 1))+1)
	}
	posterior.list <- lapply(dapc.list,function(x){x$posterior})
	# Number of discriminant functions retained at each value of K
	Ks.n.da        <- sapply(dapc.list,function(x){x$n.da})
	# Number of principle components retained at each value of K
	Ks.n.pca       <- sapply(dapc.list,function(x){x$n.pca})
	# Data frame holding assignment probability of each individual in each cluster at each K.
	q.df           <- do.call(rbind,lapply(X=1:length(dapc.list),FUN=function(x){data.frame(indv=rep(rownames(posterior.list[[x]]),ncol(posterior.list[[x]])), pop=rep(colnames(posterior.list[[x]]),each=nrow(posterior.list[[x]])), assignment=c(posterior.list[[x]]),K=(x+1))}))
	if(".Qlog" %in% include.out){
		write.table(x=q.df,file=paste0(tools::file_path_sans_ext(save.as),".Qlog"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
	}
	
	#####
	## Density Plots section still produces errors under some scenarios.
	if(plot.components){
		#######
		## Density plots
		#######
		## List with two data frames with combinations of K and DF, or K and PC, respectively, for which density can be plotted because more than two individuals per cluster.
		plottable   <- plottable.dapc(dapc.list)
		if(!is.null(plottable)){
			K.plottable <- as.numeric(rownames(plottable$DF))
			#layout.da <- plottable$DF
			#layout.da[plottable$DF]  <- which(plottable$DF)
			#layout.da[!plottable$DF] <- NA
			### Creating plottable matrices for biplots
			# if(any(Ks.n.da>2)){
			# 	pairs <- lapply(X=apply(X=plottable$DF[apply(plottable$DF,MARGIN=1,FUN=function(x){length(which(x))>2}),], MARGIN=1, FUN=which),FUN=pset,min.length=2,max.length=2)
			#}
			#layout.pc <- plottable$PC
			#layout.pc[plottable$PC]  <- which(plottable$PC)
			#layout.pc[!plottable$PC] <- NA
			#indexmat.da <- matrix(1:length(plottable$DF),ncol=ncol(plottable$DF))
			#indexmat.pc <- matrix(1:length(plottable$PC),ncol=ncol(plottable$PC))
			## plots should be added by column, from bottom to top of the gtable.
			indexmat0.da <- matrix(1:length(plottable$DF),ncol=ncol(plottable$DF))
			indexmat0.pc <- matrix(1:length(plottable$PC),ncol=ncol(plottable$PC))
			if(nrow(indexmat0.da)>1){
				indexmat.da <- apply(indexmat0.da,2,rev)
			} else {
				indexmat.da <- indexmat0.da
			}
			if(nrow(indexmat0.pc)>1){
				indexmat.pc <- apply(indexmat0.pc,2,rev)
			} else {
				indexmat.pc <- indexmat0.pc
			}
			### Matrices to hold row and colnames of the da and pca gtables that will be produced below.
			emptymat.da  <- matrix(data=rep("",length(plottable$DF)),ncol=ncol(plottable$DF))
			emptymat.pca <- matrix(data=rep("",length(plottable$PC)),ncol=ncol(plottable$PC))
			left.mat.da  <- bottom.mat.da  <- right.mat.da  <- top.mat.da  <-emptymat.da
			left.mat.pca <- bottom.mat.pca <- right.mat.pca <- top.mat.pca <-emptymat.pca
			left.mat.da[,1]  <- left.mat.pca[,1] <- rev(paste0("K",K.plottable))
			right.mat.da[,ncol(plottable$DF)] <- right.mat.pca[,ncol(plottable$PC)] <- rev(paste0("K",K.plottable))
			bottom.mat.da[nrow(plottable$DF),]  <- colnames(plottable$DF)
			bottom.mat.pca[nrow(plottable$PC),] <- colnames(plottable$PC)
			top.mat.da[1,]  <- colnames(plottable$DF)
			top.mat.pca[1,] <- colnames(plottable$PC)
			
			#### Generate DF density plots.
			df.plots.list <- list(); length(df.plots.list) <- length(plottable$DF)
			ctr <- 0
			for(j in 1:ncol(plottable$DF)){
				for(i in 1:nrow(plottable$DF)){
					ctr <- ctr+1
					if(!plottable$DF[i,j] | any(table(dapc.list[[(K.plottable[i]-1)]]$assign) < 2)){
					#	next
						plot.ij  <- grid::rectGrob(gp=grid::gpar(col=NA))
						
					} else {
					#	ctr <- ctr+1
						plot.ij  <- ggplot2::ggplotGrob(ggscatter.dapc(dapc.list[[(K.plottable[i]-1)]], vartype="df", xax=j, yax=j, col=myCols, legend=F, show.title=F, hideperimeter=T))
					}
					df.plots.list[[ctr]]        <- plot.ij
					names(df.plots.list)[[ctr]] <- paste0("K",K.plottable[i],".DF",j)
				}
			}
			da.arranged0   <- lapply(X=1:length(df.plots.list),FUN=function(x){gridExtra::arrangeGrob(df.plots.list[[x]],left=left.mat.da[indexmat.da[x]],right=right.mat.da[indexmat.da[x]],bottom=bottom.mat.da[indexmat.da[x]],top=top.mat.da[indexmat.da[x]])})
			da.arranged    <- gridExtra::arrangeGrob(grobs=da.arranged0,layout_matrix=indexmat.da,respect=TRUE)
			vp             <- grid::viewport(height=grid::unit(0.9,"npc"),width=grid::unit(0.9,"npc"))
			pdf(file=paste0(tools::file_path_sans_ext(save.as),"_densityPlots_DF.pdf"), height=(nrow(indexmat.da)*3),width=(ncol(indexmat.da)*3))
			grid::grid.draw(da.arranged)
			dev.off()
			
			#### Generate PC density plots.
			pc.plots.list <- list(); length(pc.plots.list) <- length(plottable$PC)
			ctr <- 0
			for(j in 1:ncol(plottable$PC)){
				for(i in 1:nrow(plottable$PC) ){
					ctr <- ctr+1
					if(!plottable$PC[i,j]){
					#	next
						plot.ij  <- grid::rectGrob(gp=grid::gpar(col=NA))
						
					} else {
					#	ctr <- ctr+1
						plot.ij  <- ggplot2::ggplotGrob(ggscatter.dapc(dapc.list[[(K.plottable[i]-1)]], vartype="pc", xax=j, yax=j, col=myCols, legend=F, show.title=F, hideperimeter=T))
					}
					pc.plots.list[[ctr]]        <- plot.ij
					names(pc.plots.list)[[ctr]] <- paste0("K",K.plottable[i],".PC",j)
				}
			}
			pc.arranged0   <- lapply(X=1:length(pc.plots.list),FUN=function(x){gridExtra::arrangeGrob(pc.plots.list[[x]],left=left.mat.pca[indexmat.pc[x]],right=right.mat.pca[indexmat.pc[x]],bottom=bottom.mat.pca[indexmat.pc[x]],top=top.mat.pca[indexmat.pc[x]])})
			pc.arranged    <- gridExtra::arrangeGrob(grobs=pc.arranged0,layout_matrix=indexmat.pc,respect=TRUE)
			vp             <- grid::viewport(height=grid::unit(0.9,"npc"),width=grid::unit(0.9,"npc"))
			pdf(file=paste0(tools::file_path_sans_ext(save.as),"_densityPlots_PC.pdf"), height=(nrow(indexmat.pc)*3),width=(ncol(indexmat.pc)*3))
			grid::grid.draw(pc.arranged)
			dev.off()
		
		}
		#######
		## Bi-component plots
		#######
		if(any(Ks.n.da>2)){
			#distmats.da <- lapply(apply(plottable$DF[apply(plottable$DF,MARGIN=1,FUN=function(x){length(which(x))>2}),],MARGIN=1,FUN=which), FUN=function(y){dist(y)})
			distmats.da <- lapply(X=Ks.n.da[which(Ks.n.da>1)], FUN=function(x){dist(table(paste0("DF",c(1:x))))})
			distlabs.da <- lapply(distmats.da,FUN=function(x){gtools::mixedsort(attributes(x)$Labels)})
			#damats0  <- lapply(apply(plottable$DF[apply(plottable$DF,MARGIN=1,FUN=function(x){length(which(x))>2}),],MARGIN=1,FUN=which), FUN=function(y){lower.tri(dist(y))})
			damats0  <- lapply(distmats.da, FUN=lower.tri)
			davals   <- lapply(damats0,FUN=function(x){1:length(which(x))})
			damats1  <- damats0
			for(i in 1:length(damats1)){
				damats1[[i]][damats0[[i]]]  <- davals[[i]]
				damats1[[i]][!damats0[[i]]] <- NA
				damats1[[i]] <- damats1[[i]] # [-c(1),-c(ncol(damats1[[i]]))]
				rownames(damats1[[i]]) <- distlabs.da[[i]] # [-1]
				colnames(damats1[[i]]) <- distlabs.da[[i]] # [-c(lengths(distlabs.da)[i])]
			}
			damats2 <- damats1
			for(i in 1:length(damats2)){
				damats2[[i]][!is.na(damats2[[i]])] <- paste0("bp",damats2[[i]][!is.na(damats2[[i]])],".da")
				damats2[[i]] <- t(damats2[[i]])
			}
			if(any(Ks.n.pca>2)){
				#distmats.pc <- lapply(apply(plottable$PC[apply(plottable$PC,MARGIN=1,FUN=function(x){length(which(x))>2}),],MARGIN=1,FUN=which), FUN=function(y){dist(y)})
				distmats.pc <- lapply(X=Ks.n.pca[which(Ks.n.pca>1)], FUN=function(x){dist(table(paste0("DF",c(1:x))))})
				distlabs.pc <- lapply(distmats.pc,FUN=function(x){gtools::mixedsort(attributes(x)$Labels)})
				#pcmats0  <- lapply(apply(plottable$PC[apply(plottable$PC,MARGIN=1,FUN=function(x){length(which(x))>2}),],MARGIN=1,FUN=which), FUN=function(y){lower.tri(dist(y))})
				pcmats0  <- lapply(distmats.pc,lower.tri)
				pcvals   <- lapply(pcmats0,FUN=function(x){1:length(which(x))})
				pcmats1  <- pcmats0
				for(i in 1:length(pcmats1)){
					pcmats1[[i]][pcmats0[[i]]]  <- pcvals[[i]]
					pcmats1[[i]][!pcmats0[[i]]] <- NA
					pcmats1[[i]] <- pcmats1[[i]]   # [-c(1),-c(ncol(pcmats1[[i]]))]
					rownames(pcmats1[[i]]) <- distlabs.pc[[i]]  #[-1]
					colnames(pcmats1[[i]]) <- distlabs.pc[[i]]  # [-c(lengths(distlabs.pc)[i])]
				}
				pcmats2 <- pcmats1
				for(i in 1:length(pcmats2)){
					pcmats2[[i]][!is.na(pcmats2[[i]])] <- paste0("bp",pcmats2[[i]][!is.na(pcmats2[[i]])],".pc")
				}
				dapc.mats <- pcmats2
				unionKs <- union(names(pcmats2),names(damats2))
				### The value of dapc.mats mostly portrays the desired plotting layout.
				for(i in 1:length(dapc.mats)){
					if(unionKs[i] %in% names(pcmats2) & unionKs[i] %in% names(damats2)){
						j.pca <- which(unionKs[i] == names(pcmats2))
						j.da  <- which(unionKs[i] == names(damats2))
						dapc.mats[[i]][which(is.na(pcmats2[[j.pca]]))[1:length(c(damats2[[j.da]][!lower.tri(damats2[[j.da]])]))]] <- c(damats2[[j.da]][!lower.tri(damats2[[j.da]])])
					}
					if( (unionKs[i] %in% names(pcmats2)) & (!(unionKs[i] %in% names(damats2)))){
						j.pca <- which(unionKs[i] == names(pcmats2))
						dapc.mats[[i]] <- pcmats2[[j.pca]]
					}
					if((!(unionKs[i] %in% names(pcmats2))) & (unionKs[i] %in% names(damats2))){
						j.da  <- which(unionKs[i] == names(damats2))
						dapc.mats[[i]] <- damats2[[j.da]]
					}
					colnames(dapc.mats[[i]]) <- 1:ncol(dapc.mats[[i]])
					# which(is.na(pcmats2[[i]])) [1:length(c(damats2[[i]] [!lower.tri(damats2[[i]])]))]
				}
				## Row m of matrix in da.psets2 contains the dimensions to of the mth da biplot to generate.
				da.psets  <- lapply(X=1:length(Ks.n.da),FUN=function(x){do.call(rbind, pset(x=1:Ks.n.da[x], min.length=2, max.length=2))})
				da.psets2 <- lapply(which(lengths(da.psets)>0), FUN=function(x){da.psets[[x]][order(c(da.psets[[x]][,1]),c( da.psets[[x]][,2])),]})
				names(da.psets2) <- which(lengths(da.psets)>0)+1
				### Row m of matrix in pca.psets2 contains the dimensions to of the mth pca biplot to generate.
				pca.psets <- lapply(X=1:length(Ks.n.pca),FUN=function(x){ do.call(rbind, pset(x=1:Ks.n.pca[x], min.length=2, max.length=2))})
				pca.psets2 <- lapply(which(lengths(pca.psets)>0), FUN=function(x){pca.psets[[x]][order(c(pca.psets[[x]][,1]),c( pca.psets[[x]][,2])),]})
				names(pca.psets2) <- which(lengths(pca.psets)>0)+1
				#### Generate DF biplots.
				df.biplots.list <- list(); length(df.biplots.list) <- length(da.psets2)
				for(i in 1:length(da.psets2)){
					K <- as.numeric(names(da.psets2))[i]
					df.biplots.list[[i]] <- lapply(X=1:nrow(rbind(da.psets2[[i]])), FUN=function(x){ggplot2::ggplotGrob(ggscatter.dapc(dapc.list[[(K-1)]], vartype="df", xax= rbind(da.psets2[[i]]) [x,1], yax=rbind(da.psets2[[i]])[x,2], col=myCols, legend=F, show.title=F, hideperimeter=T,cellipse=0))})
					# Names matching those used in 'dapc.mats'
					names(df.biplots.list[[i]]) <- c(t(damats2[[i]]))[!is.na(c(t(damats2[[i]])))]
				}
				#### Generate PC biplots.
				pc.biplots.list <- list(); length(pc.biplots.list) <- length(pca.psets2)
				for(i in 1:length(pca.psets2)){
					K <- as.numeric(names(pca.psets2))[i]
					pc.biplots.list[[i]] <- lapply(X=1:nrow(rbind(pca.psets2[[i]])), FUN=function(x){ggplot2::ggplotGrob(ggscatter.dapc(dapc.list[[(K-1)]], vartype="pc", xax=rbind(pca.psets2[[i]])[x,1], yax=rbind(pca.psets2[[i]])[x,2], col=myCols, legend=F, show.title=F, hideperimeter=T,cellipse=0))})
					names(pc.biplots.list[[i]]) <- c(pcmats2[[i]])[!is.na(c(pcmats2[[i]]))]
				}
				### Reverse the order of rows so that axis labels (DF and PC components) start at lower left.
				dapc.mats2 <- lapply(dapc.mats,apply,2,rev)
				### Emptry matrices
				empty.bi <- lapply(1:length(dapc.mats2), FUN=function(x){matrix(data=rep("",length(dapc.mats2[[x]])), ncol=ncol(dapc.mats2[[x]]))})
				### Arrange the biplots into a gtable for each K
				Ks.n.da.bi  <- Ks.n.da[which(Ks.n.da>1)]
				Ks.n.pca.bi <- Ks.n.pca[which(Ks.n.pca>1)]
				# initializing lists of matrices that will hold the names along axes
				left.mat.bi <- right.mat.bi <- bottom.mat.bi <- top.mat.bi <- empty.bi
				for(i in 1:length(dapc.mats2)){
					if(unionKs[i] %in% names(pcmats2) & unionKs[i] %in% names(damats2)){
						j.pca <- which(unionKs[i] == names(pcmats2))
						j.da  <- which(unionKs[i] == names(damats2))
						left.mat.bi[[i]][,1]   <- rev(c("",paste0("PC",2:Ks.n.pca.bi[j.pca])))
						right.mat.bi[[i]][,ncol(right.mat.bi[[i]])]  <- rev(c(paste0("DF", 1:(Ks.n.da.bi[j.da]-1)),rep("",(nrow(right.mat.bi[[i]])-(Ks.n.da.bi[j.da]-1)))))
						bottom.mat.bi[[i]][nrow(bottom.mat.bi[[i]]),] <- c(rep("",(ncol(bottom.mat.bi[[i]])- (Ks.n.da.bi[j.da]-1) )),paste0("DF",2:(Ks.n.da.bi[j.da])))
						top.mat.bi[[i]][1,]  <- c(paste0("PC", (1:(Ks.n.pca.bi[j.pca]-1))), rep("",(ncol(top.mat.bi[[i]])-(Ks.n.pca.bi[j.pca]-1))))
					}
					if((unionKs[i] %in% names(pcmats2)) & (!(unionKs[i] %in% names(damats2)))){
						j.pca <- which(unionKs[i] == names(pcmats2))
						left.mat.bi[[i]][,1]   <- rev(c("",paste0("PC",2:Ks.n.pca.bi[j.pca])))
						top.mat.bi[[i]][1,]  <- c(paste0("PC", (1:(Ks.n.pca.bi[j.pca]-1))), rep("",(ncol(top.mat.bi[[i]])-(Ks.n.pca.bi[j.pca]-1))))
					}
					if((!(unionKs[i] %in% names(pcmats2))) & (unionKs[i] %in% names(damats2))){
						j.da  <- which(unionKs[i] == names(damats2))
						right.mat.bi[[i]][,ncol(right.mat.bi[[i]])]  <- rev(c(paste0("DF", 1:(Ks.n.da.bi[j.da]-1)),rep("",(nrow(right.mat.bi[[i]])-(Ks.n.da.bi[j.da]-1)))))
						bottom.mat.bi[[i]][nrow(bottom.mat.bi[[i]]),] <- c(rep("",(ncol(bottom.mat.bi[[i]])- (Ks.n.da.bi[j.da]-1) )),paste0("DF",2:(Ks.n.da.bi[j.da])))
					}
				}
				### Removes rows and columns that are completely filled with NAs
				dapc.mats3 <- dapc.mats2
				for(i in 1:length(dapc.mats2)){
					#dapc.mats3[[i]] <- unname(dapc.mats2[[i]][!apply(dapc.mats2[[i]],MARGIN=1,FUN=function(x){all(is.na(x))}),!apply(dapc.mats2[[i]],MARGIN=2,FUN=function(x){all(is.na(x))})])
					rowskeep <- !apply(dapc.mats2[[i]],MARGIN=1,FUN=function(x){all(is.na(x))})
					colskeep <- !apply(dapc.mats2[[i]],MARGIN=2,FUN=function(x){all(is.na(x))})
					dapc.mats3[[i]]    <- unname(dapc.mats2[[i]][rowskeep,colskeep])
					left.mat.bi[[i]]   <- left.mat.bi[[i]][rowskeep,colskeep]
					right.mat.bi[[i]]  <- right.mat.bi[[i]][rowskeep,colskeep]
					bottom.mat.bi[[i]] <- bottom.mat.bi[[i]][rowskeep,colskeep]
					top.mat.bi[[i]]    <- top.mat.bi[[i]][rowskeep,colskeep]
				}
				### Index matrix
				index.bi <- lapply(1:length(dapc.mats3), FUN=function(x){apply(matrix(data=c(1:length(dapc.mats3[[x]])) ,ncol=ncol(dapc.mats3[[x]])),2,rev)})
				
				bi.arranged0 <- list(); length(bi.arranged0) <- length(dapc.mats3)
				bi.arranged  <- list(); length(bi.arranged) <- length(dapc.mats3)
				
				### In this loop, i= table of plots for a specific K, and j=the index number in the matrix 'index.bi'.
				for(i in 1:length(dapc.mats3)){
					if(unionKs[i] %in% names(pcmats2)){
						x.pca <- which(unionKs[i] == names(pcmats2))
					} else {
						x.pca <- NULL
					}
					if(unionKs[i] %in% names(damats2)){
						x.da  <- which(unionKs[i] == names(damats2))
					} else {
						x.da  <- NULL
					}
					for(j in 1:length(dapc.mats3[[i]])){
						plot.name <- dapc.mats3[[i]][which(index.bi[[i]]==j)]
						if(is.na(plot.name)){
							plot.ij <- grid::rectGrob(gp=grid::gpar(col=NA))
						} else {
							if(!is.null(x.pca)){
								if(plot.name %in% names(pc.biplots.list[[x.pca]])){
									plot.ij <- pc.biplots.list[[x.pca]][plot.name][[1]]
								} else {
									if(!is.null(x.da)){
										if(plot.name %in% names(df.biplots.list[[x.da]])){
											plot.ij <- df.biplots.list[[x.da]][plot.name][[1]]
										}
									}
								}
							}
						}
						bi.arranged0[[i]][[j]] <- gridExtra::arrangeGrob(plot.ij,left=left.mat.bi[[i]][which(index.bi[[i]]==j)],right=right.mat.bi[[i]][which(index.bi[[i]]==j)],bottom=bottom.mat.bi[[i]][which(index.bi[[i]]==j)],top=top.mat.bi[[i]][which(index.bi[[i]]==j)])
					}
					bi.arranged[[i]] <- gridExtra::arrangeGrob(grobs=bi.arranged0[[i]],layout_matrix=index.bi[[i]],respect=TRUE)
				}
				#pc.biplots.list[[3]]["bp1.pc"][[1]]
				#pc.arranged0[[i]]   <- lapply(X=1:length(dapc.mats3[[i]]),FUN=function(x){gridExtra::arrangeGrob(dapc.mats3[[i]][],left=left.mat.pca[indexmat.pc[x]],right=right.mat.pca[indexmat.pc[x]],bottom=bottom.mat.pca[indexmat.pc[x]],top=top.mat.pca[indexmat.pc[x]])})
				#pc.arranged    <- gridExtra::arrangeGrob(grobs=pc.arranged0,layout_matrix=indexmat.pc,respect=TRUE)
				vp             <- grid::viewport(height=grid::unit(0.95,"npc"),width=grid::unit(0.95,"npc"))
				pdf(file=paste0(tools::file_path_sans_ext(save.as),"_BiPlots.pdf"), height=(max(sapply(dapc.mats3,nrow))*3),width=(max(sapply(dapc.mats3,ncol))*3))
				for(i in 1:length(dapc.mats3)){
					grid::grid.draw(bi.arranged[[i]])
					if(i < length(dapc.mats3)){
						grid::grid.newpage()
					}
				}
				dev.off()
			}
		}
	}
	if(debug) message("step 3")
	for(K in 2:max.clusters){
		if(debug) message(paste0("K=",K," step 3.1"))
		i=(K-1)
		q.matrix           <- posterior.list[[i]]
		rownames(q.matrix) <- samplenames
		colnames(q.matrix) <- paste0("cluster",1:ncol(q.matrix))
		indv.pop            <- apply(X=q.matrix, MARGIN=1, FUN=function(x){which(x==max(x))[1]})
		posterior.df        <- q.df[q.df$K==K,]
		#posterior.df$indv   <- factor(posterior.df$indv, levels = names(sort(posterior.df$assignment)))
		posterior.df$indv  <- factor(posterior.df$indv, levels = names(sort(indv.pop)))
		if(debug) message(cat("\r",paste0("K=",K," step 3.3")))
		if(debug) message(cat("\r",paste0("K=",K," step 3.5")))
		posterior.gg        <- ggplot2::ggplot(posterior.df, ggplot2::aes(fill= pop, x= assignment, y=indv)) + ggplot2::geom_bar(position="stack", stat="identity") + ggplot2::theme_classic() + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank()) + ggplot2::labs(x = "Membership Probability",y="",fill="Cluster",title=paste0("K = ",K,"; PCs retained = ",best.npca[i])) + ggplot2::scale_fill_manual(values=myCols[1:K])
		admixturePlot[[i]]  <- posterior.gg
		if(debug) message(cat("\r",paste0("K=",K," step 3.7")))
		indv.maxPosterior  <- apply(X=q.matrix, MARGIN=1, FUN=function(x){max(x)})
		labels             <- rep("",nrow(posterior.df))
		labels[posterior.df[,"assignment"] %in% indv.maxPosterior] <- "+"
		assignment.K        <- ggplot2::ggplot(data=posterior.df, ggplot2::aes(x= pop, y=indv,fill=assignment)) + ggplot2::geom_tile(color="gray") + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), axis.text.y = ggplot2::element_text(size = label.size), panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), legend.position = "none", ) + ggplot2::labs(title = paste0("K = ",K,"; PCs retained = ", best.npca[i]), x="Clusters", y="") + ggplot2::scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint = 0.5) + ggplot2::geom_text(label=labels)
		assignmentPlot[[i]] <- assignment.K
		if(debug) message(cat("\r",paste0("K=",K," step 3.8")))
		# if(debug) message(cat("\r",paste0("smallest cluster has ",minsize.grp," individuals")))
		if(!is.null(coords)){
			my.palette      <- tess3r::CreatePalette(myCols[1:K], 9)
			tess3r.qmat     <- suppressWarnings(tess3r::as.qmatrix(q.matrix))
			#coords.mat      <- as.matrix(coords)
			coords.mat      <- as.matrix(coords)[,1:2]
			mapplot.i       <- tess3r::ggtess3Q(tess3r.qmat,coords.mat, interpolation.model = tess3r::FieldsKrigModel(10),resolution = c(500,500), col.palette = my.palette, window=c(x.min,x.max,y.min,y.max),background=TRUE,map.polygon=world_sp)
			coords.df <- coords
			#colnames(coords.df[,1:2]) <- c("Lon","Lat")
			names(coords.df)[1:2] <- c("Lon","Lat")
		#	message("Note to self: change the next line to use geom_sf and coord_sf instead of 'current.gg.sf'")
		#	mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + current.gg.sf + ggplot2::geom_point(data = coords.df, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			### Use coord_sf and geom_sf instead of 
			# ggplot2::coord_sf(xlim =,ylim = )
			#mapplot[[i]]    <- mapplot.i + ggplot2::theme_classic() + ggplot2::labs(title=paste0("Ancestry coefficients; K=",K), x="latitude", y="longitude") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::geom_point(data = coords, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
		}
		if(debug) message(cat("\r",paste0("K=",K," step 3.9")))
	} ### END OF LOOP
	
	if(debug) message("step 10.0")
	### Defining the viewport size and creating lists of plots
	vp                  <- grid::viewport(height=grid::unit(0.95,"npc"),width=grid::unit(0.95,"npc"))
	#resultsA.grobs.list <- lapply(list(BICPlot, grp.plot2), FUN=ggplot2::ggplotGrob)
	resultsA            <- lapply(X=lapply(list(BICPlot, grp.plot2), FUN=ggplot2::ggplotGrob), FUN=gridExtra::arrangeGrob, vp=vp)
	if(!is.null(coords)){
		gg.plots.list0 <- c(admixturePlot, assignmentPlot, mapplot)
		gg.plots.list  <- gg.plots.list0[!sapply(gg.plots.list0,is.null)]
		resultsB       <- lapply(X= lapply(c(admixturePlot, assignmentPlot, mapplot), FUN=ggplot2::ggplotGrob), FUN=gridExtra::arrangeGrob,vp=vp)
		
	} else {
		resultsB <- lapply(X=lapply(c(admixturePlot, assignmentPlot), FUN=ggplot2::ggplotGrob), FUN=gridExtra::arrangeGrob,vp=vp)
	}
	results1   <- c(resultsA, resultsB)
	if(debug) message("step 11")
	if(".pdf" %in% include.out){
		pdf(height=6,width=10,file=save.as, onefile=TRUE)
		for(i in 1:length(results1)){
			grid::grid.draw(results1[[i]])
			if(i < length(results1)){
				grid::grid.newpage()
			}
		}
		dev.off()
	}
	results1
}
#' @examples
#' library(misc.wrappers)
#' # Define path to input VCF file containing similated data for 500 SNPs from 50 individuals in three populations.
#' example_vcf_path <- file.path(system.file("extdata", package = "misc.wrappers"),"simK4.vcf.gz")
#' # DAPC analyses for K=2–10 and 30 replicates of find.clusters (adegenet). Saves graphs to "DAPC_example.pdf" in your current directory.
#' run_DAPC(x=example_vcf_path, kmax=10, reps=30, save.as="DAPC_example.pdf", include.out=c(".pdf"))
#' Performs 30 replicates adegenet::find.clusters and plot boxplots of the K vs. BIC;  on the simulated dataset; 

#' @title Arrange DAPC density plots.
#' 
#' Function to arrange plots of density for each DF or PC and each K
#' 
#' @param x A list of ggplots produced by ggscatter.dapc function
#' @param variable Either "DF" or "PC"
#' @param layout.mat Custom layout matrix to use. Default is NULL, in which case the layout will be generated automatically as a function of the number of plots.
#' @param pos.x.labs Where the default labels x axis labels (i.e., column labels) should be drawn. Default = 1 (below table of plots); alternatively 3 (above plot).
#' @param pos.y.labs Where the default labels y axis labels (i.e., row labels) should be drawn. Default = 2 (left side); alternatively 3 (right side).
#' @param row.labels.left Labels to use the left of the first column of plots. Default NULL (no labels).
#' @param col.labels.top Text labels to use above the first row of plots. Default NULL (no labels).
#' @param row.labels.right Text labels to to the right of the last column of plots. Default NULL (no labels).
#' @param col.labels.bottom Text labels to use below the bottom row of plots. Default NULL (no labels).
#' @param use.diag NULL (the default) or numerical vector with values in 1:4, indicating which sides of the gtable (1=bottom, 2=left, 3=top, 4=right) should have labels applied to the diagonal rather than the table margin plots. Ignored if layout.mat is not square. For example, if use.diag=1, then the values of col.labels.bottom will appear below the plots on the diagonal rather than below the plots of the bottom row.
#' @param pad Amount of space between plots, in units of line widths (Default 0.1).
#' @param K Which set of biplots to use. If NULL (the default), the function will attemp to draw all biplots for all K (max 25 plots).
#' @param outer.text A list with length=4, with each entry either NULL (the default) or character string to use as labels below, left, above, and to the right of the arrangement of plots.
#' @param maxMat Numerical vector with length 2 that specifies the maximum number of rows and columns of plots, respectively, in the output gtable. Default = c(7,8).
#' @return A gtable object
#' @export dapc.plot.arrange
dapc.plot.arrange <- function(x,variable="DF",layout.mat=NULL,pos.x.labs=1,pos.y.labs=2,row.labels.left=NULL,col.labels.top=NULL,row.labels.right=NULL,col.labels.bottom=NULL,use.diag=NULL,pad=0.1,K=NULL,outer.text=list(NULL,NULL,NULL,NULL),maxMat=c(7,2)){
	if(length(x)>maxMat[1]){
		x <- x[1:maxMat[1]]
	}
	if(any(lengths(x)>maxMat[2])){
		x <- lapply(x,function(i){keep=min(length(i),maxMat[2]);i[1:keep]})
	}
	numplots   <- lengths(x)
	stat.max   <- max(numplots)
	#layout.mat0 <- matrix(data=NA,nrow=length(numplots), ncol=stat.max)
	#for(i in 1:length(numplots)){
	#	if(numplots[i]==0){
	#		next
	#	} else {
	#		layout.mat0[i,1:numplots[i]] <- rep(1,numplots[i])
	#	}
	#}
	#vals       <- c(t(layout.mat0))
	#ent.update <- which(vals!=0)
	#vals[ent.update] <- 1:length(ent.update)
	vals <- do.call(c,lapply(1:length(numplots),FUN=function(x){c(rep(1, numplots[x]),rep(NA,c(stat.max-numplots)[x]))}))
	vals[which(vals!=0)] <- 1:length(which(vals!=0))
	gg.list     <- do.call(c, x)
	grobs.list  <- lapply(gg.list, FUN=ggplot2::ggplotGrob)
	if(is.null(layout.mat)){
		layout.mat  <- matrix(data=vals,nrow=length(numplots), ncol=stat.max,byrow=TRUE)
	}
	#layout.mat  <- matrix(data=1:length(vals),nrow=length(numplots), ncol=stat.max,byrow=TRUE)
	#layout.mat0 <- matrix(data=1:length(vals),nrow=length(numplots), ncol=stat.max,byrow=TRUE)
	#layout.mat  <- layout.mat0[rev(1:nrow(layout.mat0)),,drop=FALSE]
	kmax        <- nrow(layout.mat)+1
	stat.max    <- ncol(layout.mat)
	#if(FALSE){
	#	col1.vals  <- c(layout.mat[,1])[-1]
	#	row1.vals  <- c(layout.mat[1,])[-1]
	#	m1n1.names <- c("K=2",paste0(variable,"1"))
	#	if(stat.max>1){
	#		col1.names      <- paste0(variable,2:stat.max)
	#		row.labels.left <- row1.names
	#		
	#	} else {
	#		col1.names <- NULL
	#	}
	#	if(kmax>2){
	#		row1.names      <- paste0("K=",3:kmax)
	#		col.labels.top <- col1.names
	#	} else {
	#		row1.names <- NULL
	#	}
	#}
	if(pos.x.labs==1){
		col.labels.bottom  <- paste0(variable,1:stat.max)
	}
	if(pos.x.labs==3){
		col.labels.top  <- paste0(variable,1:stat.max)
	}
	if(pos.y.labs==2){
		row.labels.left <- rev(paste0("K=",2:kmax))
	}
	if(pos.y.labs==4){
		row.labels.right <- paste0("K=",2:kmax)
	}
	nm   <- nrow(layout.mat)
	nn   <- ncol(layout.mat)
	len  <- nm*nn
#	vals <- c(t(layout.mat))
	#numplots <- length(vals)
	if(is.null(row.labels.left)){
		row.labels.left <- rep("",nrow(layout.mat))
	}
	if(is.null(row.labels.right)){
		row.labels.right <- rep("",nrow(layout.mat))
	}
	if(is.null(col.labels.top)){
		col.labels.top <- rep("",ncol(layout.mat))
	}
	if(is.null(col.labels.bottom)){
		col.labels.bottom <- rep("",ncol(layout.mat))
	}
	### character matrix of empty strings
	empty.matrix <- matrix(data="",nrow=nm,ncol=nn)
	### Index assigned to each plot of the table of plots
	index.matrix0 <- matrix(1:len,ncol=nn,byrow=TRUE)
	index.matrix  <- index.matrix0[rev(1:nrow(index.matrix0)),,drop=FALSE]
	### Matrices holding the bottom, left, top, and right labels, respectively, for each plot.
#	if(nm>maxMat[1]){
#		keep.rows    <- rev(rev(1:nm)[1:maxMat[1]])
#		index.matrix <- index.matrix[keep.rows,,drop=FALSE]
#		layout.mat   <- layout.mat[keep.rows,,drop=FALSE]
#		empty.matrix <- empty.matrix[keep.rows,,drop=FALSE]
#		nm  <- maxMat[1]
#	}
#	if(nn>maxMat[2]){
#		keep.cols    <- 1:maxMat[2]
#		index.matrix <- index.matrix[,keep.cols,drop=FALSE]
#		layout.mat   <- layout.mat[,keep.cols,drop=FALSE]
#		empty.matrix <- empty.matrix[,keep.cols,drop=FALSE]
#		nn  <- maxMat[2]
#	}
#	len  <- nm*nn
#	vals <- c(t(layout.mat))
	if(nm>1){
		if(1 %in% use.diag & nm==nn){
			bottom.mat <- empty.matrix
			bottom.mat[diag(index.matrix)] <- col.labels.bottom
		} else {
			bottom.mat <- unname(rbind(matrix(data="",nrow=(nm-1),ncol=nn),col.labels.bottom))
		}
		if(3 %in% use.diag & nm==nn){
			top.mat <- empty.matrix
			top.mat[diag(index.matrix)] <- col.labels.top
		} else {
			top.mat    <- unname(rbind(col.labels.top,matrix(data="",nrow=(nm-1),ncol=nn)))
		}
	} else {
		bottom.mat <- matrix(data=col.labels.bottom,nrow=1)
		top.mat    <- matrix(data=col.labels.top,nrow=1)
	}
	if(nn>1){
		if(2 %in% use.diag & nm==nn){
			left.mat   <- empty.matrix
			left.mat[diag(index.matrix)] <- row.labels.left
		} else {
			left.mat   <- unname(cbind(row.labels.left,matrix(data="",nrow=nm,ncol=(nn-1))))
		}
		if(4 %in% use.diag & nm==nn){
			right.mat  <- empty.matrix
			right.mat[diag(index.matrix)] <- row.labels.right
		} else {
			right.mat  <- unname(cbind(matrix(data="",nrow=nm,ncol=(nn-1)),row.labels.right))
		}
	} else {
		left.mat   <- matrix(data=row.labels.left,ncol=1)
		right.mat  <- matrix(data=row.labels.right,ncol=1)
	}
	grobsTable.list <- list(); length(grobsTable.list) <- length(vals)
	for(i in 1:length(vals)){
		#### Either creating an empty grob or getting a grob from grobs.list
		if(is.na(vals[i])){
			grob.i          <- grid::rectGrob(gp=grid::gpar(col=NA))
		} else {
			grob.i          <- grobs.list[[vals[i]]]
		}
		z <- which(index.matrix == i, arr.ind=TRUE)
		labels.i.list  <- list(bottom.mat[z],left.mat[z],top.mat[z],right.mat[z])
		if(is.null(use.diag)){
			labels.i.list2 <- list(); length(labels.i.list2) <- 4
			for(j in 1:4){
				if(labels.i.list[[j]]!=""){
					labels.i.list2[[j]] <- labels.i.list[[j]]
				}
			}
		} else {
			labels.i.list2 <- labels.i.list
		}
	#	grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=bottom.mat[z],left=left.mat[z],top=top.mat[z],right=right.mat[z])
		grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=labels.i.list2[[1]],left=labels.i.list2[[2]],top=labels.i.list2[[3]],right=labels.i.list2[[4]])
	}
	grobs.arranged0 <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=index.matrix,padding=unit(pad,"line"),respect=TRUE)
	### trim the figure if too many rows
	if(nrow(grobs.arranged0)>maxMat[1]){
		#rows.keep       <- rev(1:nrow(grobs.arranged0))[c(1:maxMat[1])]
		#grobs.arranged0 <- grobs.arranged0[rows.keep,]
		grobs.arranged0 <- grobs.arranged0[c(1:maxMat[1]),]
	}
	### trim the figure if too many columns
	if(ncol(grobs.arranged0)>maxMat[2]){
		grobs.arranged0 <- grobs.arranged0[,c(1:maxMat[2])]
	}
	vp <- grid::viewport(height=grid::unit(0.9,"npc"),width=grid::unit(0.95,"npc"))
	if(!is.null(unlist(outer.text))){
		grobs.arranged <- gridExtra::arrangeGrob(grobs.arranged0,vp=vp,bottom=outer.text[[1]],left=outer.text[[2]],top=outer.text[[3]],right=outer.text[[4]])
	} else {
		grobs.arranged <- gridExtra::arrangeGrob(grobs.arranged0,vp=vp)
	}
	#if(FALSE){
	#	grobsTable.list <- list(); length(grobsTable.list) <- length(vals)
	#	for(i in 1:length(vals)){
	#		if(is.na(vals[i])){
	#			grob.i          <- grid::rectGrob(gp=grid::gpar(col=NA))
	#		} else {
	#			grob.i          <- grobs.list[[vals[i]]]
	#		}
	#		if(i==1){
	#			grobTable.i <- gridExtra::arrangeGrob(grob.i,left=m1n1.names[1],top=m1n1.names[2])
	#		}
	#		if(i %in% col1.vals){
	#			z <- which(col1.vals %in% i)
	#			grobTable.i <- gridExtra::arrangeGrob(grob.i,left=row1.names[z])
	#		}
	#		if(i %in% row1.vals){
	#			z <- which(row1.vals %in% i)
	#			grobTable.i <- gridExtra::arrangeGrob(grob.i,top=col1.names[z])
	#		}
	#		if(!(i %in% c(1,col1.vals,row1.vals))){
	#			grobTable.i <- gridExtra::arrangeGrob(grob.i)
	#		}
	#		grobsTable.list[[i]] <- grobTable.i
	#	}
	#	grobs.arranged <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=layout.mat)
	#}
	grobs.arranged
}

#' @title Arrange DAPC density plots.
#' 
#' Function to arrange plots of density for each DF or PC and each K
#' 
#' @param x A list of ggplots produced by ggscatter.dapc function
#' @param variable Either "DF" or "PC"
#' @param layout.mat Custom layout matrix to use. Default is NULL, in which case the layout will be generated automatically as a function of the number of plots.
#' @param pos.x.labs Where the default labels x axis labels (i.e., column labels) should be drawn. Default = 1 (below table of plots); alternatively 3 (above plot).
#' @param pos.y.labs Where the default labels y axis labels (i.e., row labels) should be drawn. Default = 2 (left side); alternatively 3 (right side).
#' @param row.labels.left Labels to use the left of the first column of plots. Default NULL (no labels).
#' @param col.labels.top Text labels to use above the first row of plots. Default NULL (no labels).
#' @param row.labels.right Text labels to to the right of the last column of plots. Default NULL (no labels).
#' @param col.labels.bottom Text labels to use below the bottom row of plots. Default NULL (no labels).
#' @param use.diag NULL (the default) or numerical vector with values in 1:4, indicating which sides of the gtable (1=bottom, 2=left, 3=top, 4=right) should have labels applied to the diagonal rather than the table margin plots. Ignored if layout.mat is not square. For example, if use.diag=1, then the values of col.labels.bottom will appear below the plots on the diagonal rather than below the plots of the bottom row.
#' @param pad Amount of space between plots, in units of line widths (Default 0.1).
#' @param K Which set of biplots to use. If NULL (the default), the function will attemp to draw all biplots for all K (max 25 plots).
#' @param outer.text A list with length=4, with each entry either NULL (the default) or character string to use as labels below, left, above, and to the right of the arrangement of plots.
#' @param maxMat Numerical vector with length 2 that specifies the maximum number of rows and columns of plots, respectively, in the output gtable. Default = c(7,8).
#' @return A gtable object
#' @export dapc.plot.arrange2
dapc.plot.arrange2 <- function(x,variable="DF",layout.mat=NULL,pos.x.labs=1,pos.y.labs=2,row.labels.left=NULL,col.labels.top=NULL,row.labels.right=NULL,col.labels.bottom=NULL,use.diag=NULL,pad=0.1,K=NULL,outer.text=list(NULL,NULL,NULL,NULL),maxMat=c(4,4)){
	x <- lapply(x,function(i){i[1]})
	#x3 <- do.call(c,x2)
	gg.list  <- do.call(c, x)
	if(any(lengths(gg.list))>0){
		A <- 2:(length(gg.list)+1)
		B <- which(lengths(gg.list)>0)
		Krange  <- A[B]
		gg.list <- gg.list[B]
	}
	grobs.list <- suppressWarnings(lapply(gg.list, FUN=ggplot2::ggplotGrob))
	if(is.null(layout.mat)){
		range.list   <- list(c(1),c(2:4),c(5:9),c(10:16),c(17:25))
		layout.test  <- sapply(1:length(range.list),function(x){length(grobs.list) %in% range.list[[x]]})
		if(any(layout.test)){
			layout.index  <- which(layout.test)
			layout.vector <- rep(NA,layout.index^2)
			layout.vector[1:length(grobs.list)] <- 1:length(grobs.list)
			layout.mat    <- matrix(data=layout.vector,ncol=layout.index,byrow=TRUE)
		}
		na.rowcheck <- unlist(apply(layout.mat,MARGIN=1,FUN=function(x){all(is.na(x))}))
		if(any(na.rowcheck)){
			layout.mat <- layout.mat[which(!na.rowcheck),]
		}
	}
	layout.mat <- as.matrix(layout.mat)
	kmax       <- length(layout.mat)+1
	if(pos.x.labs==1){
		col.labels.bottom  <- paste0("K=",Krange)
	}
	if(pos.x.labs==3){
		col.labels.top  <- paste0("K=",Krange)
	}
	if(pos.y.labs==2){
		row.labels.left <- paste0("K=",Krange)
	}
	if(pos.y.labs==4){
		row.labels.right <- paste0("K=",Krange)
	}
	nm   <- nrow(layout.mat)
	nn   <- ncol(layout.mat)
	len  <- nm*nn
	### character matrix of empty strings
	empty.matrix <- matrix(data="",nrow=nm,ncol=nn)
	### Index assigned to each plot of the table of plots
	index.matrix  <- matrix(1:len,ncol=nn,byrow=TRUE)
	#bottom.mat    <- matrix(data=col.labels.bottom,nrow=nm,ncol=nn)
	#top.mat   <- empty.matrix
	#right.mat <- empty.matrix
	#left.mat  <- empty.matrix
	vals <- c(layout.mat)
	grobsTable.list <- list(); length(grobsTable.list) <- length(vals)
	for(i in 1:length(vals)){
		if(is.na(vals[i])){
			grob.i          <- grid::rectGrob(gp=grid::gpar(col=NA))
			grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i)
		} else {
			grob.i          <- grobs.list[[vals[i]]]
			grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=col.labels.bottom[i])
		}
		#z <- which(index.matrix == i, arr.ind=TRUE)
		#labels.i.list  <- list(bottom.mat[z],left.mat[z],top.mat[z],right.mat[z])
		#grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=labels.i.list[[1]],left=labels.i.list[[2]],top=labels.i.list[[3]],right=labels.i.list[[4]])
		#grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=col.labels.bottom[i])
		#col.labels.bottom
	}
	grobs.arranged0 <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=index.matrix,padding=unit(pad,"line"),respect=TRUE)
	vp <- grid::viewport(height=grid::unit(0.9,"npc"),width=grid::unit(0.95,"npc"))
	if(!is.null(unlist(outer.text))){
		grobs.arranged <- gridExtra::arrangeGrob(grobs.arranged0,vp=vp,bottom=outer.text[[1]],left=outer.text[[2]],top=outer.text[[3]],right=outer.text[[4]])
	} else {
		grobs.arranged <- gridExtra::arrangeGrob(grobs.arranged0,vp=vp)
	}
	grobs.arranged
}

#' @title Arrange list of ggplot biplots
#' 
#' Function to arrange ggplot biplots for each pair of DFs or PCs of each K
#' Max number of plots = 25.
#' 
#' @param x A list of ggplots produced by ggscatter.dapc function
#' @param layout.mat Custom layout matrix to use. Default is NULL, in which case the layout will be generated automatically as a fucnton of the number of plots.
#' @param row.labels.left Text labels to use to the left of the first column of plots. Default NULL (no labels).
#' @param col.labels.top Text labels to use above the first row of plots. Default NULL (no labels).
#' @param row.labels.right Text labels to to the right of the last column of plots. Default NULL (no labels).
#' @param col.labels.bottom Text labels to use below the bottom row of plots. Default NULL (no labels).
#' @param use.diag NULL (the default) or numerical vector with values in 1:4, indicating which sides of the gtable (1=bottom, 2=left, 3=top, 4=right) should have labels applied to the diagonal rather than the table margin plots. Ignored if layout.mat is not square. For example, if use.diag=1, then the values of col.labels.bottom will appear below the plots on the diagonal rather than below the plots of the bottom row.
#' @param pad Amount of space between plots, in units of line widths (Default 0.1).
#' @param K Which set of biplots to use. If NULL (the default), the function will attemp to draw all biplots for all K (max 25 plots).
#' @param outer.text A list with length=4, with each entry either NULL (the default) or character string to use as labels below, left, above, and to the right of the arrangement of plots.
#' @param maxMat Numerical vector with length 2 that specifies the maximum number of rows and columns of plots, respectively, in the output gtable. Default = c(7,8).
#' @return A gtable object
#' @export dapc.biplot.arrange
dapc.biplot.arrange <- function(x,layout.mat=NULL,row.labels.left=NULL,col.labels.top=NULL,row.labels.right=NULL,col.labels.bottom=NULL,use.diag=NULL,pad=0.1,K=NULL,outer.text=list(NULL,NULL,NULL,NULL), maxMat=c(4,4)){
	if(length(x)>maxMat[1]){
		x <- x[1:maxMat[1]]
	}
	if(any(lengths(x)>maxMat[2])){
		x <- lapply(x,function(i){keep=min(length(i),maxMat[2]);i[1:keep]})
	}
	### Reset outer.text argument to default if it is not supplied properly, and show warning.
	if(length(outer.text)!=4 | !is(outer.text,"list")){
		outer.text <- rep(list(NULL),4)
		warning("'outer.text' ignored because not a list or length(outer.text)!=4.")
	}
	if(is.null(K)){
		gg.list  <- do.call(c, x)
	} else {
		gg.list  <- do.call(c, x[K-1])
	}
	gg.list2 <- gg.list[which(lengths(gg.list)!=0)]
	## convert a list of lists of ggplots into a list of ggplots
	grobs.list <- suppressWarnings(lapply(gg.list2, FUN=ggplot2::ggplotGrob))
	#if(length(grobs.list) > 25){
	#	grobs.list <- grobs.list[1:25]
	#}
	if(is.null(layout.mat)){
		range.list   <- list(c(1),c(2:4),c(5:9),c(10:16),c(17:25))
		layout.test  <- sapply(1:length(range.list),function(x){length(grobs.list) %in% range.list[[x]]})
		if(any(layout.test)){
			layout.index  <- which(layout.test)
			layout.vector <- rep(NA,layout.index^2)
			layout.vector[1:length(grobs.list)] <- 1:length(grobs.list)
			layout.mat    <- matrix(data=layout.vector,ncol=layout.index,byrow=TRUE)
		}
		na.rowcheck <- unlist(apply(layout.mat,MARGIN=1,FUN=function(x){all(is.na(x))}))
		if(any(na.rowcheck)){
			layout.mat <- layout.mat[which(!na.rowcheck),]
		}
	}
	layout.mat <- as.matrix(layout.mat)
	nm       <- nrow(layout.mat)
	nn       <- ncol(layout.mat)
	len      <- nm*nn
	#numplots <- length(vals)
	if(is.null(row.labels.left)){
		row.labels.left <- rep("",nm)
	}
	if(is.null(row.labels.right)){
		row.labels.right <- rep("",nm)
	}
	if(is.null(col.labels.top)){
		col.labels.top <- rep("",nn)
	}
	if(is.null(col.labels.bottom)){
		col.labels.bottom <- rep("",nn)
	}
	### Character matrix of empty strings
	empty.matrix <- matrix(data="",nrow=nm,ncol=nn)
	### Index assigned to each plot of the table of plots
	index.matrix  <- matrix(1:length(layout.mat),ncol=nn,byrow=FALSE)
	#index.matrix2 <- index.matrix[rev(1:nrow(index.matrix)),]
#	if(nm>maxMat[1]){
#		keep.rows    <- rev(rev(1:nm)[1:maxMat[1]])
#		index.matrix <- index.matrix[keep.rows,,drop=FALSE]
#		layout.mat   <- layout.mat[keep.rows,,drop=FALSE]
#		empty.matrix <- empty.matrix[keep.rows,,drop=FALSE]
#		nm <- maxMat[1]
#	}
#	if(nn>maxMat[2]){
#		keep.cols    <- 1:maxMat[2]
#		index.matrix <- index.matrix[,keep.cols,drop=FALSE]
#		layout.mat   <- layout.mat[,keep.cols,drop=FALSE]
#		empty.matrix <- empty.matrix[,keep.cols,drop=FALSE]
#		nn <- maxMat[2]
#	}
	vals     <- c(t(layout.mat))
	### Matrices holding the bottom, left, top, and right labels, respectively, for each plot.
	if(nm>1){
		if(1 %in% use.diag & nm==nn){
			bottom.mat <- empty.matrix
			bottom.mat[diag(index.matrix)] <- col.labels.bottom
		} else {
			bottom.mat <- unname(rbind(matrix(data="",nrow=(nm-1),ncol=nn),col.labels.bottom))
		}
		if(3 %in% use.diag & nm==nn){
			top.mat <- empty.matrix
			top.mat[diag(index.matrix)] <- col.labels.top
		} else {
			top.mat  <- unname(rbind(col.labels.top,matrix(data="",nrow=(nm-1),ncol=nn)))
		}
	} else {
		bottom.mat <- matrix(data=col.labels.bottom,nrow=1)
		top.mat    <- matrix(data=col.labels.top,nrow=1)
	}
	if(nn>1){
		if(2 %in% use.diag & nm==nn){
			left.mat   <- empty.matrix
			left.mat[diag(index.matrix)] <- row.labels.left
		} else {
			left.mat   <- unname(cbind(row.labels.left,matrix(data="",nrow=nm,ncol=(nn-1))))
		}
		if(4 %in% use.diag & nm==nn){
			right.mat  <- empty.matrix
			right.mat[diag(index.matrix)] <- row.labels.right
		} else {
			right.mat  <- unname(cbind(matrix(data="",nrow=nm,ncol=(nn-1)),row.labels.right))
		}
	} else {
		left.mat   <- matrix(data=row.labels.left,ncol=1)
		right.mat  <- matrix(data=row.labels.right,ncol=1)
	}
	grobsTable.list <- list(); length(grobsTable.list) <- length(vals)
	for(i in 1:length(vals)){
		#### Either creating an empty grob or getting a grob from grobs.list
		if(is.na(vals[i])){
			grob.i          <- grid::rectGrob(gp=grid::gpar(col=NA))
		} else {
			grob.i          <- grobs.list[[vals[i]]]
		}
		z <- which(index.matrix == i, arr.ind=TRUE)
		labels.i.list  <- list(bottom.mat[z],left.mat[z],top.mat[z],right.mat[z])
		if(is.null(use.diag)){
			labels.i.list2 <- list(); length(labels.i.list2) <- 4
			for(j in 1:4){
				if(labels.i.list[[j]]!=""){
					labels.i.list2[[j]] <- labels.i.list[[j]]
				}
			}
		} else {
			labels.i.list2 <- labels.i.list
		}
	#	grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=bottom.mat[z],left=left.mat[z],top=top.mat[z],right=right.mat[z])
		grobsTable.list[[i]] <- gridExtra::arrangeGrob(grob.i,bottom=labels.i.list2[[1]],left=labels.i.list2[[2]],top=labels.i.list2[[3]],right=labels.i.list2[[4]])
	}
	#return(grobsTable.list)
	#grobs.arranged <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=layout.mat,padding=unit(pad,"line"))
#	grobs.arranged0 <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=index.matrix,padding=unit(pad,"line"),respect=TRUE)
#	grobs.arranged0 <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=index.matrix[rev(1:nrow(index.matrix)),],padding=unit(pad,"line"),respect=TRUE)
	grobs.arranged0 <- gridExtra::arrangeGrob(grobs=grobsTable.list,layout_matrix=index.matrix,padding=unit(pad,"line"),respect=TRUE)
	### trim top rows if the figure has too many rows
	if(nrow(grobs.arranged0)>maxMat[1]){
		#rows.keep       <- rev(1:nrow(grobs.arranged0))[c(1:maxMat[1])]
		#grobs.arranged0 <- grobs.arranged0[rows.keep,]
		grobs.arranged0 <- grobs.arranged0[c(1:maxMat[1]),]
	}
	### trim rightmost columns if the figure has too many columns
	if(ncol(grobs.arranged0)>maxMat[2]){
		grobs.arranged0 <- grobs.arranged0[,c(1:maxMat[2])]
	} 
	vp <- grid::viewport(height=grid::unit(0.9,"npc"),width=grid::unit(0.95,"npc"))
	if(!is.null(unlist(outer.text))){
		grobs.arranged <- gridExtra::arrangeGrob(grobs.arranged0,vp=vp, bottom=outer.text[[1]],left=outer.text[[2]],top=outer.text[[3]],right=outer.text[[4]])
	} else {
		grobs.arranged <- gridExtra::arrangeGrob(grobs.arranged0,vp=vp)
	}
	grobs.arranged
}
#' da.biPlot; pca.biPlot
#' x=da.biPlot

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
	#cols10 <- c("#98ACD2","#6197ED","#7A478F","#986D65","#D9C452","#22553B","#B4993A","#0DF1FA","#E8971B","#36D8B9") # "#51BE1A" "#34A6A2" "#86D0E8" "#BDFD81" "#21594F" "#F1F9D6" "#678534" "#6B49C8" "#57191C" "#EC906B"
	cols10 <- c("#98ACD2","#E24911","#7A478F","#986D65","#D9C452","#22553B","#B4993A","#0DF1FA","#E8971B","#36D8B9") # "#51BE1A" "#34A6A2" "#86D0E8" "#BDFD81" "#21594F" "#F1F9D6" "#678534" "#6B49C8" "#57191C" "#EC906B"
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
#' @param varname Character string with name to use for the type of variable. Determined from the value of vartype when NULL (the default), such that 'varname'='Discriminant function' when vartype="df", and 'varname'="Principle component" when vartype="pc".
#' @param axis.title.cex Font cex size of x and y axis titles.
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
#' @param cellipse A positive coefficient for the inertia ellipse size. Default 1.5. Setting to zero removes ellipses. Ignored if axesell is TRUE.
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
ggscatter.dapc <- function (x, xax = 1, yax = 2, vartype="df", varname=NULL,axis.title.cex=1, grp = x$grp , cpoint=2, col = adegenet::seasun(length(levels(grp))), txt.leg = levels(grp), label = TRUE, pch = 20, solid = 0.9, hideperimeter=FALSE, show.title=TRUE,scree.da = TRUE, scree.pca = FALSE, posi.da = "bottomright", posi.pca = "bottomleft",bg="white", bg.inset = "white", ratio.da = 0.25, ratio.pca = 0.25, inset.da = 0.02, inset.pca = 0.02, inset.solid = 0.5, onedim.filled = TRUE, mstree = FALSE, lwd = 0.25, lty = 1, segcol = "black", legend = FALSE, posi.leg = "topright", cleg = 1, cstar = 1, cellipse = 1.5, axesell = FALSE, clabel = 1, xlim = NULL, ylim = NULL, grid = FALSE, addaxes = TRUE, ltyaxes=2, lwdaxes=0.5, origin = c(0,0), include.origin = TRUE, sub = "", csub = 1, possub = "bottomleft", cgrid = 1, pixmap = NULL, contour = NULL, area = NULL, label.inds = NULL, new.pred=NULL){
	if(vartype=="df"){
		ind.vals <- x$ind.coord
		if(is.null(varname)){
			varname  <- "Discriminant function"
		}
	} else {
		if(vartype=="pc"){
			ind.vals <- x$tab
			if(is.null(varname)){
				varname  <- "Principle component"
			}
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
		if(cellipse>0){
			ggscatter.tempA      <- ggplot2::ggplot(coords.df, ggplot2::aes(x=x.coords, y=y.coords,color=Cluster,shape=Cluster,fill=Cluster)) + ggplot2::scale_x_continuous(name=paste(varname,xax)) + ggplot2::scale_y_continuous(name=paste(varname,yax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = bg))  + ggplot2::stat_ellipse(color="white") + ggplot2::theme_classic() + ggplot2::geom_blank()
		} else {
			ggscatter.tempA      <- ggplot2::ggplot(coords.df, ggplot2::aes(x=x.coords, y=y.coords,color=Cluster,shape=Cluster,fill=Cluster)) + ggplot2::scale_x_continuous(name=paste(varname,xax)) + ggplot2::scale_y_continuous(name=paste(varname,yax)) + ggplot2::theme(panel.background = ggplot2::element_rect(fill = bg)) + ggplot2::theme_classic() + ggplot2::geom_blank()
		}
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
			ggscatter.tempC  <- ggscatter.tempB + ggplot2::theme(axis.title.x=ggplot2::element_text(size=(axis.title.cex*12)), axis.text.x=ggplot2::element_text(size=(axis.title.cex*10)), axis.title.y=ggplot2::element_text(size=(axis.title.cex*12)), axis.text.y=ggplot2::element_text(size=(axis.title.cex*10)))
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
			ggscatter.temp3 <- ggscatter.temp2 + suppressWarnings(ggplot2::geom_segment(data = coords.df, ggplot2::aes(x = x3, y = y3, xend = grp.center.x, yend = grp.center.y, color = Cluster),show.legend=FALSE))
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
			ggscatter.temp4 <- ggscatter.temp3 + suppressWarnings(ggplot2::geom_segment(data = coords.df, ggplot2::aes(x = tree.x0, y = tree.y0, xend = tree.x1, yend = tree.y1), color = segcol, size=lwd,linetype=lty,show.legend=FALSE))
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
		
		#xpoints <- ind.vals[grp == levels(grp)[i], pcLab]
		#ypoints <- rep(0, sum(grp == levels(grp)[i]))
		
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
#' @param min.n Integer >= 1 specifying the minimum number of non-missing alleles required to keep a site. Default = 4. If set to "all", only complete-data sites kept.
#' @param max.fMD Number in the range (0,1) specifying the maximum fraction of missing alleles at a site. Default = 1.
#' @param min.n0 Minimum number of individuals required to have at least one copy of the major allele to keep a site. Default = 2.
#' @param min.n1 Minimum number of individuals required to have at least one copy of the minor allele to keep a site. Default = 1.
#' @param which.site Character string indicating the method for choosing a site to keep for each locus (or chromosome). Default = "best", which is considered the one with the least missing data, or the first among sites tied for least missing data. Other options are "all.passing", which retains all sites (positions) that pass variation filters (min.n, min.0.n.0, min.1.n), "first" (first site kept at each locus), or "random".
#' @return List with [[1]] path to vcftools, [[2]] dataframe with input and output values for VCF filepaths, number of loci (chromosomes), sites (positions), and individuals (samples).
#' @export vcf_getSNP
vcf_getSNP      <- function(vcftools.path,vcf,out,indv.keep=NULL,which.site="best",min.n=4,max.fMD=1,min.n0=2,min.n1=1){
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
	test.sample <- unlist(gt.mat)[!is.na(unlist(gt.mat))][1]
	ploidy      <- length(unlist(strsplit(test.sample,split="[/,|]")))
	# If min.n = "all", update min.all to equal the maximum possible number of alleles at a site.
	if(min.n=="all"){
		min.n  <- ploidy*ncol(gt.mat)
	}
	n.ind  <- ploidy*ncol(gt.mat)
	min.n2 <- ceiling(n.ind-(max.fMD*n.ind))
	min.n  <- max(min.n, min.n2)
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



#' @title Filter a VCF file or vcfR object.
#' 
#' Filter or sample individuals or sites of a SNPs dataset held in either a VCF file or vcfR object. Possible filtering criteria include missingness, site variation, homozygosity (not yet implemented), and the function provides three simple methods for obtaining "unlinked" SNPs.
#' Includes options for changing some aspects of formatting, but still within withi the VCF specifications. The function depends on the vcfR. for reading and writing VCF files into R are from the package vcfR.
#' NOT QUITE FUNCTIONAL.
#' 
#' @param x Either a vcfR object or a character string to the input VCF file.
#' @param save.as Null or a character string where output VCF should be saved.
#' @param nsample A vector with two values, each NA or a number, specifying the number of individuals and sites, respectively, to randomly retain after applying the other filtering rules. Default value c(NA,NA) retains all individuals and sites that pass other filtering criteria.
#' @param specificSites Numerical vector with which sites to keep or to sample within. Default NULL.
#' @param specificIndv Character string with names of individuals to keep or to sample within. Default is NULL (all individuals kept).
#' @param only.GT Logical indicating whether or not extra, non-genotype information should be included in the genotype columns of the output VCF. Default FALSE. If TRUE, all entries of the 'FORMAT' column will be 'GT', and the genotype columns will only include genotype data.
#' @param sep.GT Character string with the separator to use for haplotypes if individuals are polyploid. Default "/", but "|" is also valid. No other separators valid.
#' @param min.n Minimum number of non-missing alleles required to keep a site. Default = 4. If set to "all", then no sites with any missing data are removed (after first filtering individuals if indv.keep is non-NULL).
#' @param min.n0 Minimum number of individuals required to have at least one copy of the major allele to keep a site. Default = 2.
#' @param min.n1 Minimum number of individuals required to have at least one copy of the minor allele to keep a site. Default = 1.
#' @param method.linked Character string indicating the method for choosing a site to keep for each locus (or chromosome). Default = "best", which is considered the one with the least missing data, or the first among sites tied for least missing data. Other options are "all.passing", which retains all sites (positions) that pass variation filters (min.n, min.0.n.0, min.1.n), "first" (first site kept at each locus), or "random".
#' @param no.indv Logical. If TRUE, specifies that the VCF does not include genotype information for specific individuals. Default FALSE. If there are no individuals and this is left FALSE then an error will be generated.
#' @return A vcfR object
#' @export filter.vcf
filter.vcf <- function(x,save.as=NULL,nsample=c(NA,NA), specificSites=NULL,specificIndv=NULL, only.GT=TRUE, sep.GT="/",min.n=4,min.n0=2,min.n1=1,method.linked="best",no.indv=FALSE){
	vcf.obj  <- vcfR::read.vcfR(x,checkFile=F)
	### Matrix of individual genotypes and potential other individual-specific information. The first column specifies the format, including if any other information is present.
	gt <- vcf.obj@gt
	### Matrix with eight columns and as many rows as sites, holding site-specific information not associated with particular individuals.
	fx <- vcf.obj@fix
	### character vector with meta data lines.
	md <- vcf.obj@meta
	if(ncol(gt)>1){
		# Data of first site of first individual. Used for detecting ploidy and the character used to separate haplotypes in the input data.
		testsite <- unname(gsub(":.+","",gt[1,2]))
		# Detect ploidy
		ploidy      <- lengths(strsplit(testsite,split="[/,|]"))
		# Determine which haplotype separator character was used for input data
		if(ploidy>1){
			testsplit1 <-  lengths(strsplit(testsite,split="[/]"))
			testsplit2 <- lengths(strsplit(testsite,split="[|]"))
			input.separator <- c("/","|")[which(c(testsplit1,testsplit2)==ploidy)]
		}
		# Genotype matrix of input file, starting with the first column with data and excluding any extra information.
		# gt0 <- gt[,-1]
		# Remove extra individual-specific info except genotype and sample names
		if(only.GT){
			gt1 <- gsub(":.+","",gt)
			fx1 <- fx
		} else {
			gt1 <- gt
			fx1 <- fx
		}
		#### Use different haplotype separator if needed.
		if(sep.GT!=input.separator){
			gt2 <- gsub(input.separator,sep.GT,gt1,fixed=TRUE)
			fx2 <- fx1
		} else {
			gt2 <- gt1
			fx2 <- fx1
		}
		### Convert nsample into a list and use NULL instead of NA
		# nsample <- lapply(nsample,FUN=function(i){if(!is.na(i))i})
		if(!is.na(nsample[1])){
			indv.keep <- sample(2:ncol(gt2),size=nsample[1])
			gt10 <- gt2[,indv.keep,drop=FALSE]
			fx10 <- fx2
		} else {
			gt10 <- gt2
			fx10 <- fx2
		}
		if(!is.na(nsample[2])){
			sites.keep <- sample(1:nrow(gt10),size=nsample[2])
			gt11 <- gt10[sites.keep,,drop=FALSE]
			fx11 <- fx10[sites.keep,,drop=FALSE]
		} else {
			gt11 <- gt10
			fx11 <- fx10
		}
	} else {
		stop("Input VCF contains no individuals. Feature not yet implemented.")
		if(!no.indv){
			message("Input VCF contains no individuals")
		}
	}
	vcf.out <- vcf.obj
	vcf.out@fix <- fx11
	vcf.out@gt  <- gt11

	if(!is.null(save.as)){
		vcfR::write.vcf(x=vcf.out,file=save.as)
	}
	vcf.out
}


#' Simulate SNPs as VCF
#' 
#' Create a simulated SNP dataset from either an input VCF file, an object of class vcfR, or by specifying characteristics of the simulated dataset using the function arguments.
#' Most of the work is done by a call to the function 'glSim' from the 'adegenet' package.
#' Note: At present I do not recommend setting LD=TRUE (or LD=NULL if the input VCF includes linked SNPs), because CHROM and POS columns of the output VCF will not reflect the presence of linkage. Check back for updates.
#' 
#' @param x Path to input VCF or a vcfR object. Default NULL. Used to extract probabilities that a site has a particular pair of Reference and Alternative nucleotides. If NULL 'probs' must be supplied.
#' @param save.as Null or a character string where the output VCF containing the simulated dataset should be saved.
#' @param RA.probs Either NULL (the default), one of the character strings "equal" or "empirical", or a numerical vector with probability of each combination of reference and alternative allele. Numerical vectors must have length 12 and sum to 1; if elements of the vector are unnamed, assumed names will be c("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"), with the first character specifying the reference allele and the second the alternative allele.
#' If RA.probs is NULL it is coerced to 'empirical' if 'x' is non-NULL, or coerced to 'equal' if 'x' is NULL.
#' Setting RA.probs to 'equal' is equivalent to the vector produced by sapply(c("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG"),assign,(1/12)).
#' If RA.probs="empirical", the frequencies of each Reference:Alternative allele combination is calculated for the input dataset and used for sampling probabilities for the simulated dataset.
#' @param n.ind Number of individuals to include in simulated dataset. Default NULL, in which case the number of simulated individuals will match the number of individuals in 'x'.
#' @param n.snps Number of SNPs. If n.snp.nonstruc is NULL and x is non-NULL, n.snp.nonstruc will equal half the number of SNPs in the input dataset.
#' @param snp.str Fraction of SNPs that are structured [WITHIN?] populations. Default 0. Increasing more than 0 seems generate at least twice as many simulated populations as specified with K.
##' @param n.snp.nonstruc Number of nonstructured SNPs. If n.snp.nonstruc is NULL and x is non-NULL, n.snp.nonstruc will equal half the number of SNPs in the input dataset.
##' @param n.snp.struc Number of structured SNPs. If n.snp.nonstruc is NULL and x is non-NULL, 'n.snp.nonstruc' will equal the number of SNPs in the input dataset minus the value of 'n.snp.nonstruc'. Meaningless if K = 1.
#' @param ploidy Number indicating ploidy of individuals. Probably can only be 1 or 2?
#' @param K Number of populations (>=2). Default 2.
#' @param include.missing Logical indicating whether or not missing data should be included in the simulated dataset.
#' @param fMD NULL or a number in the range (0,1] indicating the fraction of genotypes that should be missing in the simulated dataset. If NULL (the default) and 'include.missing'=TRUE and 'x' non-NULL, the fraction of genotypes missing in the simulated dataset will be approximately equal to the fraction of missing genotypes in the input dataset.
#' @param pMDi NULL (the default) or a numerical vector with the proportion of missing data that should be attributed to each simulated individual. If NULL, ignored unless 'x' is non-NULL and 'include.missing'=TRUE, in which case pMDi values will be sampled from missing data proportions of the input VCF, producing similar missing data structure for input and simulated datasets.
#' @param LD NULL or a logical (TRUE or FALSE) indicating whether or not snps should be simulated under linkage disequilibrium. Default FALSE. If LD='NULL' (will be the default in future versions) and 'x' is non-NULL, LD is coerced to FALSE if all sites (positions) of the input dataset are unlinked (on different blocks/chromosomes), otherwise LD is coerced to TRUE.
#' @param block.minsize Number with minimum. Ignored if 'LD' is FALSE or coerced to FALSE. Default 10.
#' @param block.maxsize Number indicating the maximum size of linkage blocks. Ignored if 'LD' is FALSE or coerced to FALSE, or if 'use.maxLDx' is TRUE and 'x' is non-NULL, in which case the max block size is equal to the max position of any snps in a linkage block). Default 1000.
#' @param use.maxLDx Logical indicating whether or not the maximum linkage block size should be set as the max position of any snp within a linkage block of the input data. Default FALSE, which means that max linkage block size should be the value of 'block.maxsize'.
#' @param sim.coords Logical indicating if a simulated coordinates (sample localities) should be produced for each simulated individual. Default FALSE.
#' @param regionsize Number between 0 and 1 specifying the fraction of the possible sampling area (the window region set by 'wnd' argument) in which points may be distributed. Default 0.25. When 'wnd' is default (entire Earth), regionsize is default (0.25), and ('over.land' condition is FALSE), then minimum convex hull of all sampled points is expected to cover 1/4 of Earth.
#' @param wnd Either a character string or vector describing one or more regions of Earth, or a length four numerical vector specifying the bounding box (longitude and latitude ranges) for the region where sampling is allowed. Default is to include all of Earth c(-180,180,-90,90). 
#' @param interactions Numeric vector with the minimum and maximum amount of overlap between each pair of groups (aka populations), calculated as (intersect area)/(minimum of non-intersected area for each group). Default c(0,1) allows for all possible scenarios. Examples: c(0,0) specifies that groups must be allopatric; c(1,1) requires complete overlap of groups, which is not realistic given the stochasticity determining region sizes; c(0.5,1) requires that at least half-overlaps between groups; c(0.2,0.25) specifies a small contact zone.
#' @param over.land, and 'interactions'
#' @param ... Additional arguments passed to 'adegenet::glSim' to control SNPs simulation, or to 'misc.wrappers::rcoords' to further control the simulation of geographic localities when 'sim.coords' is TRUE. Possible 'glSim' arguments include 'grp.size', 'pop.freq', 'alpha', 'parallel', and 'theta' (see ?adegenet::glSim for details).
#' @return An object with class vcfR (see 'vcfR' package for details regarding this class)
#' @export sim.vcf
sim.vcf <- function(x=NULL, save.as=NULL, RA.probs=NULL, n.ind=NULL, n.snps=NULL, snp.str = 0, ploidy=NULL, K=2, include.missing=FALSE, fMD=NULL, pMDi=NULL, LD=FALSE, block.minsize=10, block.maxsize=NULL, use.maxLDx=FALSE, sim.coords=FALSE , regionsize=0.7, wnd=c(-180,180,-90,90), interactions=c(0,1), over.land=TRUE, ...){
#	x=NULL; save.as=NULL; RA.probs=NULL; n.ind=NULL; n.snps=NULL; snp.str = 0; ploidy=NULL; K=2; include.missing=FALSE; fMD=NULL; LD=NULL; block.minsize=10; block.maxsize=NULL; use.maxLDx=FALSE
	# additional.args  <- list(...)
	RA.pairnames <- c("AC","AG","AT","CA","CG","CT","GA","GC","GT","TA","TC","TG")
	if(!is.null(x)){
		vcf   <- vcfR::read.vcfR(x,verbose=F,checkFile=F,convertNA=FALSE)
		# genotype matrix
		gt    <- vcf@gt
		### genotype matrix without format column or extra (non-genotype) information.
		gt.simple <- gsub(":.+","",gt[,-c(1)])
		# fixed matrix
		fx    <- vcf@fix
		# metadata strings
		met   <- vcf@meta
		# character strings with reference and alternative allele for each site
		RA    <- paste0(fx[,"REF"],fx[,"ALT"])
		# number of each combination of reference and alternative alleles
		fRA   <- table(RA)
		# Proportion of sites with each combination of reference and alternative alleles
		pRA   <- fRA/sum(fRA) ### Use this for drawing samples from names(RA)
		# chromosome (linkage block) ID at each site
		chr   <- fx[,"CHROM"]
		# unique set of chromosome (linkage block) IDs
		uchr  <- unique(chr)
		# position of each variant within its linkage block
		pos   <- fx[,"POS"]
		# maximum position value of any variant within any linkage block.
		maxLDx <- max(pos)
		if(is.null(n.ind)){
			n.ind <- (ncol(gt)-1)
		}
		if(is.null(n.snps)){
			n.snps <- nrow(gt)
		}
		if(is.null(ploidy)){
			testsite <- unname(gsub(":.+","",gt[1,2]))
			ploidy   <- lengths(strsplit(testsite,split="[/,|]"))
		}
		if(ploidy>1){
			testsplit1 <-  lengths(strsplit(testsite,split="[/]"))
			testsplit2 <- lengths(strsplit(testsite,split="[|]"))
			input.separator <- c("/","|")[which(c(testsplit1,testsplit2)==ploidy)]
			# String indicating missing genotypes in input dataset.
			xmval <- paste(rep(".",ploidy),collapse=input.separator)
		} else {
			xmval <- "."
		}
		# The total number of missing genotypes across all individuals and variants
		tnMD <- length(grep(xmval,gt.simple,fixed=TRUE))
		# If any missing data exists, calculate number and probability of missing genotypes for each individual
		if(tnMD>0){
			indv.nMD <- apply(X=gt.simple,MARGIN=2,FUN=function(xn){length(grep(xmval,xn,fixed=TRUE))})
			# The probability that a genotype selected randomly from the input VCF belongs to a particular individual conditioned on the genotype being missing.
			indv.pMD <- indv.nMD/tnMD
			# x variable of density distribution of indv.pMD; may be used to generate sim.indv.pMD for assigning missing data to simulated individuals.
			indv.dMD <- density(indv.pMD)$x
		} else {
			indv.pMD <- rep(0,n.ind)
		}
		if(tnMD==0 & is.null(fMD)){
			include.missing=FALSE
		}
		# whether or not any snps in the input dataset are on different linkage blocks.
		if(length(chr) == length(uchr)){
			LDx <- FALSE
			if(is.null(LD)){
				LD <- FALSE
			}
		} else {
			LDx <- TRUE
			if(is.null(LD)){
				LD <- TRUE
			}
			if(use.maxLDx){
				block.maxsize <- maxLDx
			}
		}
	} else {
		if(length(c(n.ind,n.snp.nonstruc,ploidy)) !=3){
			stop("'n.ind','n.snp.nonstruc', and 'ploidy' must be non-NULL if 'x' is NULL")
		}
	}
	if(is.null(RA.probs)){
		if(is.null(x)){
			RA.probs <- "equal"
		} else {
			RA.probs <- "empirical"
		}
	}
	if(is(RA.probs[1],"character")){
		if(RA.probs=="equal"){
			RA.probs <- sapply(RA.pairnames,assign,(1/12))
		} else {
			if(RA.probs=="empirical"){
				if(is.null(x)){
					stop("'x' must be non-NULL if 'RA.probs'='empirical' ")
				} else {
					RA.probs <- pRA
				}
			}
		}
	} else {
		if(is(RA.probs,"numeric") & length(RA.probs)==12 & sum(RA.probs)==1){
			if(is.null(names(RA.probs))){
				names(RA.probs) <- RA.pairnames
			} else {
				if(!all(RA.pairnames %in% names(RA.probs))){
					stop("invalid RA.probs")
				}
			}
		} else{
			stop("invalid RA.probs")
		}
	}
	
	n.snp.struc    <- round(n.snps*snp.str)
	n.snp.nonstruc <- n.snps-n.snp.struc
	### Perform simulation
#	simK2 <- adegenet::glSim(n.ind=n.ind,n.snp.nonstruc=n.snp.nonstruc,ploidy=ploidy,K=K,args)
#	sim <- adegenet::glSim(n.ind=n.ind, n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, k=K, LD=LD, block.minsize = block.minsize, block.maxsize = block.maxsize, sort.pop=TRUE,...)
	# glSim      <- adegenet::glSim
	#glSim.args <- list(n.ind=n.ind,n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, LD=LD, block.minsize = block.minsize) # , block.maxsize = block.maxsize)
	#names(...) %in% names(supplied)
	# sim   <- evalfun(fun.name="glSim", n.ind=n.ind, n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, k=K, LD=LD, block.minsize = block.minsize, block.maxsize = block.maxsize, sort.pop=TRUE)
	#sim   <- evalfun(fun.name="glSim", supplied.list=glSim.args, k=K, sort.pop=TRUE)
#	sim0   <- evalfun(fun.name="glSim", k=K, sort.pop=TRUE,n.ind=(n.ind*K),n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, LD=LD)
	K.weights <- rep(1, K)
	K.sizes   <- table(sample(1:K,size=n.ind,replace=T,prob=K.weights))
	while(length(K.sizes) < K){
		K.sizes   <- table(sample(1:K,size=n.ind,replace=T,prob=K.weights))
	}
	PASS <- FALSE
	while(!PASS){
		### The glSim function requires each simulated population to have at least 10 individuals.
		sim0   <- adegenet::glSim(k=K, sort.pop=TRUE, n.ind=(n.ind*K*10), n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, LD=LD)
		indvs.pops <- c(sim0@other$ancestral.pops)
		if(length(table(indvs.pops)) != K){
			next
		}
		if(any(table(indvs.pops) < K.sizes)){
			next
		}
		ind.keep  <- unlist(lapply(X=1:K,FUN=function(x) grep(x,indvs.pops)[1:K.sizes[x]]))
		sim       <- sim0[ind.keep]
		PASS <- TRUE
	}
	#sim   <- adegenet::glSim(k=K, sort.pop=TRUE, n.ind=(n.ind*K), n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, LD=LD)
	#return(sim)
	# sim5   <- adegenet::glSim(k=K, sort.pop=TRUE, n.ind=(n.ind*K), n.snp.nonstruc=n.snp.nonstruc, n.snp.struc=n.snp.struc, ploidy=ploidy, LD=LD)
	# pops.index.mat <- matrix(data=1:nrow(sim0), ncol=K)
	# ind.keep  <- sort(unlist(lapply(1:K,FUN=function(x){sample(pops.index.mat[,x],size=K.sizes[x],replace=F)})))
	### Ancestral population assignments of simulated samples
	anc.pops <- gsub(".", "", as.character(sim@other$ancestral.pops), fixed=T)
	### Genotype matrix of simulated dataset
	sim.gt0 <- t(as.matrix(sim))
	# Check that each site is actually a SNP. Rare for invariant sites to be generated by glSim, but it occurs at a very low frequency.
	sitediv0 <- apply(X=sim.gt0, MARGIN=1, FUN=function(x){length(unique(x))>1})
	# Randomly set each invariant site equal to some other site.
	while(!all(sitediv0)){
		sites.resample <- which(!sitediv0)
		sim.gt0[sites.resample,] <- sim.gt0[sample(1:n.snps,size=length(sites.resample)),]
		sitediv0 <- apply(X=sim.gt0, MARGIN=1, FUN=function(x){length(unique(x))>1})
	}
	### Introduce missing data
	if(include.missing){
		if(is.null(fMD)){
			### Calculate fraction of missing data at each site in the input dataset.
			sites.md     <- apply(X=gt.simple,MARGIN=1,FUN=function(j){length(grep(pattern="[.]",x=j))/length(j)})
			### Sample values from sites.md to use for sites in the simulated dataset.
			sim.sites.fFM <- sample(sites.md,size=n.snps,replace=FALSE)
			### For each simulated site, the number of individuals that should have genotype missing
			sim.sites.nMD <- floor(sim.sites.fFM*n.ind)
			# String to use for missing genotypes. Will temporarily use NA for missing data and then convert NAs to this string.
			mval <- paste(rep(".",ploidy),collapse="/")
			#sim.gt4 <- mapply(A=sim.gt3,B=sim.sites.nMD,MARGIN=1,FUN=function(a,b){res=a; res[sample(length(a),size=b)] <- mval})
			#sim.gt3, B=sim.sites.nMD
###			sim.gt4 <- lapply(X=1:length(sim.sites.nMD),FUN=function(x){res=sim.gt3[x,]; res[sample(n.ind,size=sim.sites.nMD[x])] <- NA;res})
			### Distribute missing data across individuals in a similar way as in the input dataset
			if(is.null(pMDi)){
				if(!is.null(indv.pMD)){
					### Probability weights
					pMDiw <- sample(indv.pMD,size=n.ind,replace=T)
					### Probabilities
					pMDi  <- unname(c(pMDiw/sum(pMDiw)))
					
				} else {
					pMDi <- (1/n.ind)
				}
			}
			### Introduce missing data at each site.
			sim.gt00 <- do.call(rbind,lapply(X=1:length(sim.sites.nMD),FUN=function(x){res=sim.gt0[x,]; res[sample(n.ind,size=sim.sites.nMD[x],replace=FALSE,prob=pMDi)] <- NA; res}))
			# check which sites have less than two genotypes
###			sitediv <- vapply(X=sim.gt4,FUN=function(x){length(table(x))},FUN.VALUE=1)
			sitediv <- apply(X=sim.gt00,MARGIN=1,FUN=function(x){length(unique(x))>1})
			while(!all(sitediv)){
				sites.resample <- which(!sitediv)
				sim.gt00[sites.resample,] <- sim.gt00[sample(1:n.snps,size=length(sites.resample)),]
				sitediv <- apply(X=sim.gt00,MARGIN=1,FUN=function(x){length(unique(x))>1})
			}
		}
		if(any(is.na(sim.gt00))){
			sim.gt00[is.na(sim.gt00)] <- mval
		}
	} else {
		### Sample the fraction of missing data for each site from the mean distribution.
		sim.gt00 <- sim.gt0
	}
	##
	mode(sim.gt00) <- "character"
	sim.gt1 <- gsub("^0$","0/0",sim.gt00)
	sim.gt2 <- gsub("^1$","0/1",sim.gt1)
	sim.gt3 <- gsub("^2$","1/1",sim.gt2)
	# Names to use for simulated individuals
	samplenames <- paste0("Sample",c(1:ncol(sim.gt3)),".",anc.pops)
	# Add sample names
	colnames(sim.gt3) <- samplenames
	### Change sim.gt3 to the penultimate gt version
	sim.gt <- cbind("FORMAT"=rep("GT",nrow(sim.gt3)),sim.gt3)
	### Arrange fixed matrix for simulated dataset
	# Sample Reference and Alternative allele pair for each SNP
	sim.RA <- do.call(rbind,strsplit(sample(names(RA.probs),size=n.snps,replace=TRUE,prob=RA.probs),split=""))
	colnames(sim.RA) <- c("REF","ALT")
	sim.fx <- cbind("CHROM"=c(1:n.snps),"POS"=rep(1,n.snps),"ID"=rep(NA,n.snps),sim.RA,"QUAL"=rep(NA,n.snps),"FILTER"=rep("PASS",n.snps),"INFO"=rep(NA,n.snps))
	### Construct metadata lines.
	sim.meta <- c("##fileformat=VCFv4.2",
		paste0("##fileDate=",gsub("-","",Sys.Date())),
		"##source=\"misc.wrappers\"",
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
	
	
	#return(list(n.per.pop=anc.pops,K=K))
	if(sim.coords){
		#n.per.pop    <- table(as.numeric(gsub("pop","",anc.pops)))
		n.per.pop    <- table(anc.pops)
		#simcoords <- rcoords(regionsize=0.7,samplesize=n.per.pop,n.grp=K,show.plot=F,wnd=c(-180,180,-60,60),interactions=c(0,0),return.as="DF_plot")
		simcoords <- evalfun("rcoords", regionsize=regionsize, samplesize=n.per.pop, n.grp=K, wnd=wnd, over.land=over.land, interactions=interactions, show.plot=TRUE)
		#return(simcoords)
		if(is(simcoords,"list")){
			simcoords.df <- simcoords$coords.df
		} else {
			if(is(simcoords,"data.frame")){
				simcoords.df <- simcoords
			}
		}
		rownames(simcoords.df) <- samplenames
		colnames(simcoords.df) <- c("longitude","latitude","group")
		#save.as3 <- paste0(tools::file_path_sans_ext(save.as),"_simulated_coords.txt")
		save.as3 <- file.path(dirname(save.as),paste0(gsub("[.vcf].+","",basename(save.as)),"_simulated_coords.txt"))
		write.table(x=simcoords.df,file=save.as3, row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
		if(is(simcoords,"list")){
			if("map.ggplot" %in% names(simcoords)){
				gg.map <- simcoords$map.ggplot
				save.as4 <- file.path(dirname(save.as),paste0(gsub("[.vcf].+","",basename(save.as)),"_simulated_coords_map.pdf"))
				pdf(file=save.as4,width=8,height=6,onefile=TRUE)
					grid::grid.newpage()
					grid::grid.draw(gg.map)
				dev.off()
			}
		}
		# Add coordinates as metadata lines?
	}

	#### Create vcfR object with simulated dataset
	vcf.sim <- new("vcfR", meta=sim.meta,fix=sim.fx,gt=sim.gt)
	if(!is.null(save.as)){
		### Use extension ".vcf.gz" even is ".vcf" is specified as the extension.
		save.as2 <- gsub(".vcf$",".vcf.gz",save.as)
		vcfR::write.vcf(x=vcf.sim,file=save.as2)
		vcf.sim <- vcfR::read.vcfR(save.as2)
	}

	### Return the simulated dataset as vcfR object
	# If save.as is null, then the object returned will always report "zero missing data", although missing data may exists. Writing and rereading the VCF removes this vcfR bug. May need to use NA in gt matrix of vcfR.
	if(sim.coords){
		if(exists("bothmaps")){
			return(list(simulated.vcf=vcf.sim,simulated.LonLat=simcoords.df,gg.map=bothmaps))
		} else {
			return(list(simulated.vcf=vcf.sim,simulated.LonLat=simcoords.df))
		}
	} else {
		return(vcf.sim)
	}
}
#' @examples
#' library(misc.wrappers)
#' # Define path to example input VCF containing 5000 variants and 100 individuals.
#' vcf.path <- file.path(system.file("extdata", package = "misc.wrappers"),"example.vcf.gz")
#' 
#' # Simulate a dataset of 1000 variants and 50 individuals in one population, and save the simulated dataset in the current directory as "example_simulated_K2.vcf.gz"
#' simK2 <- sim.vcf(x=vcf.path,save.as="example_simulated_K2.vcf.gz",n.ind=100,n.snps=1000,K=2)
#' 
#' # Simulate a dataset of 1000 variants and 50 individuals in one population, and save the simulated dataset in the current directory as "example_simulated_K3.vcf.gz"
#' simK3<- sim.vcf(x=vcf.path,save.as="example_simulated_K3.vcf.gz",n.ind=100,n.snps=1000,K=3)
#' 
#' # Simulate a dataset of 1000 variants and 50 individuals in one population, and save the simulated dataset in the current directory as "example_simulated_K4.vcf.gz"
#' simK4 <- sim.vcf(x=vcf.path,save.as="example_simulated_K4.vcf.gz",n.ind=100,n.snps=1000,K=4)

#' @title Components plottable
#' 
#' Creates a list with two logical data frames, the first with K vs. DF, and the second with K vs. PC, with TRUE values indicating that it is possible to make a density plot for K=rowname(m) vs. DF=colname(n).
#' 
#' @param x list of DAPC objects
#' @return Returns a list with two data frames, each with rownames indicating K for which density can be plotted (because smallest cluster size has > 1 individual). Columns show names of discriminat functions (first data frame) or principle components (second data frame). Values are TRUE or FALSE.
#' The function is useful for creating layout matrices for plotting all possible density plots.
#' @export plottable.dapc
plottable.dapc <- function(x){
	grps.mat           <- do.call(cbind,lapply(x,FUN=function(j){c(j$assign)}))
	colnames(grps.mat) <- paste0("K", (1:length(x))+1)
	grp.minsizes       <- apply(grps.mat, MARGIN=2,FUN=function(x){min(table(x))})
	if(any((grp.minsizes > 1))){
		Ks.plottable <- unname(which((grp.minsizes > 1))+1)
	} else {
		return(NULL)
	}
	result <- list(); length(result) <- 2
	for(i in 1:2){
		if(i==1){
			Ks.n <- sapply(x,function(x){x$n.da})
			stat <- "DF"
		}
		if(i==2){
			Ks.n <- sapply(x,function(x){x$n.pc})
			stat="PC"
		}
		Ks.df2plot              <- lapply(X=Ks.plottable,FUN=function(x){A=(1:Ks.n[(x-1)]); names(A)=paste0(stat,(1:Ks.n[(x-1)])); as.data.frame(t(A))})
		Ks.df2plot.df           <- as.data.frame(data.table::rbindlist(Ks.df2plot,fill=TRUE))
#		rownames(Ks.df2plot.df) <- paste0("K",Ks.plottable)
		rownames(Ks.df2plot.df) <- Ks.plottable
		result[[i]] <- !is.na(Ks.df2plot.df)
	}
	names(result) <- c("DF","PC")
	result
}


