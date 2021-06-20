#' @title plot_eems function
#' 
#' Invokes the reemsplots2 function make_eems_plots, with additional functionality such as adding island borders and points showing individual coordinates..
#' 
#' @param xdir Character string to directory of a particular mcmc chain, or, to a directory containing subdirectories "/mcmc/chain*".
#' @param save.in Path to directory where output files should be saved. Default is NULL, in which case the value of 'mcmcdir' is used.
#' @param plot.coords Whether or not to overlay the sample coordinates on EEMS maps. Default TRUE.
#' @param plot.geography Whether or not to overlay country borders on EEMS maps. Default TRUE.
#' @param mask.oceans Whether or not mask oceans on EEMS maps. Default TRUE.
#' @param include.out Character string vector indicating the type of output files to generate.
#' @param usechains Integer vector specifying which chains to make plots for; Default is NULL, in which case plots are made for all chains.
#' @return List of ggplots
#' @export plot_eems
plot_eems <- function(xdir, save.in=NULL, plot.coords=T, plot.geography=T, mask.oceans=T, include.out=c("pdf","raster"), usechains=NULL){
	#if(is.null(save.in)){
	#	save.in <- mcmcdir
	#}
	if("mcmc" %in% list.files(xdir)){
		nchains <- length(list.files(xdir,pattern="^runeems_snps_chain.*.sh$"))
		mcmcdirs <- file.path(xdir,"mcmc",paste0("chain",1:nchains))
	} else {
		nchains <- 1
		mcmcdirs <- xdir
		### overrides value of argument
		plot.coords <- FALSE
	}
	spdf_oceans_10 <- misc.wrappers::spdf_oceans_10
	oceans_sf      <- sf::st_as_sf(spdf_oceans_10)
	world_sf       <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
	if(is.null(usechains)){
		usechains <- 1:nchains
	} else {
		usechains <- intersect(usechains,1:nchains)
	}
	result <- list(); length(result) <- nchains
	if("raster" %in% include.out & length(usechains) > 1){
		bricklist <- list(); length(bricklist) <- nchains
	}
	for(i in usechains){
		mcmcdir <- mcmcdirs[[i]]
		if(is.null(save.in)){
			save.in <- mcmcdir
		}
		### Generate EEMS plots
		######
		## CHANGE the NEXT LINE TO FIRST SOURCE FUNCTIONS FROM 'make_eems_plots_JLW.R' (in 'misc.wrappers/inst/extdata'), then use the version of the make_eems_plots function from that file rather than the version in reemsplots2
		source(file.path(system.file("extdata", package = "misc.wrappers"),"make_eems_plots_JLW.R"))
		# plots.chain     <- suppressWarnings(reemsplots2::make_eems_plots(mcmcdir, longlat = TRUE))
		plots.chain     <- suppressWarnings(make_eems_plots(mcmcdir, longlat = TRUE))
		if(plot.coords){
			#coordinates <- read.table(paste0(mcmcdir,"/rdistoDemes.txt"),header=F)[,c(1:2)]
			#colnames(coordinates) <- c("Lon","Lat")
			coordinates <- read.table(file.path(xdir,"data","data.coord"),header=F)[,c(1:2)]
			colnames(coordinates) <- c("Lon","Lat")
		}
		mapplot  <- list(); length(mapplot) <- 4
		for(i in 1:4){
			mapplot.i  <- plots.chain[[i]]
			coords <- unique(mapplot.i$data[,c(1:2)])
			colnames(coords) <- c("Lon","Lat")
			x.min    <- min((coords[,1]-0.5))
			x.max    <- max((coords[,1]+0.5))
			y.min    <- min((coords[,2]-0.5))
			y.max    <- max((coords[,2]+0.5))
			# Add country borders
			if(plot.geography){
				mapplot2.i <- suppressMessages(mapplot.i + ggplot2::geom_sf(data=world_sf, colour = "black", fill = NA, inherit.aes=FALSE)) # + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
			} else {
				mapplot2.i <- mapplot.i
			}
			# Mask oceans
			if(mask.oceans){
				mapplot3.i <- mapplot2.i + ggplot2::geom_sf(data=oceans_sf, colour = NA, fill = "white", inherit.aes=FALSE) # + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
			} else {
				mapplot3.i <- mapplot2.i
			}
			# Crop to window defined by xmin, xmax, ymin, ymax
			mapplot4.i <- mapplot3.i + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max), expand=FALSE)
			# Add points for coordinates of samples
			if(plot.coords){
				mapplot[[i]] <- mapplot4.i + ggplot2::geom_point(data = coordinates, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			} else {
				mapplot[[i]] <- mapplot4.i
			}
		}
		result.gglist <- c(mapplot,plots.chain[5:8])
		if("pdf" %in% include.out){
			save.as.pdf <- file.path(save.in,"EEMS_maps.pdf")
			pdf(height=6, width=8, file=save.as.pdf, onefile=TRUE)
				print(result.gglist)
			dev.off()
		}
		if("raster" %in% include.out){
			save.as.raster  <- file.path(save.in,"EEMS_maps.tif")
			plots.dat       <- lapply(X=plots.chain[1:4],FUN=function(x) as.data.frame(x$data))
			rasterlist <- list(); length(rasterlist) <- 4
			for(j in 1:4){
				dat.temp <- plots.dat[[j]]
				e.temp   <- raster::extent(dat.temp[,1:2])
				r.temp   <- raster::raster(e.temp, ncol=length(unique(dat.temp$x)), nrow=length(unique(dat.temp$y)))
				rasterlist[[j]] <- raster::rasterize(dat.temp[, 1:2], r.temp, dat.temp[,3], fun=mean, background=-9999)
			}
			# Hold the rasters in a rasterbrick
			res.brick <- raster::brick(rasterlist)
			# Write the brick as a multi-band geotiff
			wb <- suppressWarnings(raster::writeRaster(x=res.brick, filename=save.as.raster, format="GTiff"))
			if(length(usechains) > 1){
				bricklist[[i]] <- res.brick
			}
		}
		result[[i]] <- result.gglist
	}
	### For each of the four EEMs plots, generate a raster with chain means. Hold these in a brick and write as a GeoTiff.
	if("raster" %in% include.out & length(usechains) > 1){
		bricklist2 <- bricklist[usechains]
		bp1 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[1]]}))
		bp2 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[2]]}))
		bp3 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[3]]}))
		bp4 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[4]]}))
		rp1.mean <- raster::mean(bp1)
		rp2.mean <- raster::mean(bp2)
		rp3.mean <- raster::mean(bp3)
		rp4.mean <- raster::mean(bp4)
		bmeans   <- raster::brick(list(rp1.mean,rp2.mean,rp3.mean,rp4.mean))
		bw.mean  <- suppressWarnings(raster::writeRaster(x=bmeans,filename=file.path(xdir,"EEMS_maps_chainMeans.tif",format="GTiff")))
	}
	result2 <- result[usechains]
	result2
}




