#' @title plot_eems function
#' 
#' Invokes the reemsplots2 function make_eems_plots, with additional functionality such as adding island borders and points showing individual coordinates..
#' 
#' @param mcmcdir Character string to directory of a particular mcmc chain
#' @param coords Optional character string to sample coordinates file.
#' @param output.plot.path Character string with path where to save output. If NULL, the plot is printed but not saved.
#' @return NULL; plots are printed to screen and saved to output.plot.path
#' @export plot_eems
plot_eems <- function(mcmcdir,coords=NULL,output.plot.path=NULL){
	spdf_world_10    <- rnaturalearth::ne_countries(scale=10)
	spdf_oceans_10   <- invisible(ne_download(scale = 10, type = 'ocean', category = 'physical'))
	### Generate EEMS plots
	plots.chain     <- suppressWarnings(reemsplots2::make_eems_plots(mcmcdir, longlat = TRUE))
	### extract the first eems plot (this is a ggplot)
#	eems1.plot  <- plots.chain[[1]]
	### instructions to make ggplot of first eems plot
#	eems1.build <- ggplot2::ggplot_build(eems1.plot)
	### A list of instructions for each ggplot in plots.chain
	eems.build <- lapply(plots.chain,FUN=ggplot2::ggplot_build)
	###
	eems.data <- lapply(X=eems.build,FUN=function(x){x[[1]][[1]]})
	### The data frame used to make the ggplot of the first eems plot
#	eems1.data      <- eems1.build[[1]][[1]] # a 41590x15 data frame; each entry of column one indicates the fill color of a tile as a hex string
#	eems1.df.data   <- eems1.data[,setdiff(colnames(eems1.data),c("x","y","colour","alpha","width","height"))]
#	eems1.data.spdf <- sp::SpatialPointsDataFrame(coords=eems1.data[,c("x","y")],data=eems1.df.data)
	###
	eems.df.data   <- lapply(X=eems.data[1:4],FUN=function(df){df[,c("fill","group")]})
#	eems1.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[1]][,c("x","y")],data=eems.df.data[[1]])
#	eems2.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[2]][,c("x","y")],data=eems.df.data[[2]])
#	eems3.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[3]][,c("x","y")],data=eems.df.data[[3]])
#	eems4.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[4]][,c("x","y")],data=eems.df.data[[4]])
	eems.data.spdf <- lapply(X=1:4,FUN=function(i){sp::SpatialPointsDataFrame(eems.data[[i]][,c("x","y")],data=eems.df.data[[i]],proj4string=CRS(crs.string))})
	### Setting the crs
	crs.string <- "+init=EPSG:4326"
#	lapply(X=1:4,FUN=function(x){raster::crs(eems.data.spdf[[x]])   <- crs.string})
	# suppressWarnings(raster::crs(eems1.data.spdf)   <- crs.string)
	suppressWarnings(raster::crs(spdf_world_10)     <- crs.string)
	suppressWarnings(raster::crs(spdf_oceans_10)    <- crs.string)
	intersect.eems              <- lapply(X=1:4,function(x){sp::over(eems.data.spdf[[x]],spdf_world_10)})
	### interections of each plot with land
	eems1.intersect.spdf <- intersect.eems[[1]]
	eems2.intersect.spdf <- intersect.eems[[2]]
	eems3.intersect.spdf <- intersect.eems[[3]]
	eems4.intersect.spdf <- intersect.eems[[4]]
	### unique rows of each plot
	eems1.intersect.spdf.unique <- eems1.intersect.spdf[apply(X=eems1.intersect.spdf,MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	eems2.intersect.spdf.unique <- eems2.intersect.spdf[apply(X=eems2.intersect.spdf,MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	eems3.intersect.spdf.unique <- eems3.intersect.spdf[apply(X=eems3.intersect.spdf,MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	eems4.intersect.spdf.unique <- eems4.intersect.spdf[apply(X=eems4.intersect.spdf,MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	### All four data frames held in a list
	eems.intersect.spdf.unique <- list(eems1.intersect.spdf.unique,eems2.intersect.spdf.unique,eems3.intersect.spdf.unique,eems4.intersect.spdf.unique)
	#intersect.eems.unique       <- eems1.data.spdf[apply(X=intersect.eems,MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	### Each plot
	for(i in 1:4){
		intersect.eems.unique.plot <-  sp::plot(eems.intersect.spdf.unique[[i]],col=eems.intersect.spdf.unique[[i]]$fill)
		# intersect.eems1.unique.plot  <- sp::plot(eems1.intersect.spdf.unique,col=eems1.intersect.spdf.unique$fill)
		# intersect.eems2.unique.plot  <- sp::plot(eems2.intersect.spdf.unique,col=eems2.intersect.spdf.unique$fill)
		# intersect.eems3.unique.plot  <- sp::plot(eems3.intersect.spdf.unique,col=eems3.intersect.spdf.unique$fill)
		# intersect.eems4.unique.plot  <- sp::plot(eems4.intersect.spdf.unique,col=eems4.intersect.spdf.unique$fill)
		ne_plot_10                 <- sp::plot(spdf_world_10,add=T,border="black",lwd=2)
		ne_oceans.plot_10          <- sp::plot(spdf_oceans_10,add=T,border="black",lwd=2,col="white")
		if(!is.null(coords)){
			coordinates <- data.matrix(read.table(coords,header=F))
			plot(coordinates,add=T,pch=20)
		}
	}
}






