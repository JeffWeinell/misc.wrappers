#' @title eemsgg2raster function
#' 
#' Convert one of the four maps produced by make_eems_plots from class ggplot to class RasterLayer, and optionally write the raster as a geoTIFF file.
#' 
#' @param x One of the plots produced by make_eems_plots: plots$mrates01, plots$mrates02, plots$qrates01, plots$qrates02
#' @param file.out Path where raster should be saved. Use extension '.tif' to save as geoTIFF.
#' @return Onject of class RasterLayer
#' @export eemsgg2raster
eemsgg2raster <- function(gg.obj,file.out=NULL){
	### data frame with xy coordinates and hex fill color code for each plot
	# gg.data.mrates1 <- ggplot2::ggplot_build(gg.obj$mrates01)[[1]][[1]]
	# gg.data.mrates2 <- ggplot2::ggplot_build(gg.obj$mrates02)[[1]][[1]]
	# gg.data.qrates1 <- ggplot2::ggplot_build(gg.obj$qrates01)[[1]][[1]]
	# gg.data.qrates2 <- ggplot2::ggplot_build(gg.obj$qrates02)[[1]][[1]]
	# ### data frame for each map holding colors in red,green,blue columns
	# rgb.mat.mrates1 <- t(col2rgb(gg.data.mrates1$fill))
	# rgb.mat.mrates2 <- t(col2rgb(gg.data.mrates2$fill))
	# rgb.mat.qrates1 <- t(col2rgb(gg.data.qrates1$fill))
	# rgb.mat.qrates2 <- t(col2rgb(gg.data.qrates2$fill))
	# ### rename columns
	# colnames(rgb.mat.mrates1) <- paste0(c("red","green","blue"),".m1")
	# colnames(rgb.mat.mrates2) <- paste0(c("red","green","blue"),".m2")
	# colnames(rgb.mat.qrates1) <- paste0(c("red","green","blue"),".q1")
	# colnames(rgb.mat.qrates2) <- paste0(c("red","green","blue"),".q2")
	# ### combine into one data frame
	# df <- cbind(gg.data.mrates1[,c("x","y")], rgb.mat.mrates1, rgb.mat.mrates2, rgb.mat.qrates1, rgb.mat.qrates2)
	gg.data  <- ggplot2::ggplot_build(gg.obj)[[1]][[1]]
	## data frame with red, green, blue color columns for each coordinate
	rgb.mat  <- t(col2rgb(gg.data$fill))
	df       <- cbind(gg.data[,c(2,3)],rgb.mat)
	sp::coordinates(df) <- ~ x + y
	sp::gridded(df)     <- TRUE
	# rasterDF  <- raster::raster(df)
	brickDF  <- raster::brick(df)
	if(!is.null(file.out)){
		raster::writeRaster(brickDF,filename=file.out)
	}
	rasterDF
}
