#' @title eemsgg2raster function
#' 
#' Convert one of the four maps produced by make_eems_plots from class ggplot to class RasterLayer, and optionally write the raster as a geoTIFF file.
#' 
#' @param x One of the plots produced by make_eems_plots: plots$mrates01, plots$mrates02, plots$qrates01, plots$qrates02
#' @param file.out Path where raster should be saved. Use extension '.tif' to save as geoTIFF.
#' @return Onject of class RasterLayer
#' @export eemsgg2raster
eemsgg2raster <- function(gg.obj,file.out=NULL){
	gg.data  <- ggplot2::ggplot_build(gg.obj)[[1]][[1]]
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
