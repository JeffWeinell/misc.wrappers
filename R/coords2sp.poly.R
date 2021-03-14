#' @title coords2sp.poly function
#' 
#' Convert a coordinates matrix into a SpatialPolygons object of circles centered at input coordinates
#' 
#' @param coords.mat Numerical matrix of coordinates with longitude and latitude columns, respectively; coordinates must be in decimal degrees format.
#' @param r Number indicating radius of circles, in degrees.
#' @return SpatialPolygons object (see sp package)
#' @export coords2sp.poly
coords2sp.poly <- function(coords.mat,r=0.01){
	colnames(coords.mat) <- c("x","y")
	circles.list <- apply(X=coords.mat,MARGIN=1,FUN=function(x){sampSurf::spCircle(radius=r,centerPoint=x)[[1]]})
	res          <- do.call(raster::bind,circles.list)
	res
}