#' @title sp2coords function
#' 
#' Function to extract the coordinates matrix from a class Polygon, Polygons, SpatialPolygons, or SpatialPolygonsDataFrame object that containins a single polygon.
#' Warning: Do not use this function if the inout spatial object contains more than coordinates matrix.
#' 
#' @param x sp::Spatial* class object ultimately determined by a single coordinates matrix.
#' @return Numerical matrix of coordinates with longitude and latitude columns, respectively; coordinates must be in decimal degrees format.
#' @export sp2coords
sp2coords <- function(x){
	if(is(x,"SpatialPolygonsDataFrame") | is(x,"SpatialPolygons")){
		result <- x@polygons[[1]]@Polygons[[1]]@coords
	}
	if(is(x,"Polygons")){
		result <- x@Polygons[[1]]@coords
	}
	if(is(x,"Polygon") | is(x,"SpatialPoints")){
		result <- x@coords
	}
	unname(result)
}