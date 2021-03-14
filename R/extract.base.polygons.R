#' @title extract.base.polygons function
#' 
#' Extract all lowest-level polygons of a SpatialPolygonsDataFrame
#' 
#' @param spdf SpatialPolygonsDataFrame (see sp package)
#' @return SpatialPolygons object with each lowest-level polygon of spdf held as a feature.
#' @export extract.base.polygons
extract.base.polygons <- function(spdf){
	for(i in 1:nrow(spdf)){
		sp.temp <- spdf@polygons[[i]]
		for(j in 1:length(sp.temp@Polygons)){
			id.temp <- paste0(i,".",j)
			polygons.ij <- sp::SpatialPolygons(list(sp::Polygons(list(sp.temp@Polygons[[j]]),ID=id.temp)))
			if(i==1 & j==1){
				polygons.all <- polygons.ij
			} else {
				polygons.all <- raster::bind(polygons.all,polygons.ij)
			}
		}
	}
	### Preserve crs definition
	raster::crs(polygons.all) <- raster::crs(spdf)
	### object returned
	polygons.all
}