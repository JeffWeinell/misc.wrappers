#' @title fracD function
#' 
#' Calculates the fractal dimension of perimeter-area ratio, a metric of polygon complexity.
#' More info: http://www.umass.edu/landeco/teaching/landscape_ecology/schedule/chapter9_metrics.pdf
#' 
#' @param x SpatialPolygons object (see sp package)
#' @return Number in [1,2), with larger values returned for increasing complex polygons; a value of 1 is returned for simple euclidean shapes likes squares and circles.
#' @export fracD
fracD <- function(x){
	# area in meters squared
	x.area      <- geosphere::areaPolygon(x)
	# perimeter in meters
	x.perimeter <- geosphere::perimeter(x)
	(2*log10(x.perimeter))/(log10(x.area))
}