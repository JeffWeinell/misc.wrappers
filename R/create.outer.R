#' @title create.outer function
#' 
#' This function automatically generates the '*.outer' file needed by EEMS, given the '*.coords' file also needed by EEMS.
#' 
#' @param coords Two column numeric matrix containing longitude and latitude of the samples in decimal degree format, or a character string with path to the input .coords file used for eems; the input file must be space separated and without a header.
#' @param method A number (1, 2, or 3) determining how the output coordinates will be generated.
#' If method=1 (the default) produces a potentially more useful result, and uses the alpha hull of geographic regions that intersect with input coordinates. This method may be slow if input coordinates include large sampling gaps. In such cases increases the value of alpha.start.
#' If method=2 is similar to method 1, except that the output coordinates are only determined from the input coordinates (and buffer setting), rather than from the regions of the base map (political geography from natural earth dataset) that intersect the input coordinates. Either method should be fine for EEMS, but method 1 tends to look nicer.
#' If method=3 is fast and uses the minimum convex polygon of input coodinates as a starting point.
#' @param buffer.adj This argument is not yet implemented. A number between 0 and 1 that determines the size of the buffer region surrounding the points supplied by coords.
#' @param coords.radius Number specifying the radius in degrees to use around input coordinates when determining spatial intersections (default 0.01). Increasing this value can be useful if coordinates are located offshore of regions outlined in the geography basemap.
#' @param max.fractal.dimension Value between 1 and 2 contolling the shape complexity of the polygon with the output coordinates. Values near 1 have low complexity (e.g., circles, squares), and shape complexity increases as max.fractal.dimension approaches 2. The default value 1.1 seems to work well for most cases. Ignored if method=3.
#' @param plot.outer Logical indicating if the points defining the output ".outer" table should be plotted.
#' @param ask.use Logical indicating if, after plotting, the user should be asked if the output should be written to file. Default is FALSE. Ignored if output.path is NULL or if plot.outer is FALSE.
#' @param counter.clockwise Should the output cooordinates be ordered counterclockwise. Default is TRUE and is the required format for EEMS.
#' @param output.path Character string with path where to save coordinates. Default is NULL.
#' @param plot.output.path Optional character string with path to save plot. Default is NULL. Ignored if plot.outer is FALSE.
#' @return A two column numerical matrix containing the longitude and latitude of the output polygon. The matrix written to output.path can be used as the ".outer" polygon used by EEMS. A map is plotted to visualize the results.
#' @export create.outer
create.outer <- function(coords,method=1,buffer.adj=0,coords.radius=0.01,max.fractal.dimension=1.1,plot.outer=TRUE,ask.use=FALSE,counter.clockwise=TRUE,output.path=NULL,plot.output.path=NULL){
	input.coords <- coords
	### Check which type of object is being supplied to coords and define coords accordingly
	if(is(input.coords,"character")){
		coords   <- data.matrix(read.table(input.coords,header=F))
	}
	if(is(input.coords,"matrix") | is(input.coords,"data.frame")){
		coords       <- data.matrix(input.coords)
		mode(coords) <- "numeric"
	}
	## High resolution global map of political regions.
	spdf_world_10             <- rnaturalearth::ne_countries(scale=10)
	## Character string defining the Coordinate Reference System (CRS) in WKT format
	# crs.string <- rgdal::showWKT(sp::proj4string(spdf_world_10))
	# crs.string <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
	crs.string <- "+init=EPSG:4326"
	# wkt.string <- rgdal::showWKT(crs.string)
	### Define CRS using wkt format
	suppressWarnings(raster::crs(spdf_world_10)   <- crs.string)
	### SpatialPoints object holding coords
	points.sp <- sp::SpatialPoints(coords)
	suppressWarnings(raster::crs(points.sp) <- crs.string)
	### SpatialPolygons object holding circles with centerpoints at coords
	circles.poly <- coords2sp.poly(coords.mat=coords,r=coords.radius)
	suppressWarnings(raster::crs(circles.poly) <- crs.string)
	### Spatial data frame for the map features that intersect with the input coords
	unique.features.at.circles <- do.call(rbind,unique(sp::over(circles.poly, spdf_world_10,returnList=T)))
	### spdf row indices for the unique regions in unique.features.at.circles
	indices_unique.features.at.circles <- which(spdf_world_10$sov_a3 %in% unique(unique.features.at.circles$sov_a3))
	### Corresponding SpatialPolygonsDateFrame object
	spdf_unique.features.at.circles <- spdf_world_10[indices_unique.features.at.circles,]
	### Restructuring the heirarchy of lowest-level polygons, treating them as separate features of one SpatialPolygons object.
	## The purpose of this is to be able to use sp::over() to search for polygons that intersect with coords.
	base.polygons.sp <- extract.base.polygons(spdf_unique.features.at.circles)
	### Get the feature indices for the polygons that intersect with circles representing input coordinates
	indices_polygons.at.circles <- sp::over(circles.poly,base.polygons.sp)
	### Check if any NA values in indices_unique.polygons.at.circles, which occurs when one more of the input circles does not intersect any polygons (i.e. they fall in the ocean).
	### If NAs exist, try increasing the value of r parameter of coords2sp.poly function.
	if(any(is.na(indices_polygons.at.circles))){
		problem.points <- unname(which(is.na(indices_polygons.at.circles)))
		stop(paste0("sample(s) coords[",paste(problem.points,collapse=" "),",]  does not intersect map. Try increasing coords.radius or remove sample(s) from input"))
	}
	### Extract polygons in base.polygons.sp that intersect circles.poly
	sp_unique.polygons.at.circles      <- base.polygons.sp[unique(indices_polygons.at.circles)]
	sp_unique.polygons.without.circles <- base.polygons.sp[-unique(indices_polygons.at.circles)]
	if(method %in% c(1,2)){
		#### Method 1 for getting maps.coords.list
		### Two column numerical matrix of the coordinates of polygons that intersect input points
		coords.polygons.with.points <- unique(round(do.call(rbind,lapply(X=1:length(sp_unique.polygons.at.circles),FUN=function(x){sp2coords(sp::polygons(sp_unique.polygons.at.circles)[x])})),digits=3))
		###
		if(method==1){
			coords.for.alphahull <- coords.polygons.with.points
			#buffer.width <- ((-0.75)+buffer.adj)
		}
		if(method==2){
			coords.for.alphahull <- unique(coords)
			#buffer.width <- buffer.adj
		}
		### initiallizing the while loop
		alpha.lower <- 1
		alpha.upper <- 50
		i=1; areas.equal <- FALSE; alpha <- alpha.upper; alpha.tol <- 0.1
		### Run while loop to find a low value of alpha that also generates a simple polygon (not self-intersecting or filled with holes) when input coordinates are supplied to ahull.
		# A polygon is considered simple if supplying the ahull polygon to gBuffer with buffer.width=0 produces a polygon with the same area as the input polygon. This works because gBuffer produces simple polygons even if the result clips the input polygon.
		while(((alpha.upper-alpha.lower)>alpha.tol) & (i < 100)){
			hull.alpha         <- alphahull::ahull(coords.for.alphahull,alpha=alpha)
			hull.alpha.coords  <- rbind(hull.alpha$arcs[,1:2],hull.alpha$arcs[1,1:2])
			hull.alpha.sp      <- sp::SpatialPoints(hull.alpha.coords)
			hull.alpha.polygon <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(hull.alpha.sp)),ID="alpha.hull")))
			if(rgeos::gIsSimple(hull.alpha.polygon)){
				### Test hull.alpha.polygon is simple
				# geosphere::areaPolygon(hull.alpha.polygon) == geosphere::areaPolygon(gBuffer(hull.alpha.polygon,width=0))
				### Create a polygon encompassing hull.alpha.polygon, with size determined by the value of buffer.width
				# Note: Do not set the crs of hull.alpha.polygon because gBuffer function would introduce distortions associated with the projection. In reality this might not be be a big deal but I'm not sure.
				buffer.poly        <- rgeos::gBuffer(hull.alpha.polygon,width=((-0.75)+buffer.adj))
				### Same as buffer.poly except points are added at midpoints of edges longer than max_distance; process repeated until all edges less than. This can make it easier to visualize the border of the output polygon, but otherwise this isnt useful.
				outer.poly         <- smoothr::densify(buffer.poly,max_distance=(0.25*alpha))
				### Run this test before switching back to Euclidean space
				### Set CRS of outer.poly (now we are back in Euclidean space); this must be done before using function sp::over()
				suppressWarnings(raster::crs(outer.poly) <- crs.string)
				### Test that the output coordinates define a simple polygon.
				test.simple        <- all(c(rgeos::gIsSimple(outer.poly), rgeos::gIsSimple(hull.alpha.polygon)))
				### Test if all of the input points are in the output polygon
				test.coverage      <- (!any(is.na(sp::over(points.sp,outer.poly))))
			} else {
				test.simple   <- FALSE
				test.coverage <- FALSE
			}
			if(test.simple & test.coverage){
				alpha.upper <- min(alpha.upper,alpha)
			} else {
				alpha.lower <- max(alpha.lower,alpha)
			}
			### save a copy of the current alpha value to initialize the next while loop
			alpha.i <- alpha
			alpha   <- median(c(alpha.lower,alpha.upper))
			i=(i+1)
		}
		fractal.dimension  <- fracD(outer.poly)
		alpha          <- alpha.i
		while(fractal.dimension > max.fractal.dimension){
			alpha              <- (alpha+0.1)
			hull.alpha         <- alphahull::ahull(coords.for.alphahull,alpha=alpha)
			hull.alpha.coords  <- rbind(hull.alpha$arcs[,1:2],hull.alpha$arcs[1,1:2])
			hull.alpha.sp      <- sp::SpatialPoints(hull.alpha.coords)
			hull.alpha.polygon <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(hull.alpha.sp)),ID="alpha.hull")))
			buffer.poly        <- rgeos::gBuffer(hull.alpha.polygon,width=((-0.75)+buffer.adj))
			outer.poly         <- smoothr::densify(buffer.poly,max_distance=(0.25*alpha))
			fractal.dimension  <- fracD(outer.poly)
		}
		suppressWarnings(raster::crs(outer.poly) <- crs.string)
	} else {
		if(method==3){
			### Minimum convex polygon that contains the maps with points.
			polygons.mcp <- suppressWarnings(adehabitatHR::mcp(points.sp,percent=100))
			mcp.coords   <- unique(sp2coords(polygons.mcp))
			### Generate the alpha-hull of the polygons with points
			hull.alpha         <- alphahull::ahull(mcp.coords,alpha=1)
			hull.alpha.coords  <- rbind(hull.alpha$arcs[,1:2],hull.alpha$arcs[1,1:2])
			hull.alpha.sp      <- sp::SpatialPoints(hull.alpha.coords)
			hull.alpha.polygon <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(hull.alpha.sp)),ID="alpha.hull")))
			buffer.width       <- (0.5+buffer.adj)
			### Create a polygon encompassing hull.alpha.polygon, with size determined by the value of buffer.width
			buffer.poly        <- rgeos::gBuffer(hull.alpha.polygon,width=buffer.width)
			### Same as buffer.poly except points are added at midpoints of edges longer than max_distance; process repeated until all edges less than. This can make it easier to visualize the border of the output polygon, but otherwise this isnt useful.
			outer.poly         <- smoothr::densify(buffer.poly,max_distance=(0.25*alpha))
			suppressWarnings(raster::crs(outer.poly) <- crs.string)
		} else {
			stop("set method to 1, 2, or 3")
		}
	}
	outer     <- unname(sp2coords(outer.poly))
	if(counter.clockwise){
		outer <- cbind(rev(outer[,1]),rev(outer[,2]))
	}
	### Show the relationship between the outer coordinates, the input coordinates, and geographic regions
	if(plot.outer){
		### Turns off current graphics devices before plotting. This seems to help prevent the new plots from inheriting graphical parameters from the current device.
		#if(!is.null(dev.list())){
		#	dev.off()
		#}
		### Plots in slightly different ways depending depending on the method used.
		if(!is.null(plot.output.path)){
			pdf(plot.output.path)
		}
		if(method %in% c(1,2)){
			### Plot the outer coordinates in white to define the extent
			sp::plot(sp::SpatialPoints(outer),col="white")
			#plot(submaps.with.points.sp,add=T,col="gray90")
			sp::plot(spdf_world_10,add=T)
			if(method==1){
			### Add a layer containing the regions of the basemap that include at least one point (regions in gray)
				sp::plot(spdf_unique.features.at.circles,add=T,col="gray90")
				sp::plot(sp_unique.polygons.without.circles,add=T,col="white")
			}
			### Add green crosshairs for output coordinates (the outer file for eems)
			sp::plot(sp::SpatialPoints(outer),col="green",add=T)
			### Add black circles for locations of the input coordinates
			sp::plot(points.sp,add=T,pch=20)
			### Add text/legend to show what the green crosshairs are
			mtext("+ = coordinates written to output.path (if set)",col="green",adj=0,line=1.5,cex=0.8)
			mtext(paste0("P:A fractal dimension of area inscribed by output coordinates",round(fractal.dimension,digits=3)),adj=0,line=0.5,cex=0.8)
			### Add map axes
			maps::map.axes()
		} else {
			sp::plot(sp::SpatialPoints(outer),col="white")
			extent.sp <- sp::SpatialPoints(rbind(outer,coords))
			sp::plot(spdf_world_10,add=T)
			sp::plot(sp::SpatialPoints(outer),col="green",add=T)
			sp::plot(points.sp,add=T,pch=20)
			mtext("+ outer (coordinates written to output.path)",col="green",adj=0.1,line=0.25,cex=0.8)
			maps::map.axes()
		}
		if(!is.null(plot.output.path)){
			dev.off()
		}
	}
	### Write outer file to output.path
	if(!is.null(output.path)){
		# If ask.use and plot.outer are both TRUE, then prompt the user to confirm that they want to write the outer coordinates to output.path.
		if(ask.use & plot.outer){
			user.response <- readline(prompt="Use outer habitat coordinates? (yes/no): ")
			if(user.response %in% c("yes","no")){
				if(user.response == "yes"){
					write.table(outer,output.path,row.names=F,col.names=F,sep=" ",quote=F)
				}
			} else {
				return(stop("aborting. unrecognized response to promt; expected (yes/no) response"))
			}
		} else {
			write.table(outer,output.path,row.names=F,col.names=F,sep=" ",quote=F)
		}
	}
	### return a list including a matrix of the coordinates in the output file and a SpatialPolygons object of the output coordinates
	list(outer,outer.poly)
}




