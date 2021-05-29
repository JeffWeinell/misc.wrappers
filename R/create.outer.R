#' @title create.outer function
#' 
#' This function generates an '*.outer' file needed by EEMS, given the '*.coords' file also needed by EEMS.
#' 
#' @param coords Two column numeric matrix containing longitude and latitude of the samples in decimal degree format, or a character string with path to the input .coords file used for eems; the input file must be space separated and without a header.
#' @param method A number (1, 2, or 3) determining how the output coordinates will be generated.
#' If method=1  generates an alpha hull surrounding geographic regions geographic regions that intersect with input coordinates. This method tends to works best for island arichpelagos, but may be slow if input coordinates include large sampling gaps. In such cases increase the value of alpha.start.
#' If method=2 (the default) is similar to method 1, except that the output coordinates are only determined from the input coordinates (and buffer setting), rather than from the regions of the base map (political geography from natural earth dataset) that intersect the input coordinates. Either method should be fine for EEMS, but method 1 tends to look nicer for island archipelagos.
#' If method=3 is fast and uses the minimum convex polygon of input coodinates as a starting point.
#' @param buffer.adj This argument is not yet implemented. A number between 0 and 1 that determines the size of the buffer region surrounding the points supplied by coords.
#' @param coords.radius Number specifying the radius in degrees to use around input coordinates when determining spatial intersections (default 0.01). Increasing this value can be useful if coordinates are located offshore of regions outlined in the geography basemap.
#' @param max.fractal.dimension Value between 1 and 2 contolling the shape complexity of the polygon with the output coordinates. Values near 1 have low complexity (e.g., circles, squares), and shape complexity increases as max.fractal.dimension approaches 2. The default value 1.1 seems to work well for most cases. Ignored if method=3.
#' @param plot.outer Logical indicating if the points defining the output ".outer" table should be plotted. Default TRUE.
#' @param ask.use Logical indicating if, after plotting, the user should be prompted to ask if the output should be saved. Default is FALSE. Ignored if output.path is NULL or if plot.outer is FALSE.
#' @param counter.clockwise Should the output cooordinates be ordered counterclockwise. Default is TRUE and is the required format for EEMS.
#' @param output.path Character string with path where to save coordinates. Default is NULL.
#' @param plot.output.path Optional character string with path to save plot. Default is NULL, in which case the value is equal to appending '.pdf' to 'output.path'. Ignored if plot.outer is FALSE.
#' @return A two column numerical matrix containing the longitude and latitude of the output polygon. The matrix written to output.path can be used as the ".outer" polygon used by EEMS. A map is plotted to visualize the results.
#' @export create.outer
create.outer <- function(coords,method=2,buffer.adj=0,coords.radius=0.01,max.fractal.dimension=1.1,plot.outer=TRUE,ask.use=FALSE,counter.clockwise=TRUE,output.path=NULL,plot.output.path=NULL){
	usegg <- TRUE
	input.coords <- unname(coords)
	colnames(input.coords) <- c("latitude","longitude")
	### Check which type of object is being supplied to coords and define coords accordingly
	if(is(input.coords,"character")){
		coords   <- data.matrix(read.table(input.coords,header=F))
	}
	if(is(input.coords,"matrix") | is(input.coords,"data.frame")){
		coords       <- data.matrix(input.coords)
		mode(coords) <- "numeric"
	}
	### Data frame copy of coords
	coords.df <- data.frame(X=coords[,1],Y=coords[,2])
	### SpatialPoints object holding coords
	points.sp <- sp::SpatialPoints(coords)
	### SpatialPolygons object holding circles with centerpoints at coords
	circles.poly <- coords2sp.poly(coords.mat=coords,r=coords.radius)
	### coordinate reference system to use for all spatial objects.
	crs.string <- "+init=EPSG:4326"
	## High resolution global map of political regions.
	spdf_world_10  <- rnaturalearth::ne_countries(scale=10)
	## Same high resolution global map of political regions, but held as an SF object.
	sf_world_10  <- rnaturalearth::ne_countries(scale=10,returnclass="sf")
	## Character string defining the Coordinate Reference System (CRS) in WKT format
	# crs.string <- rgdal::showWKT(sp::proj4string(spdf_world_10))
	# crs.string <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
	# wkt.string <- rgdal::showWKT(crs.string)
	### Define CRS using wkt format
	suppressWarnings(raster::crs(spdf_world_10) <- raster::crs(points.sp) <- raster::crs(circles.poly) <- crs.string)
	#suppressWarnings(raster::crs(points.sp) <- crs.string)
	#suppressWarnings(raster::crs(circles.poly) <- crs.string)
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
	outer.df <- data.frame(X=outer[,1],Y=outer[,2])
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
			### ggplot version for Methods 1/2
			if(usegg){
				#gg.baseplot  <- ggplot2::ggplot(data = sf_world_10) + ggplot2::theme_classic() + ggplot2::geom_sf(fill="white",color="black") + ggplot2::coord_sf(xlim =xlim, ylim = ylim) + ggplot2::xlab("longitude") + ggplot2::ylab("latitude") + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1), plot.title=ggplot2::element_text(size=10, hjust = 0.5),legend.title=ggplot2::element_blank())
				gg.baseplot  <- ggplot2::ggplot(data = sf_world_10) + ggplot2::theme_classic() + ggplot2::geom_sf(fill="gray90",color="black") + ggplot2::xlab("longitude") + ggplot2::ylab("latitude") + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1), plot.title=ggplot2::element_text(size=10, hjust = 0.5),legend.title=ggplot2::element_blank())
				if(method==1){
					xlim <- rangeBuffer(x=outer[,1],f=0.15)
					ylim <- rangeBuffer(x=outer[,2],f=0.15)
					## Convert spatial polygons to SF object
					sf_unique.features.at.circles      <- sf::st_as_sf(spdf_unique.features.at.circles)
					#sf_unique.polygons.without.circles <- sf::st_as_sf(sp_unique.polygons.without.circles)
					### Add a layer containing the regions of the basemap that include at least one point (regions in gray)
					# May be a problem with inheritance of coordinate systems...
					gg.baseplot2 <- gg.baseplot + ggplot2::geom_sf(data=sf_unique.features.at.circles,color="black",fill="gray90") + ggplot2::coord_sf(xlim =xlim, ylim = ylim, expand=FALSE)
				} else {
					xlim <- rangeBuffer(x=outer[,1],f=0.25)
					ylim <- rangeBuffer(x=outer[,2],f=0.25)
					gg.baseplot2 <- gg.baseplot + ggplot2::coord_sf(xlim =xlim, ylim = ylim)
				}
				### Data frames holding input coordinates, outer coordinates, and both the input and outer coordinates.
				gg.coords.df  <- cbind(coords.df,grp= "input sample points",col="black",shp=20,strk=1,sz=3,order=1)
				gg.outer.df   <- cbind(outer.df,grp="outer points",col="green",shp=3,strk=1.25,sz=3,order=2)
				gg.df         <- rbind(gg.outer.df,gg.coords.df)
				gg.df         <- rbind(gg.coords.df,gg.outer.df)
				titletext     <- paste0("misc.wrappers::create.outer(method=",method,")")
				#gg.outerplot  <- gg.baseplot2 + ggplot2::geom_point(data = gg.df, ggplot2::aes(x = X, y = Y, fill=grp,shape=grp,color=grp,size=grp,stroke=grp)) + ggplot2::scale_discrete_manual(aesthetics = "stroke", values = c(1,1.25)) + ggplot2::scale_size_manual(values=c(3,3)) + ggplot2::scale_fill_manual(values=c("black","green")) + ggplot2::scale_color_manual(values=c("black","green")) + ggplot2::scale_shape_manual(values=c(20,3)) + ggplot2::ggtitle(titletext)
				gg.outerplot  <- gg.baseplot2 + ggplot2::geom_point(data = gg.df, ggplot2::aes(x=X, y=Y, fill=grp, shape=grp, color=grp, size=grp, stroke=grp,alpha=grp)) + ggplot2::scale_discrete_manual(aesthetics = "stroke", values = c(`input sample points`=1,`outer points`=1.25)) + ggplot2::scale_alpha_manual(values=c(`input sample points`=1,`outer points`=0.7)) + ggplot2::scale_size_manual(values=c(`input sample points`=3,`outer points`=3)) + ggplot2::scale_fill_manual(values=c(`input sample points`="black",`outer points`="green")) + ggplot2::scale_color_manual(values=c(`input sample points`="black",`outer points`="green")) + ggplot2::scale_shape_manual(values=c(`input sample points`=20,`outer points`=3)) + ggplot2::ggtitle(titletext)
				### Reverses the order of the items in the legend without affecting the map.
				# gg.outerplot2 <- gg.outerplot + guides(color = guide_legend(reverse=T), fill = guide_legend(reverse=T), shape=guide_legend(reverse=T), size=guide_legend(reverse=T), stroke=guide_legend(reverse=T))
				print(gg.outerplot)
			}
			### base plot version for Methods 1/2
			if(!usegg){
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
				mtext("+ = habitat (.*outer) coordinates",col="green",adj=0,line=2.5,cex=0.8)
				mtext(paste0("fractal dimension of habitat shape: ",round(fractal.dimension,digits=3)),adj=0,line=1.5,cex=0.8)
				mtext(paste0("input settings: method = ",method,", max.fractal.dimension = ",max.fractal.dimension,", buffer.adj = ",buffer.adj),adj=0,line=0.5,cex=0.8)
				### Add map axes
				maps::map.axes()
			}
		} else {
			### ggplot version for Method 3
			if(usegg){
				xlim <- rangeBuffer(x=outer[,1],f=0.4)
				ylim <- rangeBuffer(x=outer[,2],f=0.4)
				gg.baseplot  <- ggplot2::ggplot(data = sf_world_10) + ggplot2::theme_classic() + ggplot2::geom_sf(fill="gray90",color="black") + ggplot2::coord_sf(xlim =xlim, ylim = ylim) + ggplot2::xlab("longitude") + ggplot2::ylab("latitude") + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1), plot.title=ggplot2::element_text(size=10, hjust = 0.5),legend.title=ggplot2::element_blank())
				gg.coords.df  <- cbind(coords.df,grp= "input sample points")
				gg.outer.df   <- cbind(outer.df,grp="outer points")
			#	gg.df         <- rbind(gg.outer.df,gg.coords.df)
				gg.df         <- rbind(gg.outer.df,gg.coords.df)
				gg.outerplot  <- gg.baseplot + ggplot2::geom_point(data = gg.df, ggplot2::aes(x = X, y = Y, fill=grp,shape=grp,color=grp,size=grp,stroke=grp,alpha=grp)) + ggplot2::scale_size_manual(values=c(`input sample points`=3,`outer points`=3)) + ggplot2::scale_alpha_manual(values=c(`input sample points`=1,`outer points`=0.7)) + ggplot2::scale_discrete_manual(aesthetics = "stroke", values = c(`input sample points`=1,`outer points`=1.25)) + ggplot2::scale_fill_manual(values=c(`input sample points`="black",`outer points`="green")) + ggplot2::scale_color_manual(values=c(`input sample points`="black",`outer points`="green")) + ggplot2::scale_shape_manual(values=c(`input sample points`=20,`outer points`=3)) + ggplot2::ggtitle("misc.wrappers::create.outer(method=3)")
				# gg.outerplot2 <- gg.outerplot + guides(color = guide_legend(reverse=T), fill = guide_legend(reverse=T),shape=guide_legend(reverse=T),size=guide_legend(reverse=T))
				print(gg.outerplot)
			}
			### base plot version for Method 3
			if(!usegg){
				sp::plot(sp::SpatialPoints(outer),col="white")
				#return(rbind(outer,coords))
				extent.sp <- sp::SpatialPoints(unname(rbind(outer,coords)))
				#return(extent.sp)
				sp::plot(spdf_world_10,add=T)
				sp::plot(sp::SpatialPoints(outer),col="green",add=T)
				sp::plot(points.sp,add=T,pch=20)
				mtext("+ outer (coordinates written to output.path)",col="green",adj=0.1,line=1.5,cex=0.8)
				mtext(paste0("input settings: method = ",method),adj=0,line=0.5,cex=0.8)
				maps::map.axes()
			}
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
	#return(gg.outerplot)
	### return a list including a matrix of the coordinates in the output file and a SpatialPolygons object of the output coordinates
	#list(outer,outer.poly)
	list(outer,gg.outerplot)
} 
#' @examples
#' library(misc.wrappers)
#' # Sample 50 points from 10-degree radius area with center located on land somewhere between -50 and 50 degrees latitude.
#' coords50 <- rcoords(r=10,size=50,limits=c(-180,180,-50,50))
#' # Use method 1 to create the outer set of points and generate a figure with points on map. Method 1 may not work well for points on continents vs. island archipelagos.
#' coords50_outer1 <- create.outer(coords=coords50, method=1, output.path="coords50_outer1.txt", plot.output.path="coords50_outer1.pdf")
#' # Use method 2
#' coords50_outer2 <- create.outer(coords=coords50, method=2, output.path="coords50_outer2.txt", plot.output.path="coords50_outer2.pdf")
#' # Use method 3
#' coords50_outer3 <- create.outer(coords=coords50, method=3, output.path="coords50_outer3.txt", plot.output.path="coords50_outer3.pdf")

#' @title Return coordinates over land
#' 
#' Returns the set of input coordinates that intersect with land. Land is considered to incude all areas covered by a feature in the rnaturalearth::ne_countries 10 meter resolution dataset. Input coordinates are are assumed to be WGS84 (EPSG:4326)
#' 
#' @param x Matrix or data frame with longitude and latitude [decimal degree] in the first two columns, respectively.
#' @param return.as Character string with class to use for object returned. Default "data.frame". Can also be "matrix" or "SP" (SpatialPoints).
#' @return An object with class equal to the value of 'return.as', which contains the subset of input points that occur on land, or more specifically, not in an ocean.
#' @export points.on.land
points.on.land <- function(x,return.as="data.frame"){
	world.land <- rnaturalearth::ne_countries(scale=10)
	crs.string <- "+init=EPSG:4326"
	if(is(x,"matrix")){
		x2 <- as.data.frame(x)
	} else {
		if(is(x,"data.frame")){
			x2 <- x
		} else {
			stop("'x' must be a matrix or data.frame containing columns with longitude and latitude")
		}
	}
	x3 <- x2[,1:2]
	colnames(x3) <- c("longitude","latitude")
	mode(x3[,"longitude"]) <- "numeric"
	mode(x3[,"latitude"]) <- "numeric"
	xsp <- sp::SpatialPoints(x3)
	suppressWarnings(raster::crs(world.land) <- raster::crs(xsp) <- crs.string)
	matches.df <- sp::over(xsp, world.land)
	landtest   <- apply(matches.df,MARGIN=1,FUN=function(m){!all(is.na(m))})
	if(any(landtest)){
		res <- x2[which(landtest),]
	} else {
		return(NULL)
	}
	if(return.as=="data.frame"){
		return(res)
	}
	if(return.as=="matrix"){
		return(data.matrix(res))
	}
	if(return.as=="SP"){
		return(xsp[which(landtest)])
	}
}

#' @title Sample points in a random-centered area or tethered group of areas.
#' 
#' This function provides a way to randomly generate one or more groups of points, such that groups are non-overlapping, and options for controlling proximity, group sizes, and optional requirement that points occur over land.
#' At present areas are circles projected onto a euclidean plane.
#' 
##' @param r Number specifying the radius of group areas, in decimal degrees, or a numerical vector with min and max values of uniform sampling distribution from which radii of areas will be drawn.
#' @param regionsize Number between 0 and 1 specifying the fraction of the possible sampling area (the window region set by 'wnd' argument) in which points may be distributed. Default 0.25. When 'wnd' is default (entire Earth), regionsize is default (0.25), and ('over.land' condition is FALSE), then minimum convex hull of all sampled points is expected to cover 1/4 of Earth.
#' @param samplesize Total number of points to sample. Default 100.
#' @param n.grp Number of groups to sample (default 1). Each group is defined by its own spatial polygon from which samples are drawn. The center of each group's sampling polygon falls within a polygon with size = ('regionsize' * area of 'wnd').
#' @param grp.n.weights Probability weights for group sample sizes of groups. Default is equal weights. The sum of points sampled for each group = 'samplesize' argument.
#' @param grp.area.weights Probability weights determing group region sizes. Default is equal weights.
#' @param wnd Either a character string or vector describing one or more regions of Earth, or a length four numerical vector specifying the bounding box (longitude and latitude ranges) for the region where sampling is allowed. Default is to include all of Earth c(-180,180,-90,90). 
#' FUTURE OPTION: 'wnd' is a matrix with four columns, with each row describing the bounding box for one of the populations.
#' Character strings can be one of the following (NOT ENTIRELY IMPLEMENTED): "Earth", "tropics", "middleSouth", "middleNorth" "middle", "Arctic", "Antarctic", "Polar", "NorthernHemisphere", "WesternHemisphere", or "EasternHemisphere".
#' @param over.land Logical indicating whether or not all returned points must occur over land. Default TRUE. Note that this condition is only applied to samples returned as output. Therefore, 'samplesize' is really a 'sample until' rule. This is faster than performing clipping operations on proposed sample polygons to conform to geography.
#' @param interactions Numeric vector with the minimum and maximum amount of overlap between each pair of groups (aka populations), calculated as (intersect area)/(minimum of non-intersected area for each group). Default c(0,1) allows for all possible scenarios. Examples: c(0,0) specifies that groups must be allopatric; c(1,1) requires complete overlap of groups, which is not realistic given the stochasticity determining region sizes; c(0.5,1) requires that at least half-overlaps between groups; c(0.2,0.25) specifies a small contact zone.
##' @param d.grp Number controlling the distance between the centers of a pair of areas, as a function of the pair's radii. Default 1, which would allow sample regions to nearly coincide. Future option may allow for a pairwise distance matrix.
##' @param expf Number that affects dispersion relative to regionsize. Higher numbers increase dispersion. Default 8.
#' @param min.grp.size Minimum number of samples required for each group. Default 2.
##' @param grp.scaler Number > 0 that scales geographical areas of all groups
#' @param return.as Character string with class to use for object returned. Default "data.frame". Can also be "matrix" or "SP" (SpatialPoints).
#' @param show.plot Whether or not the points should be plotted on a low-resolution land map. The map is used is the rnaturalearth countries map, 110 meter resolution.
#' @return An object with class equal to the value of 'return.as' and containing the set of points that meet the specified sampling requirements. If return.as='matrix' or 'data.frame', the columns are 'X' (for longitude), 'Y' (for latitude), and 'group' (all 1 if 'n.grp'=1).
#' @export rcoords
rcoords <- function(regionsize=0.25, samplesize=100, n.grp=1, grp.n.weights=NULL, grp.area.weights=NULL, min.grp.size=2, wnd= c(-180, 180,-90, 90), over.land=TRUE, interactions=c(0,1), show.plot=FALSE, return.as="data.frame"){
	if(is.null(grp.n.weights)){
		grp.n.weights <- rep(1,n.grp)
	}
	if(is.null(grp.area.weights)){
		grp.area.weights <- rep(1,n.grp)
	}
	# Number > 0 that scales geographical areas of all groups relative to 'regionsize'*(area of 'wnd')
	grp.scaler=(1/n.grp)
	# Adding this soon. Plan is to convert each bounding box to an extent object and then Polygon object, and then combine into SpatialPolygons for multi-part areas.
	if(FALSE){
		earth.bb         <- c(-180, 180,-90, 90)
		tropics.bb       <- c(-180,180,-23.43651,23.43651)
		south.middle.wnd <- c(-180,180,-23.439444,-66.560833)
		north.middle.wnd <- c(-180,180,0,66.560833)
		arctic.wnd    <- c()
		antarctic.wnd <- c()
		western.hemisphere  <- c()
		eastern.hemisphere  <- c()
		northern.hemisphere <- c()
		southern.hemisphere <- c()
		# Continents and Bioregions...
	}
	# Defaults: n.grp=2; samplesize=100; regionsize=0.25; grp.n.weights=rep(1,n.grp); grp.area.weights=rep(1,n.grp); wnd= c(-180, 180,-90, 90); over.land=TRUE; interactions=c(0,1); expf=8; show.plot=FALSE;  return.as="data.frame"; grp.scaler=1
	grp.scaler.start <- grp.scaler
	#result.temp        <- data.frame(NULL)
	PASS=FALSE
	miter <- 1000
	ctr   <- 0
	if(length(samplesize)==1){
		if((n.grp*min.grp.size) > samplesize){
			stop("Total sample size 'samplesize' must be greater than the number of groups 'n.grp' times minimum group size 'min.grp.size'.",call.=FALSE)
		}
	} else {
		if(length(samplesize)!=n.grp){
			stop("If 'samplesize' is a numerical vector, then it must have a length equal to the number of groups 'n.grp'.",call.=FALSE)
		}
	}
	# Find group centers over land
	#if(over.land){
	#	intial <- data.frame(NULL)
	#	while(nrow(intial) < size){
	#		intial <- data.frame(NULL)
	#		### Sample a point on earth to use as the center of the sampling circle
	#		sample.center  <- c(x=sample(seq(limits[1],to=limits[2],by=0.01),size=1), y=sample(seq(limits[3],to=limits[4],by=0.01),size=1))
	#		
	#	}
	#}
	if(is(wnd,"matrix")){
			stop("not yet implemented")
			### in the future, each row will describe the bounding box of a group
		} else {
			if(is(wnd,"SpatialPolygons")){
				if(FALSE){
					stop("not yet implemented")
					### In the future, a spatial polygons object can be supplied
					world.land.sp    <- rnaturalearth::ne_countries(scale=110, returnclass="sp")
					ph.land10.sp     <- world.land10.sp[which(world.land10.sp[,"sovereignt"][[1]]=="Philippines"),]
					ph.island.areas  <- sapply(X=1:length(ph.land10.sp@polygons[[1]]@Polygons),FUN=function(x){ph.land10.sp@polygons[[1]]@Polygons[[x]]@area})
					Luzon.polygon    <- ph.land10.sp@polygons[[1]]@Polygons[[which(ph.island.areas == max(ph.island.areas))]]
					Luzon.sp         <- sp::SpatialPolygons(list(sp::Polygons(list(Luzon.polygon),1)))
					# Would set 'wnd' argument to Luzon.sp
					wnd <- Luzon.sp
					# wnd.sp is simply wnd
				}
				wnd.sp <- wnd
				raster::crs(wnd.sp)  <- sp::CRS("EPSG:4326")
				wnd.area             <- wnd.sp@polygons[[1]]@Polygons[[1]]@area
				regionsize.dd0       <- (regionsize*wnd.area)
				#regionsize.dd        <- (regionsize.dd0*0.25)
				regionsize.dd        <- (regionsize.dd0)
			} else {
			if(is(wnd,"numeric")){
				# rectangular SpatialPolygons object with extent formed by 'wnd' argument.
				wnd.sp <- as(raster::extent(wnd), "SpatialPolygons")
				# area of the region bounded by 'wnd', in decimal degrees
				wnd.sp2 <- wnd.sp
				#sp::proj4string(wnd.sp2) <- sp::CRS("+proj=longlat +datum=WGS84")
				raster::crs(wnd.sp2) <- sp::CRS("EPSG:4326")
				#wnd.area <- attributes(attributes(wnd.sp2)$polygons[[1]])$area
				wnd.area <- wnd.sp2@polygons[[1]]@Polygons[[1]]@area
				# area of region spanning all group sampling areas, in decimal degrees
				regionsize.dd0 <- (regionsize*wnd.area)
				regionsize.dd  <- (regionsize.dd0*0.25)
			}
		}
	}
	# wnd.sp3              <- sp::spTransform(wnd.sp2,sp::CRS("ESRI:54009"))
	# xymeters    <- dfTransform(df=data.frame(attributes(attributes(attributes(wnd.sp3)$polygons[[1]])[[1]][[1]])$coords),CRS2="ESRI:54009")
	# xy.sp       <- as(extent(xymeters),"SpatialPolygons")
	# raster::crs(xy.sp) <- sp::CRS("ESRI:54009")
	# wnd.area.m2 <- attributes(attributes(xy.sp)$polygons[[1]])$area
	if(FALSE){
		# area of the region defined by 'wnd'; when 'wnd' includes entire earth, seems to be underestimated by about 20% (measure=geodesic or haversine), 
		pts.corners <- bbox2points(rbind(x=c(min=-180, max=180),y=c(min=-90,max=90)))
		sp::Polygon(pts.corners)
		ydist.km    <- geodist::geodist(x=t(pts.corners)[,1], y=t(pts.corners)[,2],measure="geodesic")/1000
		xdist.km    <- geodist::geodist(x=t(pts.corners)[,3], y=t(pts.corners)[,4],measure="geodesic")/1000
		area.km2    <- (ydist.km*xdist.km)
	}
	
	while(!PASS){
		ctr <- ctr + 1
		if(ctr > miter){
			stop(paste("after",miter,"attemps, failed to find points that pass conditions. ctr=",ctr,"ctr2=",ctr2,"; grp.scaler=",grp.scaler))
		}
		#if(ctr %in% seq(from=10,to=miter,by=10)){
		if(ctr>1){
			grp.scaler <- (grp.scaler*0.95)
		}
		### Clear the data frame if too few points on land during previous attempt
		#result.temp         <- data.frame(NULL)
		initial <- FALSE
		### Find group centers that pass conditions
		ctr2 <- 0
		while(!initial){
			ctr2 <- ctr2 + 1
			if(ctr2 > (miter/10)){
				break
				# stop(paste("failed to initialize after",miter,"attemps. ctr=",ctr,"ctr2=",ctr2))
			}
			if(ctr2>1){
				grp.scaler <- (grp.scaler*0.99)
			}
			# Random point within wnd.sp to determine subregion within which group centers will occur
			center.sp <- suppressWarnings(sp::spsample(wnd.sp, n=1, type="random"))
			# Length two vector with "x" and "y" entries from sample.center, to pass to spCircle
			center.xy <- c(sp::coordinates(center.sp)); names(center.xy) <- c("x","y")
			if(is(wnd,"SpatialPolygons")){
				sample.area.sp <- wnd.sp
			} else {
				sample.area.sp <- sampSurf::spCircle(radius=sqrt((regionsize.dd/pi)), centerPoint=center.xy)[[1]]
			}
			if(n.grp>1){
				# Sample group sampling-center locations from within sample.area.sp
				grp.pts <- suppressWarnings(sp::spsample(sample.area.sp,n=n.grp,type="random"))
				# coordinates extracted from grp.pts
				grp.pts.mat <- sp::coordinates(grp.pts)
				colnames(grp.pts.mat) <- c("x","y")
				grp.areas <- (regionsize.dd*(grp.area.weights/(sum(grp.area.weights)))*grp.scaler)
				# create polygons centered at grp.pts, from which each groups points will be drawn. First need to determine group sizes...
				grp.areas.sp <- lapply(X=1:n.grp, FUN=function(x){sampSurf::spCircle(radius=sqrt((grp.areas/pi))[x],centerPoint=c(grp.pts.mat[x,]))[[1]]})
				# Distance between group centers. Not sure if this is necessary.
				# geodist(grp.pts.mat,measure="geodesic")
				grp.center.dists0 <- geodist::geodist(grp.pts.mat,measure="geodesic")
				grp.center.dists  <- grp.center.dists0[lower.tri(grp.center.dists0)]
				# check that each group centers occur over land, if over.land is true
			} else {
				grp.pts.mat  <- matrix(center.xy,nrow=1,dimnames=list(c(NULL),c("x","y")))
				grp.areas.sp <- list(sample.area.sp)
			}
			if(over.land){
				if(length(points.on.land(x=grp.pts.mat)[,1]) < n.grp){
				#if(length(points.on.land(x=sp::coordinates(grp.pts))[,1]) < n.grp){
					next
				}
			}
			initial <- TRUE
		}
		
		### For each group, the number of samples that will be drawn (or saved, if over.land=TRUE) from the group's sample region.
		if(length(samplesize)==1){
			if(n.grp==1){
				grp.sizes <- samplesize
			} else {
				# Initial values for entering while loop
				grp.sizes <- rep(0,n.grp)
				while(any(grp.sizes < min.grp.size)){
					grp.sizes  <- c(table(sample(x=c(1:n.grp),size=samplesize,prob=grp.n.weights,replace=TRUE)))
					# verify that at least each group has at least min.grp.size.
				}
			}
		} else {
			grp.sizes <- samplesize
		}
		# grp.sizes <- round(samplesize * (grp.n.weights/sum(grp.n.weights)))
		### Creates a SpatialPolygons object (sampling region) for each group
		#sample.area.sp <- lapply(1:n.grp, FUN=function(x){sampSurf::spCircle(radius= sqrt((grp.areas/pi))[x] , centerPoint=c(x=group.centers[x,1],y=group.centers[x,2]))[[1]]})
		### Check if groups meet conditions set by 'interactions' argument.
		if(!is.null(interactions) & n.grp > 1){
			if(!(interactions[1]==0 & interactions[2]==1)){
				# Much less complex list holding the sample sampling polygons
				# polygons.list <- lapply(X=1:n.grp, function(x){attributes(attributes(grp.areas.sp[[x]])$"polygons"[[1]])[[1]][[1]]})
				# Possibly equivalent:
				polygons.list <- lapply(X=1:n.grp, function(x) grp.areas.sp[[x]]@polygons[[1]]@Polygons[[1]])
				# All possible pairwise group comparisons
				grp.pairs <- do.call(rbind,pset(1:n.grp,2,2))
				# Calculate intersect area yet and sum areas for each pair of groups. (Doesnt work yet).
				grp.int.area <- grp.sum.area <- grp.ratio <- c()
				for(i in 1:nrow(grp.pairs)){
					polygon.pair   <- polygons.list[c(grp.pairs[i,1],grp.pairs[i,2])]
					grp1.area.temp <- polygon.pair[[1]]@area
					grp2.area.temp <- polygon.pair[[2]]@area
					pairsum.temp   <- sum(c(grp1.area.temp, grp2.area.temp))
					#grp.int <- raster::intersect(sample.area.sp[[grp.pairs[i,1]]], sample.area.sp[[grp.pairs[i,2]]])
					# returns NULL if polygons do not intersect; returns SpatialPolygons object if they do intersect
					grp.int  <- rgeos::gIntersection(spgeom1=grp.areas.sp[[grp.pairs[i,1]]],spgeom2=grp.areas.sp[[grp.pairs[i,2]]])
					#Nonintersecting portions of areas 1 and 2 (faster than getting diff of each in the other and then adding)
					#grp.diff <- rgeos::gSymdifference(spgeom1=grp.areas.sp[[grp.pairs[i,1]]],spgeom2=grp.areas.sp[[grp.pairs[i,2]]])
					# returns a spatial polygons object with the non-intersecting portions of the group pair polygons each as their own polygons. If no overlap the returned polygons are the same the input polygons.
					grp.diff1 <- rgeos::gDifference(spgeom1=grp.areas.sp[[grp.pairs[i,1]]],spgeom2=grp.areas.sp[[grp.pairs[i,2]]])
					grp.diff2 <- rgeos::gDifference(spgeom1=grp.areas.sp[[grp.pairs[i,2]]],spgeom2=grp.areas.sp[[grp.pairs[i,1]]])
					if(is.null(grp.int)){
						int.area.temp <- 0
					} else {
						#int.area.temp <- attributes(attributes(attributes(grp.int)$"polygons"[[1]])[[1]][[1]])$area
						int.area.temp <- grp.int@polygons[[1]]@Polygons[[1]]@area
					}
					# Area of the region of group 1 that does not intersect group 2 (for this pair of groups).
					if(is.null(grp.diff1)){
						diff1.area.temp <- 0
					} else {
						#diff1.area.temp <- attributes(attributes(attributes(grp.diff1)$"polygons"[[1]])[[1]][[1]])$area
						diff1.area.temp <- grp.diff1@polygons[[1]]@Polygons[[1]]@area
					}
					if(is.null(grp.diff2)){
						diff2.area.temp <- 0
					} else {
						#diff2.area.temp <- attributes(attributes(attributes(grp.diff2)$"polygons"[[1]])[[1]][[1]])$area
						diff2.area.temp <- grp.diff2@polygons[[1]]@Polygons[[1]]@area
					}
					#if(is.null(grp.diff)){
					#	diff.area.temp <- 0
					#} else {
					#	diff1.area.temp <- attributes(attributes(attributes(grp.diff)$"polygons"[[1]])[[1]][[1]])$area
					#	diff2.area.temp <- attributes(attributes(attributes(grp.diff)$"polygons"[[1]])[[1]][[2]])$area
					#}
					diff.area.temp  <- diff1.area.temp + diff2.area.temp
					ratio.temp <- int.area.temp/diff.area.temp
					#ratio.temp <- (int.area.temp/pairsum.temp)
					#grp.int <- raster::intersect(spgeom1=sample.area.sp[[1]],spgeom2=sample.area.sp[[2]])
					grp.int.area <- c(grp.int.area,int.area.temp)
					grp.sum.area <- c(grp.sum.area,pairsum.temp)
					grp.ratio    <- c(grp.ratio,ratio.temp)
				}
				area.df <- data.frame(group.i=grp.pairs[,1],group.j=grp.pairs[,2],area.intersection=grp.int.area,area.sumpair=grp.sum.area,area.intersection.sum.ratio=grp.ratio)
				# Check that the required conditions for interactions are met
				if(!all(area.df$area.intersection.sum.ratio >= interactions[1] & area.df$area.intersection.sum.ratio <= interactions[2])){
					# If conditions not met, try again from scratch. In future versions it would be useful to slide the polygons along their center axes to make minor adjustments.
					next
				}
			}
		}
		
		### Method using st_sample from SF package. Probably not as useful as the spsample from SP.
		if(FALSE){
			### Convert each SpatialPolygons object to an SF object
	#		sample.area.sf <- sf::st_as_sf(sample.area.sp)
			grp.areas.sf <- lapply(X=1:n.grp,FUN=function(x){sf::st_as_sf(grp.areas.sp[[x]])})
			### Sample within each area. If points must be over land, sample four times as many points as requested by 'size' argument.
			if(over.land){
				nsamp <- (grp.sizes*4)
				samples.sf.temp  <- lapply(X=1:n.grp,FUN=function(x){sf::st_sample(x=grp.areas.sf[[x]], size= nsamp[x])})
			} else {
				samples.sf.temp  <- lapply(X=1:n.grp,FUN=function(x) {sf::st_sample(x= grp.areas.sf[[x]],size=grp.sizes[x])})
			}
			### Hold sampled coordinates as a list of data matrices
			samples.mat.temp <- lapply(X=1:n.grp,FUN=function(x){sf::st_coordinates(samples.sf.temp[[x]])})
		}
		### Method using spsample from SP package.
		if(over.land){
			nsamp <- (grp.sizes*4)
			samples.sp.temp  <- lapply(X=1:n.grp,FUN=function(x){sp::spsample(x=grp.areas.sp[[x]], n= nsamp[x], type="random")})
			for(i in 1:length(samples.sp.temp)){
				raster::crs(samples.sp.temp[[i]]) <- sp::CRS("EPSG:4326")
			}
		} else {
			samples.sp.temp  <- lapply(X=1:n.grp,FUN=function(x) {sp::spsample(x= grp.areas.sp[[x]],n=grp.sizes[x], type="random")})
			for(i in 1:length(samples.sp.temp)){
				raster::crs(samples.sp.temp[[i]]) <- sp::CRS("EPSG:4326")
			}
		}
		samples.mat.temp <- lapply(X=1:n.grp,FUN=function(x){sp::coordinates(samples.sp.temp[[x]])})
		### Hold sampled coordinates as a list of data frames, and include a column in each data frame to indicate group assignment
		samples.df.temp <- do.call(rbind,lapply(X=1:n.grp,FUN=function(x){data.frame(X=samples.mat.temp[[x]][,1],Y=samples.mat.temp[[x]][,2],group=x)}))
		
		# Filtering samples outside of 'wnd'; 1=sample is within wnd; NA=sample is outside of wnd
		samples.in.wnd <- unlist(lapply(1:n.grp,function(x) over(samples.sp.temp[[x]],wnd.sp)))
		#return(samples.in.wnd)
		if(any(samples.in.wnd==1,na.rm=TRUE)){
			samples.df.temp <- samples.df.temp[which(samples.in.wnd==1),]
		}
		### Hold sampled coordinates as a list of data frames, and include a column in each data frame to indicate group assignment
#		samples.df.temp <- lapply(X=1:n.grp,FUN=function(x){data.frame(X=samples.mat.temp[[x]][,1],Y=samples.mat.temp[[x]][,2],group=x)})
		####### PAUSED HERE ##########

		
		# ggplot2::ggplot(data = world.land) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = samples.df.temp, ggplot2::aes(x = X, y = Y), size = 2, shape = 21, fill = "darkred") # + coord_sf(xlim = c(-90, -78), ylim = c(24.5, 40), expand = FALSE)
		if(over.land){
			### Get subset of sampled coordinates that fall on land
			#samples.temp <- lapply(X=1:n.grp,FUN=function(x){points.on.land(x=samples.mat.temp[[x]])})
			# filtering samples that are not on land
			samples.temp <- points.on.land(samples.df.temp)
		#	samples.temp <- lapply(X=1:n.grp,FUN=function(x){points.on.land(x=samples.df.temp[[x]])})
			#samples.land.temp <- points.on.land(x=samples.df.temp)
		} else {
			samples.temp           <- samples.df.temp
			# colnames(samples.temp) <- c("longitude","latitude")
		}
		# Number of points sampled for each group after filtering points
	#	n.grp.pass  <- table(samples.temp[,"group"])
		n.grp.pass  <- sapply(X=1:n.grp, FUN=function(x) length(grep(x,samples.temp[,"group"])))
		#n.grp.pass  <- sapply(samples.temp,nrow)
		#sapply(n.grp.pass,funct)
		if(all(n.grp.pass >= grp.sizes)){
			PASS=TRUE
		} else {
			PASS=FALSE
		}
		# result.temp <- rbind(result.temp,samples.temp)
	}
	# For each group, keep corresponding 'grp.sizes' number of samples
	result <- samples.temp[unlist(lapply(1:n.grp,function(x) grep(x,samples.temp[,"group"])[1:grp.sizes[x]])),]
	
	#result0 <- lapply(1:n.grp,FUN=function(x){ samples.temp[[x]][sample(x=1:nrow(samples.temp[[x]]), size=grp.sizes[x], replace=FALSE),]}) # samples.temp[[x]]
	#result  <- do.call(rbind,result0)
	# Make "group" column a character variable
	mode(result[,"group"]) <- "character"
	# Convert data frame to spatial points object.
	result.sp <- sp::SpatialPoints(result[,c(1,2)])
	# Minimum convex hull polygon of all points.
	mcp.sp         <- suppressWarnings(adehabitatHR::mcp(result.sp,percent=100))
	#mcp.sp.area <- attributes(attributes(attributes(mcp.sp)$"polygons"[[1]])[[1]][[1]])$area
	mcp.sp.area    <- mcp.sp@polygons[[1]]@Polygons[[1]]@area
	absolute.areas <- c(windowarea=wnd.area,regionsize.dd0=regionsize.dd0,regionsize.dd=regionsize.dd,min.convex.hull.result=mcp.sp.area)
	relative.areas <- round((absolute.areas/wnd.area),digits=4)
	areas.df       <- data.frame(absolute =absolute.areas, relative.to.wnd=relative.areas)
	#return(show.plot)
	if(show.plot){
		world.land     <- rnaturalearth::ne_countries(scale=110, returnclass="sf")
		world.land10   <- rnaturalearth::ne_countries(scale=10, returnclass="sf")
		result.df2plot <- result
		xdist <- geodist::geodist(result.df2plot[order(result.df2plot[,1]),][c(1,nrow(result.df2plot)),1:2],measure="geodesic",sequential=T)
		ydist <- geodist::geodist(result.df2plot[order(result.df2plot[,2]),][c(1,nrow(result.df2plot)),1:2],measure="geodesic",sequential=T)
		plot.margin <- which(c(xdist,ydist)==max(c(xdist,ydist)))
		zoommaparea.km <- (xdist/1000)*(ydist/1000)
		if(zoommaparea.km <10e6){
			zoom.map <- world.land10
		} else {
			zoom.map <- world.land
		}
		#####
		spbb   <- sp::bbox(result.sp)
		spbb2  <- apply(spbb,MARGIN=1,FUN=rangeBuffer,f=0.1)
		dfbb   <- data.frame(x1=spbb[1,1], x2=spbb[1,2], y1=spbb[2,1], y2=spbb[2,2])
		ggrect <- ggplot2::geom_rect(data=dfbb, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),fill="black",color="black",alpha=0.1)
		zoom.plot   <- ggplot2::ggplot(data = zoom.map) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = result.df2plot, ggplot2::aes(x = X, y = Y, fill=group), size = 2, shape = 21) + ggplot2::coord_sf(xlim =range(result.df2plot[,1]),ylim = range(result.df2plot[,2])) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::theme(legend.position = "none") + ggplot2::xlab("longitude") + ggplot2::ylab("latitude")
		global.plot <- ggplot2::ggplot(data = world.land) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = result.df2plot, ggplot2::aes(x = X, y = Y, fill=group), size = 2, shape = 21) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank()) + ggrect
		###### Thus far failed attempt to reproject points and basemap using a different CRS
		grobs.list <- list(ggplot2::ggplotGrob(global.plot),ggplot2::ggplotGrob(zoom.plot))
		if(plot.margin==2){
			bothmaps   <- gridExtra::arrangeGrob(grobs=grobs.list, layout_matrix=matrix(c(1,1,2),ncol=3))
		} else {
			bothmaps   <- gridExtra::arrangeGrob(grobs=grobs.list, layout_matrix=matrix(c(1,2),ncol=1))
		}
		#if(show.plot){
		grid::grid.newpage()
		grid::grid.draw(bothmaps)
		#}
		if(FALSE){
			# 'result' data frame converted to SpatialPointsDataFrame
			# result.spdf <- sp::SpatialPointsDataFrame(coords=result[,c(1,2)],data=data.frame(group=result[,3]))
			### Set the 'original' CRS of the points
			# sp::proj4string(result.spdf) <- sp::CRS("EPSG:4326")
		#	raster::crs(result.spdf) <- sp::CRS("EPSG:4326")
			# Name of CRS to use
			crs.name <- "ESRI:54009"
			# Transformed coordinates
			result_proj <- dfTransform(df=result,CRS2=crs.name)
			# Reprojecting the 110m resolution map from latlon to Moller
			sf.world110_proj <- sf::st_transform(x=world.land,crs=sf::st_crs(crs.name))
			### ggplot with Moller project 110m resolution basemap
			gg.world110_proj <- ggplot2::ggplot(data = sf.world110_proj) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())
			### Adding the reprojected points to the ggplot
			gg.map <- gg.world110_proj +  ggplot2::geom_point(data = result_proj, ggplot2::aes(x = X, y = Y, fill=group), size = 2, shape = 21) # + ggplot2::coord_sf(crs= sf::st_crs(crs.name), datum = sf::st_crs(crs.name) )
		}
	} else {
		bothmaps <- NULL
	}
	#result <- result.temp[sample(x=1:nrow(result.temp),size=size),]
	result.df <- result
	if(show.plot){
		return(list(coords.df=result.df,map.ggplot=bothmaps))
	} else {
		return(result.df)
	}
	#if(return.as=="matrix"){
	#	result <- data.matrix(result)
	#}
	#if(return.as=="SP"){
	#	result.df <- result
	#	result <- list(result.sp=result.sp,sample.area.sp=sample.area.sp,result.df=result.df,zoom.gg=zoom.plot,global.gg=global.plot,bothmaps.gg=bothmaps,mcp.hull=mcp.sp,areas.df=areas.df,group.center.distances= grp.center.dists,grp.scaler.info=c(grp.scaler.start,grp.scaler),counters=c(ctr,ctr2))
	#	#result <- sp::SpatialPoints(result)
	#}
	#if(return.as=="DF_plot"){
	#	result.df <- result
	#	result <- list(coords.df=result.df,map.ggplot=bothmaps)
	#}
	#result
}
#' @examples
#' library(misc.wrappers)
#' # Sample 50 points from 10-degree radius area with center located on land somewhere between -50 and 50 degrees latitude.
#' coords50 <- rcoords(r=10,size=50,limits=c(-180,180,-50,50))
#' 
#' # Sample 100 points each from two allopatric groups (populations), both with centers somewhere on land between -50 and 50 latitude and centers are between 5 and 20 degrees from each other. 
#' coords50.K2.allopatric <- rcoords(r=5,size=c(50,100),limits=c(-180,180,-50,50),n.grp=2,interactions=c(0,0))
#'
#' # Use default settings to generate 100 points over land for one group
#' coords100 <- rcoords(show.plot=T)
#' 

#' @title Buffered range
#' 
#' Calculate the limits of a buffered range for a numerical vector.
#' 
#' @param x Numerical vector
#' @param f Number greater than zero indicating the amount of expansion. Default 0.1. Zero is the same as range(x).
#' @return Numerical vector with minimum and maximum limits of the buffered range.
#' @export rangeBuffer
rangeBuffer <- function(x,f){
	Rx <- range(x)
	Dx <- diff(Rx)
	Bx <- ((Dx*f)/2)
	res <- c((Rx[1]-Bx),(Rx[2]+Bx))
	res
}

#' @title Random point conditioned on distance
#' 
#' Return x and y coordinates of random points conditional on a distance or range of distances from input point.
#' 
#' @param p1 Numerical vector with decimal longitude and latitude of first point.
#' @param d Number specifying the distance (in decimal degrees) between the new point and the point with coordinates 'x1y1', or, a numerical vector, in which case the distance from x1y1 will be randomly drawn from values in the range of d.
#' @param n Number of points to return. Default 1.
#' @return Numerical vector with longitude and latitude of new point.
#' @export rx2y2
rx2y2 <- function(p1,d,n=1){
	# direction (angle) of new point(s) relative to p1
	rtheta <- sample(seq(from=0,to=(2*pi),length.out=(n*100)),size=n,replace=FALSE)
	# distance(s) from p1
	rd     <- sample(seq(from=min(d),to=max(d),length.out=(n*100)),size=n,replace=FALSE)
	xn <- (rd*cos(rtheta)) + p1[1]
	yn <- (rd*sin(rtheta)) + p1[2]
	data.frame(x=xn,y=yn)
}


#' @title corner coordinates of area described by bounding box
#' 
#' Returns a 2x4 matrix with coordinates x,y of the four corner points for a region with an extent determined by the input bounding box.
#' 
#' @param bb 2x2 numerical matrix (bounding box matrix) with x,y rows and and min, max columns
#' @return 2x4 matrix with coordinates x y of points at the corners the region with extent defined by 'bb'
#' @export bbox2points
bbox2points <- function(bb){
	res0 <- do.call(rbind,list(bb[,1],bb[,2] ,diag(bb), diag(apply(bb,1,rev))))
	# returns points clockwise from lower left
	res0[order(order(res0[,1]),order(res0[,2],decreasing=T)),]
}


# For an input function 'infunc', the function 'assignformals2global' gets the argument names and default values of 'infunc', and assigns the values of those arguments to objects of the same name in the global environment. Useful when debugging functions.
# Doesnt really work as a function because the objects have class 'call'
#assignformals2global <- function(infunc){
#	# infunc <- rcoords
#	forms <- formals(infunc)
#	for(i in 1:length(forms)){
#		assign(names(forms)[i],forms[[i]],envir=globalenv(),inherits=TRUE)
#	}
#}

#' @title Transform CRS coordinates data frame
#' 
#' Transform coordinates held in a dataframe given the starting and final coordinate reference systems
#' 
#' @param data frame with, minimally, longitude and latitude columns
#' @param x Name of the column (variable) with longitude
#' @param y Name of the column (variable) with latitude
#' @param CRS1 CRS of the input coordinates
#' @param CRS2 CRS to use for the output coordinates
#' @return data frame with reprojected coordinates. Other columns of the input data frame are preserved.
#' @export dfTransform
dfTransform <- function(df,x,y,CRS1="EPSG:4326",CRS2){
	if(missing(x)){
		x     <- 1
		xvals <- df[,1]
	} else {
		x <- grep(x,names(df))
	}
	if(missing(y)){
		y <- 2
	} else {
		y <- grep(y,names(df))
	}
	if(missing(CRS2)){
		stop("CRS2 must be supplied")
	}
	xy.df <- data.frame(X=df[,x],Y=df[,y])
	if(ncol(df)>2){
		data.df <- data.frame(df[, setdiff(1:ncol(df),c(x,y))])
	} else {
		data.df <- NULL
	}
	### Convert to Spatial Points object
	xy.sp <- sp::SpatialPoints(xy.df)
	### Set the 'original' CRS
	raster::crs(xy.sp) <- sp::CRS(CRS1)
	### Transform CRS
	xyT.sp <-  sp::spTransform(xy.sp,sp::CRS(CRS2))
	### Extract coordinates
	xyT.df <- sp::coordinates(xyT.sp)
	if(!is.null(data.df)){
		res.df           <- cbind(xyT.df,data.df)
		colnames(res.df) <- c(names(df)[c(x,y)], names(df)[-c(x,y)])
	} else {
		res.df <- xyT.df
	}
	res.df
}

#' @title ggplot projected Natural Earth
#' 
#' This function creates a ggplot object of Earth from the 110m resolution Natural Earth dataset, and the function is useful for exploring different map projections.
#' 
#' @param CRS Character string with name of EPSG or ESRI coordinate reference system ID. Default is 'EPSG:4326'.
#' @param df Optional data frame with longitude and latitude coordinates (datum must be WGS84) in first two columns, respectively, for points that should be plotted on the map.
#' @return ggplot object of Natural Earth data plotted using the projection specified by 'CRS' argument.
#' @export gg.world110_proj
gg.world110_proj <- function(CRS="EPSG:4326",df=NULL){
	world.land     <- rnaturalearth::ne_countries(scale=110, returnclass="sf")
	# world.land10   <- rnaturalearth::ne_countries(scale=10, returnclass="sf")
	crs.name <- CRS
	# Transformed coordinates
	if(!is.null(df)){
		df_proj <- dfTransform(df=df,CRS2=crs.name)
		names(df_proj[,1:2]) <- c("X","Y")
	}
	# Reprojecting the 110m resolution map from latlon to Moller
	sf.world110_proj <- sf::st_transform(x=world.land,crs=sf::st_crs(crs.name))
	### ggplot of projected basemap
	gg.map <- ggplot2::ggplot(data = sf.world110_proj) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1)) + ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())
	if(!is.null(df)){
		### Adding the reprojected points to the ggplot
		gg.map <- gg.map +  ggplot2::geom_point(data = df_proj, ggplot2::aes(x = X, y = Y), fill="black",size = 2, shape = 21)
	}
	gg.map
}


#' Evaluate function name with lists of arguments
#' 
#' This function determines which entries in the lists 'supplied' and 'passed' are arguments of the function 'fun.name', and then evaluates the function with those arguments.
#' Entries with names that are not arguments of 'fun.name' are ignored, making it possible to include arguments for multiple functions called in a parent function with argument '...'. See examples.
#' 
#' 
#' @param fun.name Character string with name of the function to be evaluated
#' @param supplied.list Arguments and values of fun.name defined explicitely in a list
#' @param passed A list with names equal to argument names in the set of names of some arguments with values to use for arguments. This list is typically formed by capturing the "additional arguments to be passed with '...'" in a parent function list(...).
#' @param ... Arguments and values of fun.name defined explicitely and not together in a list
#' @return Value returned by evaluating the function with 'fun.name' and arguments in lists 'supplied' and 'passed'
#' @export evalfun
evalfun <- function(fun.name, supplied.list=NULL, passed=list(...), ...){
	# List with default values for all arguments of the function
	fun.args0 <- formals(get("fun.name"))
	# Keep arguments except for "..."
	fun.args1 <- fun.args0[setdiff(names(fun.args0),"...")]
	# Character vector with names of all arguments, including those without default values, but not "..."
	argnames0 <- names(fun.args1)
	# Character vector with names of arguments with defined values
	argnames <- argnames0[!sapply(fun.args1,is,"name")]
	# subset of fun.args with defined values
	fun.args2 <- fun.args1[!sapply(fun.args1,is,"name")]
	# List of arguments supplied direct
	suppllied.direct <- list(...)
	# If arguments in the supplied list were also supplied directly, the ignored/remove the list-supplied version.
	if(any((names(supplied.list) %in% names(suppllied.direct)))){
		supplied.list <- supplied.list[!(names(supplied.list) %in% names(suppllied.direct))]
	}
	# List of the arguments suppplied explicitely, either directly or as a list.
	supplied <- c(supplied.list,suppllied.direct)
	# If passed arguments are already in supplied, the supplied version has priority and the duplicated version in 'passed' is ignored/removed
	if(any((names(passed) %in% names(supplied)))){
		passed <- passed[!(names(passed) %in% names(supplied))]
	}
	# Merging 'supplied' and 'passed'.
	argslist <- c(supplied, passed)
	if(!is.null(argslist)){
		# Check if any of the arguments in argslist are for this function
		if(any( names(argslist) %in% argnames0)){
			args.use <- c(argslist[names(argslist) %in% argnames0], fun.args1[!(argnames0 %in% names(argslist))])
		} else{
			# Not sure which one of these would work
			# args.use <- fun.args2
			args.use <- fun.args1
		}
	} else {
		args.use <- fun.args1
	}
	args.expr=paste(paste0(names(args.use),"=", unname(c(args.use))), collapse=",")
	form.expr=paste0(fun.name,"(",args.expr,")")
	#return(list(supplied.direct=suppllied.direct,supplied.list=supplied.list,supplied=supplied,passed=passed,argslist=argslist,args.use=args.use,form.expr=form.expr))
	# parse the expression
	eval(parse(text= form.expr))
} 
#' @examples
#' # Consider two functions, 'testfunA' and 'testfunB', each with two arguments and no shared argument names):
#' testfunA <- function(argA1=10,argA2=10){
#' 	return(argA1*argA2)
#' }
#' 
#' testfunB <- function(argB1=10,argB2=10){
#' 	return(argB1/argB2)
#' }
#' 
#' # 'evalfun' can be used to evaluate 'testfunA' and 'testfunB' from within another function, by passing arguments to 'evalfun'. The parent function must not have any argument names shared with the daughter functions that are not but not intended for inheritance.
#' testfunAB <- function(...){
#'     #moreArgs <- list(...)
#'     result1 <- evalfun("testfunA",list(...))
#'     result2 <- evalfun("testfunB",list(...))
#'     list("result1"=result1,"result"=result2)
#' }
#' # Result
#' testfun(argA1=5, argA2=6, argB1=7, argB2=8)


