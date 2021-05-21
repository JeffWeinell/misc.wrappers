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
#' @param regionsize Number specifying the approximate longitudinal diameter spanning all sampling regions. I have not verified the accuracy, but larger values generally yield larger extent across all samples.
#' @param samplesize Total number of points to return, across all groups.
#' @param grp.n.weights Probability weights for determing sample sizes of groups. Default weights are uniformly distributed.
#' @param grp.area.weights Probability weights for determing group region sizes. Default weights are uniformly distributed.
#' @param wnd Numerical vector specifying longitudinal and latitude limits for possible sampling. Default is to include all of Earth c(-180,180,-90,90).
#' @param over.land Logical indicating whether or not all returned points must occur over land. Default TRUE. Note that this condition is only applied to samples returned as output. Therefore, 'samplesize' is really a 'sample until' rule. This is faster than performing clipping operations on proposed sample polygons to conform to geography.
#' @param n.grp Number of spatial sampling groups.
#' @param interactions Numeric vector with the maximum and minimum amount of overlap between pairs of groups, calculated as (intersect area)/(mean group areas non-intersected). The default c(0,1) allows for all possible scenarios. Examples: c(0,0) specifies that groups must be allopatric; c(1,1) requires complete overlap of groups, which is not realistic given the stochasticity determining region sizes; c(0.5,1) requires that at least half-overlaps between groups; c(0.2,0.25) specifies a small contact zone.
##' @param d.grp Number controlling the distance between the centers of a pair of areas, as a function of the pair's radii. Default 1, which would allow sample regions to nearly coincide. Future option may allow for a pairwise distance matrix.
#' @param expf Number that affects dispersion relative to regionsize. Higher numbers increase dispersion. Default 8.
#' @param return.as Character string with class to use for object returned. Default "data.frame". Can also be "matrix" or "SP" (SpatialPoints).
#' @param show.plot Whether or not the points should be plotted on a low-resolution land map. The map is used is the rnaturalearth countries map, 110 meter resolution.
#' @return An object with class equal to the value of 'return.as' and containing the set of points that meet the specified sampling requirements. If return.as='matrix' or 'data.frame', the columns are 'X' (for longitude), 'Y' (for latitude), and 'group' (all 1 if 'n.grp'=1).
#' @export rcoords
rcoords <- function(regionsize, samplesize, n.grp=1, grp.n.weights=rep(1,n.grp), grp.area.weights=rep(1,n.grp), wnd=c(-180,180,-90,90), over.land=TRUE, interactions=c(0,1), expf=8, show.plot=FALSE, return.as="data.frame"){
	initial.list <- list(regionsize,grp.area.weights,grp.areas0)
	#result.temp        <- data.frame(NULL)
	PASS=FALSE
	miter <- 100
	ctr   <- 0
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
	while(!PASS){
		ctr <- ctr + 1
		if(ctr > miter){
			stop(paste("after",miter,"attemps, failed to find points that pass conditions. ctr=",ctr,"ctr2=",ctr2))
		}
		### Clear the data frame if too few points on land during previous attempt
		#result.temp         <- data.frame(NULL)
		initial <- FALSE
		### Find group centers that pass conditions
		ctr2 <- 0
		while(!initial){
			ctr2 <- ctr2 + 1
			if(ctr2 > miter){
				stop(paste("failed to initialize after",miter,"attemps. ctr=",ctr,"ctr2=",ctr2))
			}
			sample.center  <- c(x=sample(seq(wnd[1],to=wnd[2],by=0.01),size=1), y=sample(seq(wnd[3],to=wnd[4],by=0.01),size=1))
			if(n.grp>1){
				# Sampling distances ('d.grps') from group centers to the sample.center. These values are the main determinant of how groups are geographically arranged.
				# Value that may generate more space between groups.
			#	expf    <- 4
				d.grpsA <- rnorm(n=n.grp, mean=(regionsize/(n.grp*expf)), sd= expf)
				d.grpsB <- runif(n=n.grp, min=(regionsize/(n.grp*10)), max=(regionsize/(n.grp*5)))
				d.grps  <- apply(rbind(d.grpsA, d.grpsB), MARGIN=2, max)
				group.centers <- rx2y2(p1=sample.center,d=d.grps,n=n.grp)
				# Geographic sizes of each group.
				grp.areas0   <- regionsize*(grp.area.weights/(sum(grp.area.weights)))
				grp.areas    <- sapply(X=1:n.grp,FUN=function(x){rnorm(n=1,mean=grp.areas0[x],sd=0.5)})
				# d.grps <- (max(r)*(2*d.grp)):(max(r)*(3*d.grp))
			} else {
				group.centers <- sample.center
				grp.areas     <- rnorm(n=1,mean=regionsize,sd=1)
			}
			# Check that all group centers occur over land.
			if(over.land){
				if(length(points.on.land(x=group.centers)[,1]) < n.grp){
					next
				}
			}
			initial <- TRUE
		}
		# 
		#
		# return(grp.areas)
#		### Sample a point on earth to use as the center of the sampling circle (if n.grp=1), or group center if n.grp>1
#		sample.center  <- c(x=sample(seq(limits[1],to=limits[2],by=0.01),size=1), y=sample(seq(limits[3],to=limits[4],by=0.01),size=1))
#		if(n.grp>1){
#			d.grps <- (max(r)*2):(max(r)*3)
#			group.centers <- rx2y2(p1=sample.center,d=d.grps,n=n.grp)
#		} else {
#			group.centers <- sample.center
#		}
#		if(over.land){
#			if(length(points.on.land(x=group.centers)) < n.grp){
#				next
#			}
#		}
		### For each group, the number of samples that will be drawn (or saved, if over.land=TRUE) from the group's sample region.
		if(n.grp==1){
			grp.sizes <- samplesize
		} else {
			grp.sizes  <- c(table(sample(x=n.grp,size=samplesize,prob=grp.n.weights,replace=TRUE)))
		}
		### Creates a SpatialPolygons object (sampling region) for each group
		sample.area.sp <- lapply(1:n.grp, FUN=function(x){sampSurf::spCircle(radius= sqrt((grp.areas/pi))[x] , centerPoint=c(x=group.centers[x,1],y=group.centers[x,2]))[[1]]})
		### Check if groups meet conditions set by 'interactions' argument.
		if(!is.null(interactions) & n.grp>1){
			if(!(interactions[1]==0 & interactions[2]==1)){
				# Much less complex list holding the sample sampling polygons
				polygons.list <- lapply(X=1:n.grp, function(x){attributes(attributes(sample.area.sp[[x]])$"polygons"[[1]])[[1]][[1]]})
				# All possible pairwise group comparisons
				grp.pairs <- do.call(rbind,pset(1:n.grp,2,2))
				# Calculate intersect area yet and sum areas for each pair of groups. (Doesnt work yet).
				grp.int.area <- grp.sum.area <- grp.ratio <- c()
				for(i in 1:nrow(grp.pairs)){
					polygon.pair <- polygons.list[c(grp.pairs[i,1],grp.pairs[i,2])]
					pairsum.temp <- sum(c(attributes(polygon.pair[[1]])$area, attributes(polygon.pair[[2]])$area))
					#grp.int <- raster::intersect(sample.area.sp[[grp.pairs[i,1]]], sample.area.sp[[grp.pairs[i,2]]])
					grp.int <- rgeos::gIntersection(spgeom1=sample.area.sp[[grp.pairs[i,1]]],spgeom2=sample.area.sp[[grp.pairs[i,2]]]) # returns NULL if polygons do not intersect; returns SpatialPolygons object if they do intersect
					if(is.null(grp.int)){
						int.area.temp <- 0
					} else {
						int.area.temp <- attributes(attributes(attributes(grp.int)$"polygons"[[1]])[[1]][[1]])$area
					}
					ratio.temp <- (int.area.temp/pairsum.temp)
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
		### Convert each SpatialPolygons object to an SF object
	#	sample.area.sf <- sf::st_as_sf(sample.area.sp)
		sample.area.sf <- lapply(X=1:n.grp,FUN=function(x){sf::st_as_sf(sample.area.sp[[x]])})
		### Sample within each area. If points must be over land, sample four times as many points as requested by 'size' argument.
		if(over.land){
			nsamp <- (grp.sizes*4)
			samples.sf.temp  <- lapply(X=1:n.grp,FUN=function(x){sf::st_sample(x=sample.area.sf[[x]], size= nsamp[x])})
		} else {
			samples.sf.temp  <- lapply(X=1:n.grp,FUN=function(x) {sf::st_sample(x= sample.area.sf[[x]],size=grp.sizes[x])})
		}
		### Alternative sampling process
	#	samples2.sf.temp  <- sf::st_sample(x= sample.area.sf,kappa=10,mu=(size*4),scale=0.2,type="Thomas")      # spatstat.core::rThomas(kappa=10,scale=0.2, mu=(size*4),win=sample.area.sf)
		### Hold sampled coordinates as a list of data matrices
		samples.mat.temp <- lapply(X=1:n.grp,FUN=function(x){sf::st_coordinates(samples.sf.temp[[x]])})
		### Hold sampled coordinates as a list of data frames, and include a column in each data frame to indicate group assignment
		samples.df.temp <- lapply(X=1:n.grp,FUN=function(x){data.frame(X=samples.mat.temp[[x]][,1],Y=samples.mat.temp[[x]][,2],group=x)})
		### Plot samples on map (for debugging)
		# world.land <- rnaturalearth::ne_countries(scale=110,returnclass="sf")
		# ggplot2::ggplot(data = world.land) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = samples.df.temp, ggplot2::aes(x = X, y = Y), size = 2, shape = 21, fill = "darkred") # + coord_sf(xlim = c(-90, -78), ylim = c(24.5, 40), expand = FALSE)
		if(over.land){
			### Get subset of sampled coordinates that fall on land
			#samples.temp <- lapply(X=1:n.grp,FUN=function(x){points.on.land(x=samples.mat.temp[[x]])})
			samples.temp <- lapply(X=1:n.grp,FUN=function(x){points.on.land(x=samples.df.temp[[x]])})
			#samples.land.temp <- points.on.land(x=samples.df.temp)
		} else {
			samples.temp           <- samples.df.temp
			# colnames(samples.temp) <- c("longitude","latitude")
		}
		n.grp.pass  <- sapply(samples.temp,nrow)
		#sapply(n.grp.pass,funct)
		if(all(n.grp.pass >= grp.sizes)){
			PASS=TRUE
		} else {
			PASS=FALSE
		}
		# result.temp <- rbind(result.temp,samples.temp)
	}
	result0 <- lapply(1:n.grp,FUN=function(x){samples.temp[[x]][sample(x=1:nrow(samples.temp[[x]]), size=grp.sizes[x], replace=FALSE),]}) # samples.temp[[x]]
	result  <- do.call(rbind,result0)
	mode(result[,"group"]) <- "character"
	if(show.plot){
		world.land  <- rnaturalearth::ne_countries(scale=110,returnclass="sf")
		result.df2plot   <- result
		exf  <- 1.5
		xlim <- rangeBuffer(result.df2plot[,1],exf*1.5)
		ylim <- rangeBuffer(result.df2plot[,2],exf)
		#zoom.plot   <- ggplot2::ggplot(data = world.land) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = result.df2plot, ggplot2::aes(x = X, y = Y, fill=group), size = 2, shape = 21) + ggplot2::coord_sf(xlim = xlim, ylim = ylim ) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1))
		zoom.plot   <- ggplot2::ggplot(data = world.land) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = result.df2plot, ggplot2::aes(x = X, y = Y, fill=group), size = 2, shape = 21) + ggplot2::coord_sf(xlim =range(result.df2plot[,1]),ylim = range(result.df2plot[,2])) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1))
		global.plot <- ggplot2::ggplot(data = world.land) + ggplot2::geom_sf() + ggplot2::theme_classic() + ggplot2::geom_point(data = result.df2plot, ggplot2::aes(x = X, y = Y, fill=group), size = 2, shape = 21) + ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill=NA, size=1))
		#print(zoom.plot)
		grid::grid.newpage()
		grobs.list <- list(ggplot2::ggplotGrob(global.plot),ggplot2::ggplotGrob(zoom.plot))
		bothmaps <- gridExtra::arrangeGrob(grobs=grobs.list,layout=matrix(c(1,2),ncol=1))
		grid::grid.draw(bothmaps)
	}
	#result <- result.temp[sample(x=1:nrow(result.temp),size=size),]
	if(return.as=="matrix"){
		result <- data.matrix(result)
	}
	if(return.as=="SP"){
		result.df <- result
		result <- list(sp::SpatialPoints(result[,c(1,2)]),sample.area.sp,result.df,zoom.plot,global.plot)
		#result <- sp::SpatialPoints(result)
	}
	result
}
#' @examples
#' library(misc.wrappers)
#' # Sample 50 points from 10-degree radius area with center located on land somewhere between -50 and 50 degrees latitude.
#' coords50 <- rcoords(r=10,size=50,limits=c(-180,180,-50,50))
#' 
#' # Sample 100 points each from two allopatric groups (populations), both with centers somewhere on land between -50 and 50 latitude and centers are between 5 and 20 degrees from each other. 
#' coords50.K2.allopatric <- rcoords(r=5,size=c(50,100),limits=c(-180,180,-50,50),n.grp=2,interactions=c(0,0))

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
# rx2y2(x1y1=centerpoint,d=d.grp)


