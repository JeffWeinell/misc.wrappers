#' @title plot_eems function
#' 
#' Invokes the reemsplots2 function make_eems_plots, with additional functionality such as adding island borders and points showing individual coordinates..
#' 
#' @param mcmcdir Character string to directory of a particular mcmc chain
#' @param coords Optional character string to sample coordinates file.
#' @param output.plot.path Character string with path where to save output. If NULL, the plot is printed but not saved.
#' @return NULL; plots are printed to screen and saved to output.plot.path
#' @export plot_eems
plot_eems <- function(mcmcdir,coords=NULL,output.plot.path=NULL,plot.geography=T,plot.oceans=T){
	spdf_world_10    <- rnaturalearth::ne_countries(scale=10)
	#spdf_oceans_10   <- invisible(rnaturalearth::ne_download(scale = 10, type = 'ocean', category = 'physical'))
	# oceans <- misc.wrappers::spdf_oceans_10
	spdf_oceans_10 <- misc.wrappers::spdf_oceans_10
	### Define CRS
	crs.string <- "+init=EPSG:4326"
	### Generate EEMS plots
	plots.chain     <- suppressWarnings(reemsplots2::make_eems_plots(mcmcdir, longlat = TRUE))
	### extract the first eems plot (this is a ggplot)
	qrates <- as.numeric(unlist(strsplit(paste0(readLines(paste0(mcmcdir,"/mcmcqrates.txt")),collapse=""),split=" ")))
	mrates <- as.numeric(unlist(strsplit(paste0(readLines(paste0(mcmcdir,"/mcmcmrates.txt")),collapse=""),split=" ")))
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
	### Two-column matrix of lat and lon coordinates of 'eems tiles'. These locations are same for all four maps, and therefore they are extracted from the first map.
	#data.xy <- eems.data[[1]][,c("x","y")]
	### Four-column character matrix containing the Hex Color codes at each Lat/Lon for each of the four maps.
	#fill.hex <- do.call(cbind,lapply(X=eems.data[1:4],FUN=function(df){df[,c("fill")]}))
	#fill.hsv <- t(matrix(rgb2hsv(col2rgb(fill.hex)),ncol=nrow(data.xy),byrow=T))
	#rownames.use <- c(paste0("h_",c("mrates01","mrates02","qrates01","qrates02")),paste0("s_",c("mrates01","mrates02","qrates01","qrates02")),paste0("v_",c("mrates01","mrates02","qrates01","qrates02")))[c(1,5,9,2,6,10,3,7,11,4,8,12)]
	#colnames(fill.hsv) <- c(paste0("h_",c("mrates01","mrates02","qrates01","qrates02")),paste0("s_",c("mrates01","mrates02","qrates01","qrates02")),paste0("v_",c("mrates01","mrates02","qrates01","qrates02")))[c(1,5,9,2,6,10,3,7,11,4,8,12)]
	### Data frame with x and y coordinates and fill color (hsv) for each of the four maps.
	#data.xy.hsv <- cbind(data.xy,fill.hsv)
	#fill.hsv.mrates01 <- cbind(fill.hex[,1],fill.hsv[,1:3])
	#fill.hsv.mrates02 <- cbind(fill.hex[,2],fill.hsv[,4:6])
	#fill.hsv.qrates01 <- cbind(fill.hex[,3],fill.hsv[,7:9])
	#fill.hsv.qrates02 <- cbind(fill.hex[,4],fill.hsv[,10:12])
	#### Colors to use for the legend
	#fill.hsv.mrates01.ordered <- unique(fill.hsv.mrates01[order(fill.hsv.mrates01[,2], fill.hsv.mrates01[,3],fill.hsv.mrates01[,4]),])
	#fill.hsv.mrates02.ordered <- unique(fill.hsv.mrates02[order(fill.hsv.mrates02[,2], fill.hsv.mrates02[,3],fill.hsv.mrates02[,4]),])
	#fill.hsv.qrates01.ordered <- unique(fill.hsv.qrates01[order(fill.hsv.qrates01[,2], fill.hsv.qrates01[,3],fill.hsv.qrates01[,4]),])
	#fill.hsv.qrates02.ordered <- unique(fill.hsv.qrates02[order(fill.hsv.qrates02[,2], fill.hsv.qrates02[,3],fill.hsv.qrates02[,4]),])
	#eems1.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[1]][,c("x","y")],data=eems.df.data[[1]])
#	eems2.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[2]][,c("x","y")],data=eems.df.data[[2]])
#	eems3.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[3]][,c("x","y")],data=eems.df.data[[3]])
#	eems4.data.spdf <- sp::SpatialPointsDataFrame(coords=eems.data[[4]][,c("x","y")],data=eems.df.data[[4]])
	eems.data.spdf <- lapply(X=1:4,FUN=function(i){sp::SpatialPointsDataFrame(eems.data[[i]][,c("x","y")],data=eems.df.data[[i]],proj4string=CRS(crs.string))})
	### Setting the crs
	
	### Coordinates where legend will go.
#	y.mrates01.legend <- seq(from=attributes(eems.data.spdf[[1]])$bbox["y","min"],to=attributes(eems.data.spdf[[1]])$bbox["y","max"],length.out=nrow(fill.hsv.mrates01.ordered))
#	x.mrates01.legend <- rep((attributes(eems.data.spdf[[1]])$bbox["x","max"]*1.05),nrow(fill.hsv.mrates01.ordered))
#	legend("topright", inset=c(-0.1,0), fill=fill.hsv.mrates01.ordered[sort(samples),1], legend= rep("",100), title="")


#	lapply(X=1:4,FUN=function(x){raster::crs(eems.data.spdf[[x]])   <- crs.string})
	# suppressWarnings(raster::crs(eems1.data.spdf)   <- crs.string)
	suppressWarnings(raster::crs(spdf_world_10)     <- crs.string)
	suppressWarnings(raster::crs(spdf_oceans_10)    <- crs.string)
	intersect.eems              <- lapply(X=1:4,function(x){sp::over(eems.data.spdf[[x]],spdf_world_10)})
	### interections of each plot with land
	# eems1.intersect.spdf <- intersect.eems[[1]]
	# eems2.intersect.spdf <- intersect.eems[[2]]
	# eems3.intersect.spdf <- intersect.eems[[3]]
	# eems4.intersect.spdf <- intersect.eems[[4]]
	### Drop NA rows of each plot
	# eems1.intersect.spdf <- eems.data.spdf[[1]][apply(X=intersect.eems[[1]],MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	# eems2.intersect.spdf <- eems.data.spdf[[2]][apply(X=intersect.eems[[2]],MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	# eems3.intersect.spdf <- eems.data.spdf[[3]][apply(X=intersect.eems[[3]],MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	# eems4.intersect.spdf <- eems.data.spdf[[4]][apply(X=intersect.eems[[4]],MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	# ### All four data frames held in a list
	# eems.intersect.spdf <- list(eems1.intersect.spdf,eems2.intersect.spdf,eems3.intersect.spdf,eems4.intersect.spdf)
	# intersect.eems.unique       <- eems1.data.spdf[apply(X=intersect.eems,MARGIN=1,FUN=function(x){any(!is.na(x))}),]
	### Each plot
	eems.intersect.spdf <- lapply(X=1:4,FUN=function(i){eems.data.spdf[[i]][apply(X=intersect.eems[[i]],MARGIN=1,FUN=function(x){any(!is.na(x))}),]})
	#if(!is.null(output.plot.path)){
	#	pdf(output.plot.path)
	#}
	layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
	par(mar=c(2,2,1,1))
	for(i in 1:4){
		intersect.eems.unique.plot <-  sp::plot(eems.intersect.spdf[[i]],col=eems.intersect.spdf[[i]]$fill)
		### Plot a legend. This doesnt quite work yet because the colors are not sorted properly.
		# intersect.eems1.unique.plot  <- sp::plot(eems1.intersect.spdf.unique,col=eems1.intersect.spdf.unique$fill)
		# intersect.eems2.unique.plot  <- sp::plot(eems2.intersect.spdf.unique,col=eems2.intersect.spdf.unique$fill)
		# intersect.eems3.unique.plot  <- sp::plot(eems3.intersect.spdf.unique,col=eems3.intersect.spdf.unique$fill)
		# intersect.eems4.unique.plot  <- sp::plot(eems4.intersect.spdf.unique,col=eems4.intersect.spdf.unique$fill)
		if(plot.geography){
			ne_plot_10                 <- sp::plot(spdf_world_10,add=T,border="black",lwd=1.5)
		}
		if(plot.oceans){
			ne_oceans.plot_10          <- sp::plot(spdf_oceans_10,add=T,border="black",lwd=1.5,col="white")
		}
		if(!is.null(coords)){
			coordinates <- data.matrix(read.table(coords,header=F))
			plot(coordinates,add=T,pch=20)
		}
		maps::map.axes()
		legend_image <- as.raster(matrix(rev(c("#994000","#CC5800","#FF8F33","#FFAD66","#FFCA99","#FFE6CC","#FBFBFB","#CCFDFF","#99F8FF","#66F0FF","#33E4FF","#00AACC","#007A99")), ncol=1))
		log.minm <- min(log10(mrates))
		log.maxm <- max(log10(mrates))
		log.minq <- min(log10(qrates))
		log.maxq <- max(log10(qrates))
		if(i==1){
			plot(c(0,2),c((log.minm*1.1),(log.maxm*1.4)),type = 'n', axes = F,xlab = "", ylab = "",main="")
			rasterImage(legend_image, 0, log.minm, 1,log.maxm)
			text(x=1.5, y = c(-2,-1,0,1,2), labels = c("-2","-1","0","1","2"),cex=1)
			text(x=-0.2,y=log.maxm+0.2,labels="log(m)",pos=4,cex=1.25)
			segments(x0=rep(c(0,0.75),5),y0=c(-2,-2,-1,-1,0,0,1,1,2,2),x1=rep(c(0.25,1),5),col="white",lwd=1.5)
		}
		if(i==2){
			plot(c(0,5),c(-1.5,1.5),type = 'n', axes = F,xlab = "", ylab = "",main="")
			rasterImage(legend_image, 0, -1, 1, 1)
			text(x=1.1, y = c(-0.9,0.9), labels = c("P{log(m) < 0} = 0.9","P{log(m) > 0} = 0.9"),pos=4,cex=0.8)
			segments(x0=rep(c(0,0.75),2),y0=c(-0.9,-0.9,0.9,0.9),x1=rep(c(0.25,1),2),col="white",lwd=1.5)
		}
		if(i==3){
			plot(c(0,2),c((log.minq*1.1),(log.maxq*1.4)),type = 'n', axes = F,xlab = "", ylab = "",main="")
			rasterImage(legend_image, 0, log.minq, 1,log.maxq)
			text(x=1.5, y = c(-0.10,-0.05,0,0.05,0.10), labels = c("-0.10","-0.05","0.00","0.05","0.10"))
			text(x=0,y=(log.maxq*1.2),labels="log(q)",pos=4)
		}
		if(i==4){
			plot(c(0,5),c(-2,2),type = 'n', axes = F,xlab = "", ylab = "",main="")
			rasterImage(legend_image, 0, -1, 1,1)
			text(x=1.1, y = c(-0.9,0.9), labels = c("P{log(q) < 0} = 0.9","P{log(q) > 0} = 0.9"),cex=0.8)
		}
	}
	#if(is.null(output.plot.path)){
	#	dev.off()
	#}
}




