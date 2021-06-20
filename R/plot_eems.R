#' @title plot_eems function
#' 
#' Invokes the reemsplots2 function make_eems_plots, with additional functionality such as adding island borders and points showing individual coordinates..
#' 
#' @param xdir Character string to directory of a particular mcmc chain, or, to a directory containing subdirectories "/mcmc/chain*".
#' @param plot.coords Whether or not to overlay the sample coordinates on EEMS maps. Default TRUE.
#' @param plot.geography Whether or not to overlay country borders on EEMS maps. Default TRUE.
#' @param mask.oceans Whether or not mask oceans on EEMS maps. Default TRUE.
#' @param include.out Character string vector indicating the type of output files to generate.
#' @param usechains Integer vector specifying which chains to make plots for; Default is NULL, in which case plots are made for all chains.
#' @return List of ggplots
#' @export plot_eems
plot_eems <- function(xdir, plot.coords=T, plot.geography=T, mask.oceans=T, include.out=c("pdf","raster"), usechains=NULL){
	### Default parameters of reemsplots2::make_eems_plots
	longlat = TRUE; dpi = 250; add_grid = FALSE; col_grid = "#BBBBBB"; add_demes = FALSE; col_demes = "#000000"; add_outline = FALSE; col_outline = "#FFFFFF"; eems_colors = NULL; prob_level = 0.9; m_colscale = NULL; q_colscale = NULL; add_abline = FALSE
	func_params <- list(add_grid = add_grid, add_demes = add_demes, add_outline = add_outline, eems_colors = eems_colors, prob_level = prob_level,col_grid = col_grid, col_demes = col_demes,col_outline = col_outline,m_colscale = m_colscale, q_colscale = q_colscale,add_abline = add_abline)
	# List of parameters defining how plots will be generated
	plot_params <- check_plot_params(pars=func_params)

	if("mcmc" %in% list.files(xdir)){
		nchains <- length(list.files(xdir,pattern="^runeems_snps_chain.*.sh$"))
		mcmcdirs <- file.path(xdir,"mcmc",paste0("chain",1:nchains))
	} else {
		nchains <- 1
		mcmcdirs <- xdir
		### overrides value of argument
		plot.coords <- FALSE
	}

	spdf_oceans_10 <- misc.wrappers::spdf_oceans_10
	oceans_sf      <- sf::st_as_sf(spdf_oceans_10)
	world_sf       <- rnaturalearth::ne_countries(scale=10,returnclass="sf")[1]
	if(is.null(usechains)){
		usechains <- 1:nchains
	} else {
		usechains <- intersect(usechains,1:nchains)
	}

	result <- list(); length(result) <- nchains
	if("raster" %in% include.out & length(usechains) > 1){
		bricklist <- list(); length(bricklist) <- nchains
	}
	for(ch in usechains){
		mcmcdir <- mcmcdirs[ch]
		save.in <- mcmcdir
		### Check that the necessary output files are present
		check_mcmcpath_contents(mcmcpath=mcmcdir)
		### Load some of the EEMs output
		dimns   <- read_dimns(mcmcpath=mcmcdir[1], longlat=longlat, nmrks = dpi)
		### Generate EEMS plots
		mplots     <- eems_contours(mcmcpath=mcmcdir, dimns=dimns, longlat=longlat, plot_params=plot_params, is_mrates = TRUE)
		qplots     <- eems_contours(mcmcpath=mcmcdir, dimns=dimns, longlat=longlat, plot_params=plot_params, is_mrates = FALSE)
		maps.chain <- list(mplots[[1]], mplots[[2]], qplots[[1]], qplots[[2]])
		# Check if more than two demes observed
		obs_demes <- read_matrix(file.path(mcmcdir[1], "rdistoDemes.txt"), ncol = 3)
		sizes     <- obs_demes[,3]
		include_disPlots <- (!sum(sizes > 1) < 2)
		if (include_disPlots) {
			dissimilarities      <- pairwise_dist(mcmcpath=mcmcdir, longlat=longlat, plot_params=plot_params)
			dissimilarity.plots  <- plot_pairwise_dissimilarities_(dissimilarities=dissimilarities, add_abline=add_abline)
			names(dissimilarity.plots) <- c("rdist01","rdist02","rdist03")
		} else {
			dissimilarity.plots <- list(NULL)
		}
		logP.plot <- list(pilogl01=plot_log_posterior(mcmcpath=mcmcdir))
		if(plot.coords){
			coordinates <- read.table(file.path(xdir,"data","data.coord"),header=F)[,c(1:2)]
			colnames(coordinates) <- c("Lon","Lat")
		}
		mapplot  <- list(); length(mapplot) <- 4
		for(i in 1:4){
			mapplot.i  <- maps.chain[[i]]
			coords <- unique(mapplot.i$data[,c(1:2)])
			colnames(coords) <- c("Lon","Lat")
			x.min    <- min((coords[,1]-0.5))
			x.max    <- max((coords[,1]+0.5))
			y.min    <- min((coords[,2]-0.5))
			y.max    <- max((coords[,2]+0.5))
			# Add country borders
			if(plot.geography){
				mapplot2.i <- suppressMessages(mapplot.i + ggplot2::geom_sf(data=world_sf, colour = "black", fill = NA, inherit.aes=FALSE)) # + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
			} else {
				mapplot2.i <- mapplot.i
			}
			# Mask oceans
			if(mask.oceans){
				mapplot3.i <- mapplot2.i + ggplot2::geom_sf(data=oceans_sf, colour = NA, fill = "white", inherit.aes=FALSE) # + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max),expand=FALSE)
			} else {
				mapplot3.i <- mapplot2.i
			}
			# Crop to window defined by xmin, xmax, ymin, ymax
			mapplot4.i <- mapplot3.i + ggplot2::coord_sf(xlim=c(x.min, x.max), ylim=c(y.min, ymax=y.max), expand=FALSE)
			# Add points for coordinates of samples
			if(plot.coords){
				mapplot[[i]] <- mapplot4.i + ggplot2::geom_point(data = coordinates, ggplot2::aes(x = Lon, y = Lat), size = 1, shape = 21, fill = "black")
			} else {
				mapplot[[i]] <- mapplot4.i
			}
		}
		names(mapplot) <- c("mrates01","mrates02","qrates01","qrates02")
		result.gglist <- c(mapplot, dissimilarity.plots,logP.plot)
		if("pdf" %in% include.out){
			save.as.pdf <- file.path(save.in,"EEMS_maps.pdf")
			pdf(height=6, width=8, file=save.as.pdf, onefile=TRUE)
				print(result.gglist)
			dev.off()
		}
		if("raster" %in% include.out){
			save.as.raster  <- file.path(save.in,"EEMS_maps.tif")
			plots.dat       <- lapply(X=maps.chain,FUN=function(x) {as.data.frame(x$data)})
			rasterlist <- list(); length(rasterlist) <- 4
			for(j in 1:4){
				dat.temp <- plots.dat[[j]]
				e.temp   <- raster::extent(dat.temp[,1:2])
				r.temp   <- raster::raster(e.temp, ncol=length(unique(dat.temp$x)), nrow=length(unique(dat.temp$y)))
				rasterlist[[j]] <- raster::rasterize(dat.temp[, 1:2], r.temp, dat.temp[,3], fun=mean, background=-9999)
			}
			# Hold the rasters in a rasterbrick
			res.brick <- raster::brick(rasterlist)
			# Write the brick as a multi-band geotiff
			wb <- suppressWarnings(raster::writeRaster(x=res.brick, filename=save.as.raster, format="GTiff", overwrite=T))
			if(length(usechains) > 1){
				bricklist[[ch]] <- res.brick
			}
		}
		result[[ch]] <- result.gglist
	}
	### For each of the four EEMs plots, generate a raster with chain means. Hold these in a brick and write as a GeoTiff.
	if("raster" %in% include.out & length(usechains) > 1){
		bricklist2 <- bricklist[usechains]
		bp1 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[1]]}))
		bp2 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[2]]}))
		bp3 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[3]]}))
		bp4 <- raster::brick(lapply(X=1:length(bricklist2), FUN=function(x) {bricklist2[[x]][[4]]}))
		rp1.mean <- raster::mean(bp1)
		rp2.mean <- raster::mean(bp2)
		rp3.mean <- raster::mean(bp3)
		rp4.mean <- raster::mean(bp4)
		bmeans   <- raster::brick(list(rp1.mean,rp2.mean,rp3.mean,rp4.mean))
		bw.mean  <- suppressWarnings(raster::writeRaster(x=bmeans, filename=file.path(xdir, "EEMS_maps_chainMeans.tif"), format="GTiff", overwrite=T))
	}
	result2 <- result[usechains]
	result2
}

### Defining some internal functions and parameters.
### This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').

#' @title Check MCMC Path Contents
#' 
#' Checks that all of the EEMs output files are present in the output directory
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @return NULL if all expected files generated by EEMs are present
#' @export check_mcmcpath_contents
check_mcmcpath_contents <- function(mcmcpath) {
	for (path in mcmcpath) {
		for (file in c("rdistJtDobsJ.txt", "rdistJtDhatJ.txt", "rdistoDemes.txt","mcmcmtiles.txt", "mcmcmrates.txt", "mcmcxcoord.txt","mcmcycoord.txt", "mcmcqtiles.txt", "mcmcqrates.txt","mcmcwcoord.txt", "mcmczcoord.txt", "mcmcpilogl.txt","outer.txt", "demes.txt", "edges.txt", "ipmap.txt","eemsrun.txt")) {
			if (!file.exists(file.path(path, file)))
				stop("Each EEMS output folder should include ", file)
		}
	}
}


#' @title Check Plot Parameters
#' 
#' Internal function used for making EEMS plots
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param pars List of paramaters
#' @return List of parameters for plotting
#' @export check_plot_params 
check_plot_params <- function(pars) {
	if (is.logical(pars$add_grid)) {
		pars$add_grid <- pars$add_grid[1]
	} else {
		pars$add_grid <- FALSE
	}
	
	if (is_color(pars$col_grid)) pars$col_grid <- pars$col_grid[1]
	else pars$col_grid <- "#BBBBBB"

	if (is.logical(pars$add_outline)) pars$add_outline <- pars$add_outline[1]
	else pars$add_outline <- FALSE
	if (is_color(pars$col_outline)) pars$col_outline <- pars$col_outline[1]
	else pars$col_outline <- "#EEEEEE"

	if (is.logical(pars$add_demes)) pars$add_demes <- pars$add_demes[1]
	else pars$add_demes <- FALSE
	if (is_color(pars$col_demes)) pars$col_demes <- pars$col_demes[1]
	else pars$col_demes <- "#000000"
	if (is.logical(pars$add_seeds)) pars$add_seeds <- pars$add_seeds[1]
	else pars$add_seeds <- TRUE

	### Color scale for eems plot 1. Default input pars$m_colscale is null (not numeric).
	if (is.numeric(pars$m_colscale)) {
		pars$m_colscale <- pars$m_colscale
	} else {
		pars$m_colscale <- c(-2.5, 2.5)
	}
	### Color scale for eems plot 3. Default input pars$q_colscale is null (not numeric).
	if (is.numeric(pars$q_colscale)) {
		pars$q_colscale <- pars$q_colscale
	} else {
		pars$q_colscale <- c(-0.1, 0.1)
	}
	### Colors to use for color scale
	if (length(pars$eems_colors) < 2 || any(!is_color(pars$eems_colors))){
		pars$eems_colors <- default_eems_colors()
	}
	if (is.null(pars$prob_level)) prob_level <- 0.9
	else prob_level <- pars$prob_level
	prob_level <- prob_level[prob_level > 0.5 & prob_level < 1]
	if (length(prob_level) != 1) prob_level <- 0.9
	pars$prob_level <- prob_level
	pars
}

#' @title Read dimns
#' 
#' Internal function used for making EEMS plots. 
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @param longlat Whether or not to use long lat coordinates. Default TRUE
#' @param nmrks Number of marks. Default 100.
#' @return A list with parameters to use for generating EEMS plots.
#' @export read_dimns
read_dimns <- function(mcmcpath, longlat=TRUE, nmrks = 100) {
	outer <- read_matrix(file.path(mcmcpath, "outer.txt"))
	ipmap <- read_vector(file.path(mcmcpath, "ipmap.txt"))
	demes <- read_matrix(file.path(mcmcpath, "demes.txt"))
	edges <- read_matrix(file.path(mcmcpath, "edges.txt"))
	dist_metric <- get_dist_metric(mcmcpath)
	if (!longlat) {
		outer <- outer[, c(2, 1)]
		demes <- demes[, c(2, 1)]
	}
	xlim <- range(outer[, 1])
	ylim <- range(outer[, 2])
	aspect <- (diff(ylim) / diff(xlim)) / cos(mean(ylim) * pi / 180)
	aspect <- abs(aspect)
	if (aspect > 1) {
		nxmrks <- nmrks
		nymrks <- round(nxmrks * aspect)
	} else {
		nymrks <- nmrks
		nxmrks <- round(nymrks / aspect)
	}
	# Construct a rectangular "raster" of equally spaced pixels/marks
	xmrks <- seq(from = xlim[1], to = xlim[2], length = nxmrks)
	ymrks <- seq(from = ylim[1], to = ylim[2], length = nymrks)
	marks <- cbind(rep(xmrks, times = nymrks), rep(ymrks, each = nxmrks))
	# Exclude pixels that fall outside the habitat outline
	outer_poly <-sp::SpatialPolygons(list(Polygons(list(Polygon(outer, hole = FALSE)), "1")))
	marks <- sp::SpatialPoints(marks)[outer_poly, ]
	marks <- marks@coords
	### set column names to x and y
	#outer <- dplyr::as_data_frame(outer) %>% setNames(c("x", "y"))
	#outer <- as.data.frame(outer) %>% setNames(c("x", "y"))
	outer <- tibble::as_tibble(outer) %>% setNames(c("x", "y"))
	### 
	#ipmap <- dplyr::data_frame(id = ipmap) %>% dplyr::count(id)
	#ipmap <- data.frame(id = ipmap) %>% dplyr::count(id)
	ipmap <- tibble::tibble(id = ipmap) %>% dplyr::count(id)
	### set column names to x and y
	#demes <- dplyr::as_data_frame(demes) %>% setNames(c("x", "y")) %>% dplyr::mutate(id = dplyr::row_number()) %>% dplyr::left_join(ipmap) %>% dplyr::arrange(id) %>% dplyr::mutate(n = dplyr::if_else(is.na(n), 0L, n))
	demes <- tibble::as_tibble(demes) %>% setNames(c("x", "y")) %>% dplyr::mutate(id = dplyr::row_number()) %>% dplyr::left_join(ipmap) %>% dplyr::arrange(id) %>% dplyr::mutate(n = dplyr::if_else(is.na(n), 0L, n))
	edges <- dplyr::bind_cols(demes[edges[, 1], ] %>% dplyr::select(x, y),demes[edges[, 2], ] %>% dplyr::select(x, y)) %>% setNames(c("x", "y", "xend", "yend"))
	list(marks = marks, nmrks = nrow(marks), xlim = xlim, ylim = ylim,outer = outer, demes = demes, edges = edges, dist_metric = dist_metric)
}

#' @title eems contours
#' 
#' Internal function used for making EEMS plots.
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @param dimns output of read_dimns function
#' @param longlat Whether or not to use long lat coordinates. Default TRUE
#' @param plot_params List of plotting parameters
#' @param is_mrates Logical indicating whether or not to generated the first two EEMs plots or the third and fourth EEMs plot
#' @return A list with parameters to use for generating EEMS plots.
#' @export eems_contours
eems_contours <- function(mcmcpath, dimns, longlat, plot_params, is_mrates) {
	if (is_mrates){
		message("Generate effective migration surface ","(posterior mean of m rates). ","See plots$mrates01 and plots$mrates02.")
	}
	else{
		message("Generate effective diversity surface ","(posterior mean of q rates). ","See plots$qrates01 and plots$qrates02.")
	}
	zrates <- rep(0, dimns$nmrks)
	pr_gt0 <- rep(0, dimns$nmrks)
	pr_lt0 <- rep(0, dimns$nmrks)
	niters <- 0
	# Loop over each directory in mcmcpath to average the contour plots
	for (path in mcmcpath) {
		voronoi <- read_voronoi(mcmcpath=path, longlat=longlat, is_mrates=is_mrates)
		# The function tiles2contours is coded in c++
		rslt    <- tiles2contours(tiles=voronoi$tiles, rates=voronoi$rates, seeds=cbind(voronoi$xseed, voronoi$yseed), marks=dimns$marks, distm=dimns$dist_metric)
		zrates  <- zrates + rslt$zrates
		niters  <- niters + rslt$niters
		pr_gt0  <- pr_gt0 + rslt$pr_gt0
		pr_lt0  <- pr_lt0 + rslt$pr_lt0
	}
	zrates <- zrates / niters
	pr_gt0 <- pr_gt0 / niters
	pr_lt0 <- pr_lt0 / niters
	p1     <- filled_eems_contour(dimns=dimns, zmean=zrates, plot_params=plot_params, is_mrates=is_mrates)
	p2     <- filled_prob_contour(dimns=dimns, probs=(pr_gt0 - pr_lt0), plot_params=plot_params, is_mrates=is_mrates)
	result=list(p1, p2)
	result
}

#' @title read vector
#' 
#' Internal function used for making EEMS plots. Reads eems output file as a vector.
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param file Path to file to read
#' @return Vector with file contents
#' @export read_vector
read_vector <- function(file) {
	stopifnot(file.exists(file))
	scan(file, what = numeric(), quiet = TRUE)
}

#' @title read matrix
#' 
#' Internal function used for making EEMS plots. Reads eems output file as a matrix.
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param file Path to file to read
#' @param ncol Number of columns in the file. Default 2.
#' @return matrix with the file contents
#' @export read_matrix
read_matrix <- function(file, ncol = 2) {
	stopifnot(file.exists(file))
	matrix(scan(file, what = numeric(), quiet = TRUE),ncol = ncol, byrow = TRUE)
}

#' @title pairwise distances between demes
#' 
#' Internal function used for making EEMS plots (plots 5–7). Calculates pairwise distances between demes.
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @param longlat Whether or not to use long lat coordinates. Default TRUE
#' @param plot_params List of plotting parameters
#' @return List with parameter values used for the function plot_pairwise_dissimilarities_
#' @export pairwise_dist
pairwise_dist <- function(mcmcpath, longlat, plot_params) {
	# List of observed demes, with number of samples taken collected Each row
	# specifies: x coordinate, y coordinate, n samples
	obs_demes <- read_matrix(file.path(mcmcpath[1], "rdistoDemes.txt"), ncol = 3)
	sizes <- obs_demes[, 3]
	if (sum(sizes > 1) < 2) {
		message("There should be at least two observed demes ","to plot pairwise dissimilarities")
		return(NULL)
	}
	npops <- nrow(obs_demes)
	demes <- seq(npops)
	diffs_obs <- matrix(0, npops, npops)
	diffs_hat <- matrix(0, npops, npops)
	for (path in mcmcpath) {
		tempi <- read_matrix(file.path(path, "rdistoDemes.txt"), ncol = 3)
		if (sum(dim(obs_demes) != dim(tempi)) || sum(obs_demes != tempi)) {
			message("EEMS results for at least two different population grids. ", "Plot pairwise dissimilarity for each grid separately.")
			#return(list(between = dplyr::data_frame(), within = dplyr::data_frame(), ibd = dplyr::data_frame()))
			return(list(between = tibble::tibble(), within = tibble::tibble(), ibd = tibble::tibble()))
		}
		diffs_obs <- diffs_obs + as.matrix(read.table(file.path(path, "rdistJtDobsJ.txt")))
		diffs_hat <- diffs_hat + as.matrix(read.table(file.path(path, "rdistJtDhatJ.txt")))
	}
	diffs_obs <- diffs_obs / length(mcmcpath)
	diffs_hat <- diffs_hat / length(mcmcpath)
	alpha <- matrix(demes, nrow = npops, ncol = npops)
	beta <- t(alpha)
	tempi <- matrix(sizes, npops, npops)
	smaller_deme <- pmin(tempi, t(tempi))
	smaller_deme <- smaller_deme[upper.tri(smaller_deme, diag = FALSE)]
	alpha <- alpha[upper.tri(alpha, diag = FALSE)]
	beta <- beta[upper.tri(beta, diag = FALSE)]
	# Under pure isolation by distance, we expect the genetic dissimilarities
	# between demes increase with the geographic distance separating them
	dist <- geo_distm(obs_demes[, 1:2], longlat, plot_params)
	bw_obs <- decompose_distances(diffs_obs, sizes)
	bw_hat <- decompose_distances(diffs_hat)
	#b_component <- dplyr::data_frame(alpha_x = obs_demes[, 1][alpha],alpha_y = obs_demes[, 2][alpha],beta_x = obs_demes[, 1][beta],beta_y = obs_demes[, 2][beta],fitted = bw_hat$between,obsrvd = bw_obs$between,size = smaller_deme)
	#w_component <- dplyr::data_frame(alpha_x = obs_demes[, 1][demes],alpha_y = obs_demes[, 2][demes],fitted = bw_hat$within,obsrvd = bw_obs$within,size = sizes)
	#g_component <- dplyr::data_frame(alpha_x = obs_demes[, 1][alpha],alpha_y = obs_demes[, 2][alpha],beta_x = obs_demes[, 1][beta],beta_y = obs_demes[, 2][beta],fitted = dist,obsrvd = bw_obs$between,size = smaller_deme)
	b_component <- tibble::tibble(alpha_x = obs_demes[, 1][alpha],alpha_y = obs_demes[, 2][alpha],beta_x = obs_demes[, 1][beta],beta_y = obs_demes[, 2][beta],fitted = bw_hat$between,obsrvd = bw_obs$between,size = smaller_deme)
	w_component <- tibble::tibble(alpha_x = obs_demes[, 1][demes],alpha_y = obs_demes[, 2][demes],fitted = bw_hat$within,obsrvd = bw_obs$within,size = sizes)
	g_component <- tibble::tibble(alpha_x = obs_demes[, 1][alpha],alpha_y = obs_demes[, 2][alpha],beta_x = obs_demes[, 1][beta],beta_y = obs_demes[, 2][beta],fitted = dist,obsrvd = bw_obs$between,size = smaller_deme)


	list(between = b_component, within = w_component, ibd = g_component)
}


#' @title plot pairwise dissimilarities between demes
#' 
#' Internal function used for making EEMS plots (plots 5–7).
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param dissimilarities Output list from the function pairwise_dist
#' @param add_abline Whether or not to add the abline.
#' @return List of ggplot objects
#' @export plot_pairwise_dissimilarities_
plot_pairwise_dissimilarities_ <- function(dissimilarities, add_abline) {
	message("Generate average dissimilarities within and between demes. ","See plots$rdist01, plots$rdist02 and plots$rdist03.")
	p1 <- ggplot2::ggplot(dissimilarities$between %>% dplyr::filter(size > 1),ggplot2::aes(fitted, obsrvd)) + ggplot2::geom_point(shape = 1) + ggplot2::theme_minimal() + ggplot2::labs(x = expression(paste("Fitted dissimilarity between demes  ", Delta[alpha * beta], " - (", Delta[alpha * alpha], " + ", Delta[beta * beta], ") / 2")), y = expression(paste("Observed dissimilarity between demes  ",D[alpha * beta], " - (",D[alpha * alpha], " + ",D[beta * beta], ") / 2")), title = expression(paste("Dissimilarities between pairs of ", "sampled demes (", alpha, ", ", beta, ")")),subtitle = paste("Singleton demes, if any, are excluded from this", "plot (but not from EEMS)"))
	p2 <- ggplot2::ggplot(dissimilarities$within %>% dplyr::filter(size > 1), ggplot2::aes(fitted, obsrvd)) + ggplot2::geom_point(shape = 1) + ggplot2::theme_minimal() + ggplot2::labs(x = expression(paste("Fitted dissimilarity within demes  ", Delta[alpha * alpha])), y = expression(paste("Observed dissimilarity within demes ",D[alpha * alpha])), title = expression(paste("Dissimilarities within sampled ","demes ", alpha)), subtitle = paste("Singleton demes, if any, are excluded from ","this plot (but not from EEMS)"))
	p3 <- ggplot2::ggplot(dissimilarities$ibd %>% dplyr::filter(size > 1), ggplot2::aes(fitted, obsrvd)) + ggplot2::geom_point(shape = 1) + ggplot2::theme_minimal() + ggplot2::labs(x = "Great circle distance between demes (km)",y = expression(paste("Observed dissimilarity between demes  ",D[alpha * beta], " - (",D[alpha * alpha], " + ",D[beta * beta], ") / 2")),title = expression(paste("Dissimilarities between pairs of ","sampled demes (", alpha, ", ", beta, ")")),subtitle = paste("Singleton demes, if any, are excluded from this","plot (but not from EEMS)"))
	if (add_abline) {
		p1 <- p1 + ggplot2::geom_smooth(method = "lm", se = FALSE)
		p2 <- p2 + ggplot2::geom_smooth(method = "lm", se = FALSE)
	}
	result <- list(p1, p2, p3)
	result
}

#' @title plot eems log posterior
#' 
#' Used for making EEMS plot 8
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @return Object of class ggplot
#' @export plot_log_posterior
plot_log_posterior <- function(mcmcpath) {
	message("Generate posterior probability trace. ","See plots$pilog01.")
	rleid <- function(x) {
		r   <- rle(x)
		rep(seq_along(r$lengths), r$lengths)
	}
	pl_df <- NULL
	for (path in mcmcpath) {
		pl    <- read_matrix(file.path(path, "mcmcpilogl.txt"))
		#pl_df <- dplyr::bind_rows(pl_df, dplyr::as_data_frame(pl) %>% dplyr::mutate(path))
		pl_df <- dplyr::bind_rows(pl_df, tibble::as_tibble(pl) %>% dplyr::mutate(path))

	}
	pl_df <- pl_df %>% setNames(c("pi", "logl", "path")) %>% dplyr::mutate(mcmcpath = factor(data.table::rleid(path))) %>% dplyr::group_by(mcmcpath) %>% dplyr::mutate(iter = dplyr::row_number(), pilogl = pi + logl)
		ggplot2::ggplot(pl_df, ggplot2::aes(x = iter, y = pilogl, color = mcmcpath)) + ggplot2::geom_path() + ggplot2::labs(x = "MCMC iteration  (after burn-in and thinning)",y = "log posterior",title = "Have the MCMC chains converged?",subtitle = "If not, restart EEMS and/or increase numMCMCIter") + ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),panel.grid.major.x = ggplot2::element_blank())
}

#' @title geo distances
#' 
#' Used for making EEMS plots
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param coord input coordinates
#' @param longlat Whether or not to use long lat coordinates. Default TRUE
#' @param plot_params List of plotting parameters
#' @return Distance between two points
#' @export geo_distm
geo_distm <- function(coord, longlat, plot_params) {
	if (!longlat) coord <- coord[, c(2, 1)]
	dist <- sp::spDists(coord, longlat = TRUE)
	dist <- dist[upper.tri(dist, diag = FALSE)]
	dist
}

#' @title geo distances
#' 
#' Used for making EEMS plots
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param diffs Differences
#' @param sizes Sizes
#' @return List with parameters
#' @export decompose_distances
decompose_distances <- function(diffs, sizes = NULL) {
	# Diffs can have NAs on the main diagonal; these elements correspond to demes
	# with a single observation. For such deme a, no dissimilarities between
	# two distinct individuals are observed. I approximate diffs(a,a) with the
	# average diffs(b,b) computed across demes b with multiple samples.
	if (!is.null(sizes))
	  diag(diffs)[sizes < 2] <- mean(diag(diffs)[sizes >= 2])
	within <- diag(diffs)
	selfsim <- matrix(within, nrow(diffs), ncol(diffs))
	between <- diffs - (selfsim + t(selfsim)) / 2
	between <- between[upper.tri(between, diag = FALSE)]
	list(within = within, between = between)
}

#' @title load required packages
#' 
#' Load packages
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param packages Character string with names of packages to load.
#' @return message stating packages loaded
#' @export load_required_packages
load_required_packages <- function(packages) {
	for (package in packages) {
		if (!requireNamespace(package, quietly = TRUE)) {
			stop("The ", package, " package is required. ","Please install it first.")
		} else {
			message("Loading ", package, ".")
		}
	}
}


#' @title tiles2contours
#' 
#' Converts tiles two contours
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param tiles voronoi tiles
#' @param rates voronoi rates
#' @param seeds Two column matrix with voronoi xseeds and voronoi yseeds
#' @param marks 'marks' named entry of output list of read_dimns function
#' @param distm 'dist_metric' named entry of output list of read_dimns function
#' @return countours
#' @export tiles2contours
tiles2contours <- function(tiles, rates, seeds, marks, distm) {
	require("reemsplots2")
	.Call('_reemsplots2_tiles2contours', PACKAGE = 'reemsplots2', tiles, rates, seeds, marks, distm)
}

#' @title read voronoi
#' 
#' Reads EEMs output and constructs parameters for voronoi interpolation
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @param longlat Whether or not to use long lat coordinates. Default TRUE
#' @param is_mrates Logical indicating if the plots are mrates, else qrates
#' @return List with parameter values used for veronoi interpolation
#' @export read_voronoi
read_voronoi <- function(mcmcpath, longlat, is_mrates) {
	if (is_mrates) {
		rates <- read_vector(file.path(mcmcpath, "mcmcmrates.txt"))
		tiles <- read_vector(file.path(mcmcpath, "mcmcmtiles.txt"))
		xseed <- read_vector(file.path(mcmcpath, "mcmcxcoord.txt"))
		yseed <- read_vector(file.path(mcmcpath, "mcmcycoord.txt"))
	} else {
		rates <- read_vector(file.path(mcmcpath, "mcmcqrates.txt"))
		tiles <- read_vector(file.path(mcmcpath, "mcmcqtiles.txt"))
		xseed <- read_vector(file.path(mcmcpath, "mcmcwcoord.txt"))
		yseed <- read_vector(file.path(mcmcpath, "mcmczcoord.txt"))
	}
	if (!longlat) {
		tempi <- xseed
		xseed <- yseed
		yseed <- tempi
	}
	list(rates = log10(rates), tiles = tiles, xseed = xseed, yseed = yseed)
}

#' @title get distance metric
#' 
#' Get the distance metric from `eemsrun.txt`; use `euclidean` by default.
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param mcmcpath Path the directory containing EEMs output files for a particular chain
#' @return Character string with distance metric to use
#' @export get_dist_metric
get_dist_metric <- function(mcmcpath) {
	dist_metric <- "euclidean"
	lines <- tolower(readLines(file.path(mcmcpath, "eemsrun.txt")))
	for (line in lines) {
		if (grepl("\\s*distance\\s*=\\s*", line)){
			dist_metric <- gsub("\\s*distance\\s*=\\s*(\\w+)", "\\1", line)
		}
	}
	if (dist_metric != "euclidean" && dist_metric != "greatcirc"){
		stop("eemsrun.txt should specify `euclidean` or `greatcirc` distance.")
	}
	dist_metric
}

#' @title default eems colors
#' 
#' Use the default DarkOrange to Blue color scheme, which combines two color schemes from the `dichromat` package. These are based
#' on a collection of color schemes for scientific graphics:
#' See http://geog.uoregon.edu/datagraphics/color_scales.htm
#' To reproduce the default eems colors, let oranges be dichromat::colorschemes$BluetoDarkOrange.12[12:7]
#' and blues be dichromat::colorschemes$BrowntoBlue.12[7:12]
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @return Character string vector with hex codes for default eems colors
#' @export default_eems_colors
default_eems_colors <- function() {
	c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC","#FBFBFB","#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99")
}

#' @title Is color
#' 
#' Check that a string, or a vector of strings, represents a hex color
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param x Character string to check
#' @return Logical
#' @export is_color
is_color <- function(x) {
	if (is.null(x)) {
		return(FALSE)
	} else {
		return(sapply(x, function(x) {tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)}))
	}
}

#' @title theme void
#' 
#' mostly blank ggplot theme
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @return ggplot object
#' @export theme_void
theme_void <- function() {
	ggplot2::theme(line = ggplot2::element_blank(), rect = ggplot2::element_blank(),axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),legend.title = ggplot2::element_text(hjust = 0), legend.text.align = 1)
}

#' @title filled contour rates
#' 
#' Generates a ggplot object with filled contours
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param z vector of zeros with length determined by the named entry 'nmrks' of the output of read_dimns function
#' @param dimns output of read_dimns function
#' @return ggplot object
#' @export filled_contour_rates
filled_contour_rates <- function(z, dimns) {
	w <- cbind(dimns$marks, z)
	colnames(w) <- c("x", "y", "z")
	#ggplot2::ggplot(dplyr::as_data_frame(w), ggplot2::aes(x = x, y = y)) + ggplot2::geom_tile(ggplot2::aes(fill = z)) + ggplot2::scale_x_continuous(limits = dimns$xlim) + ggplot2::scale_y_continuous(limits = dimns$ylim) + ggplot2::coord_quickmap() + theme_void() + ggplot2::theme(legend.text.align = 1)
	ggplot2::ggplot(as.data.frame(w), ggplot2::aes(x = x, y = y)) + ggplot2::geom_tile(ggplot2::aes(fill = z)) + ggplot2::scale_x_continuous(limits = dimns$xlim) + ggplot2::scale_y_continuous(limits = dimns$ylim) + ggplot2::coord_quickmap() + theme_void() + ggplot2::theme(legend.text.align = 1)
}

#' @title filled contour graph
#' 
#' Generates a filled countour graph ggplot object
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param p output of the function filled_contour_rates
#' @param dimns output of read_dimns function
#' @param plot_params List of plotting parameters
#' @return ggplot object
#' @export filled_contour_graph
filled_contour_graph <- function(p, dimns, plot_params) {
	if (plot_params$add_grid) {
		p <- p + ggplot2::geom_segment(data = dimns$edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), color = plot_params$col_grid)
	}
	if (plot_params$add_demes) {
		p <- p + ggplot2::geom_point(data = dimns$demes %>% dplyr::filter(n > 0), ggplot2::aes(x = x, y = y, size = n), shape = 1,color = plot_params$col_demes) + ggplot2::scale_size_continuous(guide = FALSE)
	}
	if (plot_params$add_outline) {
		p <- p + ggplot2::geom_path(data = dimns$outer, ggplot2::aes(x = x, y = y), color = plot_params$col_outline)
	}
	p
}

#' @title filled eems contour
#' 
#' generates ggplot object corresponding to eems plots 1 and 3
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param dimns output of read_dimns function
#' @param zmean vector of zeros with length determined by the named entry 'nmrks' of the output of read_dimns function
#' @param plot_params List of plotting parameters
#' @param is_mrates Logical indicating if the plots are mrates, else qrates
#' @return ggplot object
#' @export filled_eems_contour
filled_eems_contour <- function(dimns, zmean, plot_params, is_mrates) {
	if (is_mrates) {
		title <- "log(m)"
		colscale <- plot_params$m_colscale
	} else {
		title <- "log(q)"
		colscale <- plot_params$q_colscale
	}
	limits <- range(zmean, colscale, na.rm = TRUE, finite = TRUE)
	p <- filled_contour_rates(zmean, dimns)
	p <- filled_contour_graph(p, dimns, plot_params) + ggplot2::scale_fill_gradientn(colors = plot_params$eems_colors, limits = limits, name = title)
	p
}

#' @title filled probability contour
#' 
#' Generates asthetics for ggplots.
#' This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
#' 
#' @param dimns output of read_dimns function
#' @param probs Probabilities
#' @param plot_params List of plot parameters
#' @param is_mrates Logical indicating if the plots are mrates, else qrates
#' @return ggplot aesthetics
#' @export filled_prob_contour
filled_prob_contour <- function(dimns, probs, plot_params, is_mrates) {
	probs <- (probs + 1) / 2
	probs[probs < 0] <- 0
	probs[probs > 1] <- 1
	if (is_mrates) r <- "m" else r <- "q"
	breaks <- c(1 - plot_params$prob_level, plot_params$prob_level)
	labels <- c(paste0("P{log(", r, ") < 0} = ", plot_params$prob_level),paste0("P{log(", r, ") > 0} = ", plot_params$prob_level))
	p <- filled_contour_rates(probs, dimns)
	p <- filled_contour_graph(p, dimns, plot_params) + ggplot2::geom_contour(ggplot2::aes(z = z), breaks = breaks, color = "white") + ggplot2::scale_fill_gradientn(colors = default_eems_colors(), limits = c(0, 1), name = "", breaks = breaks, labels = labels)
	p
}






