#' @title exepaths_miscwrappers function
#' 
#' Function to return the paths of executables that were previously defined using the function congig_miscwrappers.
#' 
#' @param programs Character string or vector with names of programs to return path settings for, or NULL, in which case all defined path settings are returned.
#' @return Two column data.frame with program name and corresponding full path to the program executable file
#' @export exepaths_miscwrappers
exepaths_miscwrappers <- function(programs=NULL){
	settings.file <- paste0(path.expand(user_config_dir("R_misc_wrappers", version=packageVersion("misc.wrappers"))),"/exepaths.txt")
	if(!file.exists(settings.file)){
		stop("No paths set for executables. Use config_miscwrappers function to set paths")
	}
	settings.mat  <- read.table(settings.file,sep=" ",header=T)
	if(!is.null(programs)){
		if(any(settings.mat[,"program"] %in% programs)){
			paths.defined.mat    <- settings.mat[(settings.mat[,"program"] %in% programs),]
			paths.defined        <- paths.defined.mat[,"exe_path"]
			names(paths.defined) <- paths.defined.mat[,"program"]
		} else {
			paths.defined <- NULL
		}
		paths.undefined   <- setdiff(programs,settings.mat[,"program"])
		if(length(paths.undefined)!=0){
			NAstring <- rep(NA,length(paths.undefined))
			names(NAstring) <- paths.undefined
		} else {
			NAstring <- NULL
		}
		result <- c(paths.defined,NAstring)
	} else {
		result <- settings.mat[,"exe_path"]
		names(result) <- settings.mat[,"program"]
	}
	result
}
