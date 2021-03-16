#' @title config_miscwrappers function
#' 
#' This function is used to tell misc.wrappers where an EEMS executable is located. You only need to run this function once for a particular EEMS program.
#' 
#' @param exe.paths Character string or vector with the full path(s) to one or more EEMS executables.
##' @param program A character string or character vector containing the names of the program(s) with their paths defined by ```eems_exe``` argument. Default "runeems_snps"; other options could include "str2diffs", "runeems_sats", "runeems_Fsts", "bed2diffs", "rEEMSplots" or a character vector including all or a subset of these. The programs must be executable first.
#' @return Updates settings and returns NULL if exe.paths is a valid character strings or vector of path(s); returns a dataframe with the all path settings if exe.paths is NULL.
#' @export config_miscwrappers
config_miscwrappers <- function(exe.paths=NULL){
	settings.dir  <- path.expand(rappdirs::user_config_dir("R_misc_wrappers", version=packageVersion("misc.wrappers")))
	settings.file <- paste0(settings.dir,"/exepaths.txt")
	if(is(exe.paths,"character")){
		# runeems_snps.exe.path <- paste0(eems_master.dirpath,"/runeems_snps/src/runeems_snps")
		files.executable <- sapply(exe.paths,FUN=check.if.executable)
		if(!all(files.executable==0)){
			stop(paste0("not executable: ",exe.paths[which(files.executable!=0)]))
		}
		if(!dir.exists(settings.dir)){
			dir.create(settings.dir,recursive=T)
		}
		settings.mat  <- cbind(basename(exe.paths),exe.paths)
		colnames(settings.mat) <- c("program","exe_path")
		if(!file.exists(settings.file)){
			write.table(settings.mat,settings.file,col.names=T,row.names=F,quote=T)
		} else {
			current.settings.mat  <- read.table(settings.file,sep=" ",header=T)
			current.settings.mat  <- cbind(current.settings.mat,path.version=rep(0,nrow(current.settings.mat)))
			add.settings.mat      <- cbind(settings.mat,path.version=rep(1,nrow(settings.mat)))
			combined.settings.mat <- rbind(current.settings.mat,add.settings.mat)
			combined.settings.mat.order <- combined.settings.mat[order(combined.settings.mat[,"program"],combined.settings.mat[,"path.version"],decreasing=T),]
			### Gets first row (most current) of each program setting
			use.settings.mat <- combined.settings.mat.order[match(unique(combined.settings.mat.order[,"program"]),combined.settings.mat.order[,"program"]),c(1:2)]
			write.table(use.settings.mat,settings.file,col.names=T,row.names=F,quote=T)
		}
		result <- NULL
	} else {
		if(is.null(exe.paths)){
			if(!dir.exists(settings.dir)){
				dir.create(settings.dir,recursive=T)
			}
			### If the settings file doesnt exist, create a settings file containing column headers but without any other rows.
			if(!file.exists(settings.file)){
				writeLines(text="program exe_path",con=settings.file)
			}
			current.settings.mat  <- read.table(settings.file,sep=" ",header=T)
		}
		result <- current.settings.mat
	}
	result
}

