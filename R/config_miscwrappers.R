#' @title config_miscwrappers function
#' 
#' This function is used to tell misc.wrappers where an EEMS executable is located. You only need to run this function once for a particular EEMS program.
#' 
#' @param exe.paths Character string or vector with the full path(s) to one or more EEMS executables.
##' @param program A character string or character vector containing the names of the program(s) with their paths defined by ```eems_exe``` argument. Default "runeems_snps"; other options could include "str2diffs", "runeems_sats", "runeems_Fsts", "bed2diffs", "rEEMSplots" or a character vector including all or a subset of these. The programs must be executable first.
#' @return NULL
#' @export config_miscwrappers
config_miscwrappers <- function(exe.paths){
	# runeems_snps.exe.path <- paste0(eems_master.dirpath,"/runeems_snps/src/runeems_snps")
	files.executable <- sapply(exe.paths,FUN=check.if.executable)
	if(!all(files.executable==0)){
		stop(paste0("not executable: ",exe.paths[which(files.executable!=0)]))
	}
	settings.dir  <- path.expand(rappdirs::user_config_dir("R_misc_wrappers", version=packageVersion("misc.wrappers")))
	dir.create(settings.dir,recursive=T)
	settings.mat  <- cbind(basename(exe.paths),exe.paths)
	colnames(settings.mat) <- c("program","exe_path")
	settings.file <- paste0(settings.dir,"/exepaths.txt")
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
}

## eems_master.dirpath <- "/Users/alyssaleinweber/Documents/eems-master"
## config_miscwrappers("/panfs/pfs.local/home/j926w878/programs/eems-master/runeems_snps/src/runeems_snps")
## or config_miscwrappers("/panfs/pfs.local/home/j926w878/programs/eems-master/")
