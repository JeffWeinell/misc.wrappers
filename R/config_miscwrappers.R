#' @title config_miscwrappers function
#' 
#' 
#' This function is used to tell misc.wrappers where a program (executable file) is located. You only need to run this function once for a particular program.
#' 
#' @param exe.paths Character string or vector with the full path(s) to one or more executable files for programs that misc.wrappers uses.
#' @return If exe.paths is non-NULL, the function will update settings and return 0 if successful. If exe.paths is NULL, the function returns a dataframe with all of the path settings.
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
		result <- 0
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
#' @examples
#' library(misc.wrappers)
#' 
#' # Example 1: view current path settings
#' config_miscwrappers()
#' # Example 2: set path to program 'runeems_snps' (replace example path with the real path on you machine)
#' config_miscwrappers(exe.path="/FULL/PATH/TO/runeems_snps")
#' # Example 3: set path to the fastStructure file 'structure.py'
#' config_miscwrappers(exe.path="/FULL/PATH/TO/structure.py")

#' @title sh single quote
#' 
#' Calls the base function shQuote with default options; sq is faster to type.
#' 
#' @param x Character string to be quoted
#' @return Character string with leading and trailing "'" character added to input string
#' @export sq
sq <- function(x){
	shQuote(x)
}

#' @title update miscwrappers
#' 
#' Updates misc.wrappers by calling remotes::install_github function
#' 
#' @param upgrade Whether or not to update dependencies. Default FALSE.
#' @param ... Additional arguments to pass to install_github
#' @return NULL
#' @export update_miscwrappers
update_miscwrappers <- function(upgrade=FALSE,...){
	remotes::install_github("JeffWeinell/misc.wrappers",upgrade=upgrade,...)
}

#' @title sbatch setup and submit
#' 
#' Creates slurm sbatch submission string, and optionally submits the job using system
#' 
#' @param sh.path Character string with path to bash script
#' @param partition Character string with name of partition to submit the job to. Passed to sbatch argument '--partition'. Default NA.
#' @param nodes Number of nodes to split the job across. Default 1.
#' @param ntasksPerNode Number of tasks (cores) to use per node. Default 1.
#' @param memGb Total amount of memory to allocate for the job. Default 50Gb.
#' @param time Numerical vector with length 4, with the number of days, hours, minutes, and seconds, respectively, to allocate for the job.
#' @param submit Logical with whether or not to submit the jobs to the cluster.
#' @return Either vector of character strings that can be piped to system, or NULL
#' @export rsbatch
rsbatch <- function(sh.path, partition=NA, nodes=1, ntasksPerNode=1, memGb=50, time=c(0,6,0,0), submit=FALSE, intern=FALSE){
	#sprintf("%s-%s:%s:%s",list(0,6,0,0))
	#do.call(sprintf,list(0,6,0,0),"%s-%s:%s:%s")
	#do.call(sprintf,c("%s-%s:%s:%s",list(0,6,0,0)))
	time2    <- lapply(time,FUN=formatC,width=2,format="d",flag="0")
	#pars.df  <- data.frame(args=c("--partition=","--nodes=","--ntasks-per-node=","--mem=","--time="),values=c(partition,nodes,ntasksPerNode,paste0(memGb,"Gb"),paste0(days,"-",hrs,":00:00")))
	pars.df  <- data.frame(args=c("--partition=","--nodes=","--ntasks-per-node=","--mem=","--time="), values=c(partition,nodes,ntasksPerNode,paste0(memGb,"Gb"),do.call(sprintf,c("%s-%s:%s:%s",time2))))
	pars.use <- !is.na(pars.df[,"values"])
	pars.df2 <- pars.df[pars.use,]
	params   <- paste(sapply(1:nrow(pars.df2),FUN=function(x) {paste0(pars.df2[x,],collapse="")}),collapse=" ")
	# params <- sprintf("sbatch --nodes=%s --ntasks-per-node=%s --mem=%sGb --time=%s-%s:00:00 --partition=%s %s",nodes,ntasksPerNode,memGb,days,hrs,partition,sq(sh.path))
	jobstrings <- paste("sbatch",params,sq(sh.path))
	if(!submit){
		return(jobstrings)
	} else {
		lapply(jobstrings, system, intern=intern)
	}
}

#' @title write bash and Rscript files for R code
#' 
#' Creates .sh (bash) and .R (Rscript) files given lines of R code held as a character string or character string vector, and optional arguments indicating which modules and packages to load.
#' This function probably will not work on all systems.
#' 
#' @param name Character string with name of the job. Bash and Rscript files will have names 'name.sh' and 'name.R'
#' @param Rcode Character string, possibly a vector, with lines of R code. Line breaks can be coded as \n.
#' @param modules Character string vector with names of modules to load. The defaults are the current versions of programs needed to run misc.wrappers package on my system, and therefore these defaults will probably need to be changed for other systems.
#' @param libPath Character string with path where R package libraries are held. Default is .libPaths()[1]
#' @param libs Character string vector with names of R packages to load. Default is 'misc.wrappers'.
#' @param dir Character string with path to the directory where bash and Rscript files should be written. Default is the current directory.
#' @return NULL
#' @export Rsh
Rsh <- function(name="job", Rcode, modules=c("compiler/gcc/8.3","gdal/3.0.0","geos/3.7.2","R/4.0"), libPath=.libPaths()[1], libs=c("misc.wrappers"), dir=getwd()){
	shlines     <- c("#!/bin/bash",sprintf("module load %s",modules),sprintf("Rscript %s/%s.R",dir,name),"")
	Rcodelines  <- unlist(strsplit(Rcode,split="\n"))
	Rlines      <- c(sprintf(".libPaths('%s')",libPath),sprintf("library('%s')",libs),Rcodelines,"")
	writeLines(shlines,sprintf("%s/%s.sh",dir,name))
	writeLines(Rlines,sprintf("%s/%s.R",dir,name))
}



